"""Numerical helpers shared by the decomposition methods."""

from __future__ import annotations

from collections.abc import Sequence
from math import ceil, floor
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import windows


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


def _window_values(length: int, window_type: str, parameter: float | None) -> FloatArray:
    name = window_type.lower()
    if name == "delta":
        result = np.zeros(length, dtype=np.float64)
        result[length // 2] = 1.0
    elif name in {"rect", "boxcar"}:
        result = np.ones(length, dtype=np.float64)
    elif name in {"bart", "bartlett"}:
        result = windows.bartlett(length, sym=True)
    elif name in {"hamm", "hamming"}:
        result = windows.hamming(length, sym=True)
    elif name in {"hann", "hanning"}:
        result = windows.hann(length, sym=True)
    elif name == "tukey":
        result = windows.tukey(length, alpha=0.5 if parameter is None else parameter, sym=True)
    elif name in {"gauss", "gaussian"}:
        alpha = 2.5 if parameter is None else parameter
        standard_deviation = (length - 1) / (2 * alpha) if length > 1 else 1.0
        result = windows.gaussian(length, std=standard_deviation, sym=True)
    elif name == "cosh":
        coefficient = 0.01 if parameter is None else parameter
        half = int(np.fix(length / 2))
        result = np.zeros(length, dtype=np.float64)
        for m in range(-half, half + 1):
            result[m % length] = np.cosh(m) ** (-2 * coefficient)
        result = np.fft.fftshift(result)
    elif name in {"blackmanharris", "bmh"}:
        result = windows.blackmanharris(length, sym=True)
    elif name == "nuttall":
        result = windows.nuttall(length, sym=True)
    elif name in {"dolph", "chebwin"}:
        if parameter is None:
            raise ValueError("the Dolph-Chebyshev window requires an attenuation parameter")
        result = windows.chebwin(length, at=parameter, sym=True)
    else:
        raise ValueError(f"unknown window type: {window_type}")
    return np.asarray(result, dtype=np.float64)


def get_window(
    length: int,
    window_type: str,
    parameter: float | None = None,
    dft_window: bool = False,
    pad_to: int | None = None,
) -> NDArray[np.float64] | NDArray[np.complex128]:
    """Generate and shift a window using the conventions in ``get_window.m``."""
    if length < 1:
        raise ValueError("window length must be positive")
    result: NDArray[np.float64] | NDArray[np.complex128]
    result = _window_values(length, window_type, parameter)

    if dft_window:
        result = np.roll(result, ceil(length / 2))
        result = np.fft.fft(result)
        result = np.roll(result, floor(length / 2))

    # Put the samples for non-negative indices first.
    result = np.roll(result, ceil(length / 2))
    if pad_to is not None and pad_to > 0:
        result = pad_window(result, pad_to)
    return result


def pad_window(window: ArrayLike, length: int) -> NDArray[np.generic]:
    """Zero-pad a window whose non-negative indices are stored first."""
    source = np.asarray(window).reshape(-1)
    source_length = source.size
    if length < source_length:
        raise ValueError("padded length cannot be less than the window length")
    if length == source_length:
        return source.copy()

    result = np.zeros(length, dtype=source.dtype)
    half = source_length // 2
    if source_length % 2:
        result[: half + 1] = source[: half + 1]
        if half:
            result[-half:] = source[-half:]
    else:
        result[:half] = source[:half]
        result[half] = source[half] / 2
        if half > 1:
            result[-(half - 1) :] = source[-(half - 1) :]
        result[-half] = source[half] / 2
    return result


def parse_window_spec(spec: Sequence[object], default_dft: bool) -> tuple[int, str, object, bool]:
    """Read MATLAB's ``{length, type, parameter, domain_flag}`` representation."""
    if len(spec) < 2:
        raise ValueError("a window specification needs at least a length and type")
    length = int(spec[0])
    name = str(spec[1])
    parameter = spec[2] if len(spec) >= 3 else None
    dft = bool(spec[3]) if len(spec) >= 4 else default_dft
    return length, name, parameter, dft


def get_analytic_signal(
    signal: ArrayLike,
) -> tuple[ComplexArray, int, int, int]:
    """Create the zero-padded analytic signal used by the quadratic TFD code."""
    z = np.asarray(signal).reshape(-1)
    if z.size % 2:
        warnings.warn(
            "odd-length signal: dropping the final sample before analytic conversion",
            RuntimeWarning,
            stacklevel=2,
        )
        z = z[:-1]

    if np.isrealobj(z):
        source = np.asarray(z, dtype=np.float64)
        n = source.size
        n2 = 2 * n
        padded = np.concatenate((source, np.zeros(n, dtype=np.float64)))
        spectrum = np.fft.fft(padded)
        multiplier = np.zeros(n2, dtype=np.float64)
        multiplier[[0, n]] = 1.0
        multiplier[1:n] = 2.0
        circular = np.fft.ifft(spectrum * multiplier)
        analytic = np.concatenate((circular[:n], np.zeros(n, dtype=np.complex128)))
    else:
        analytic = np.asarray(z, dtype=np.complex128)
        n2 = analytic.size
        n = n2 // 2
        if n2 != 2 * n:
            raise ValueError("an analytic signal must have length 2N")

    return analytic, n2, n, ceil(n / 2)


def scale_tfd(tfd: ArrayLike, lag_window: ArrayLike, signal_length: int) -> FloatArray:
    """Scale a TFD to satisfy the MATLAB implementation's energy convention."""
    values = np.asarray(tfd, dtype=np.float64)
    lag = np.asarray(lag_window)
    ntime, nfreq = values.shape
    scale = (nfreq / ntime) * signal_length / np.sum(lag)
    return np.asarray(values * np.real_if_close(scale), dtype=np.float64)


def peak_pick(values: ArrayLike) -> NDArray[np.bool_]:
    """Locate peaks with the nested-sign difference rule used by MATLAB."""
    x = np.asarray(values).reshape(-1)
    result = np.zeros(x.size, dtype=bool)
    if x.size >= 3:
        inner = np.diff(np.sign(np.diff(x))) + 0.5
        result[1:-1] = np.sign(-np.sign(inner) + 1).astype(bool)
    return result


def get_bandwidth(
    values: ArrayLike, peak_index: int | None = None, amplitude_fraction: float = 0.5
) -> float | None:
    """Estimate peak bandwidth with linear crossing interpolation."""
    x = np.asarray(values, dtype=np.float64).reshape(-1)
    if x.size == 0:
        return None
    peak = int(np.argmax(x)) if peak_index is None else int(peak_index)
    if not 0 <= peak < x.size:
        raise IndexError("peak_index is outside the time-frequency slice")
    target = x[peak] * amplitude_fraction

    left_candidates = np.flatnonzero(x[: peak + 1] <= target)
    left = float(left_candidates[-1]) if left_candidates.size else 0.0
    right_candidates = np.flatnonzero(x[peak + 1 :] <= target)
    right = float(peak + 1 + right_candidates[0]) if right_candidates.size else float(x.size - 1)

    left_i = int(left)
    if x[left_i] < target and left_i + 1 < x.size:
        slope = x[left_i + 1] - x[left_i]
        if slope != 0:
            left = left_i + (target - x[left_i]) / slope

    right_i = int(right)
    if x[right_i] < target and right_i > 0:
        slope = x[right_i] - x[right_i - 1]
        if slope != 0:
            right = right_i + (target - x[right_i]) / slope

    bandwidth = right - left
    return bandwidth if bandwidth >= 0 else None


def preprocess_signal(signal: ArrayLike, pad_signal: bool = False) -> FloatArray:
    """Apply the positive, windowed mirror padding used by MATLAB."""
    x = np.asarray(signal, dtype=np.float64).reshape(-1)
    if not pad_signal:
        return x.copy()
    length = x.size // 16
    if length == 0:
        return x.copy()
    taper = windows.hamming(2 * length - 1, sym=True)[:length]
    left = taper * x[1 : length + 1][::-1]
    right = (taper * x[-length - 1 : -1])[::-1]
    return np.concatenate((left, x, right))
