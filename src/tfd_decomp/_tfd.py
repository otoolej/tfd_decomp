"""Quadratic and cross time-frequency distributions."""

from __future__ import annotations

from collections.abc import Sequence
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._utils import get_analytic_signal, get_window, pad_window, parse_window_spec


def generate_doppler_kernel(
    specification: Sequence[object], signal_length: int, ntime: int | None = None
) -> tuple[NDArray[np.generic], int, int]:
    """Generate the separable kernel's Doppler window."""
    length, name, parameter, in_doppler_domain = parse_window_spec(
        specification, default_dft=True
    )
    if length > signal_length:
        warnings.warn(
            "Doppler-kernel length exceeds the signal; reducing it",
            RuntimeWarning,
            stacklevel=2,
        )
        length = signal_length
    if length % 2 == 0:
        length -= 1
        warnings.warn(
            f"forcing Doppler-kernel length to the odd value {length}",
            RuntimeWarning,
            stacklevel=2,
        )
    if length < 1:
        raise ValueError("Doppler-kernel length is too short")

    if ntime is None:
        ntime = length + 1
    if ntime < length + 1:
        warnings.warn(
            "ntime is too short; increasing it to Doppler-kernel length + 1",
            RuntimeWarning,
            stacklevel=2,
        )
        ntime = length + 1
    if ntime % 2:
        ntime += 1
        warnings.warn(
            f"forcing ntime to the even value {ntime}", RuntimeWarning, stacklevel=2
        )

    if in_doppler_domain:
        kernel = get_window(length, name, parameter, dft_window=False)
    else:
        time_window = get_window(length, name, parameter, dft_window=False)
        padded = pad_window(time_window, signal_length)
        kernel = np.real(np.fft.fft(padded))
        kernel = kernel / kernel[0]
        ntime = signal_length
        length = signal_length

    if len(specification) > 4 and str(specification[4]).lower() == "y":
        kernel = kernel / np.sum(kernel)
    return np.asarray(kernel), length, ntime


def generate_lag_kernel(
    specification: Sequence[object], signal_length: int, nfreq: int | None = None
) -> tuple[NDArray[np.generic], int, int]:
    """Generate the separable kernel's lag window."""
    length, name, parameter, in_frequency_domain = parse_window_spec(
        specification, default_dft=False
    )
    if length > signal_length:
        warnings.warn(
            "lag-kernel length exceeds the signal; reducing it",
            RuntimeWarning,
            stacklevel=2,
        )
        length = signal_length
    if length % 2 == 0:
        length -= 1
        warnings.warn(
            f"forcing lag-kernel length to the odd value {length}",
            RuntimeWarning,
            stacklevel=2,
        )
    if length < 1:
        raise ValueError("lag-kernel length is too short")

    if nfreq is None:
        nfreq = length + 1
    if nfreq < length + 1:
        warnings.warn(
            "nfreq is too short; increasing it to lag-kernel length + 1",
            RuntimeWarning,
            stacklevel=2,
        )
        nfreq = length + 1
    if nfreq % 2:
        nfreq += 1
        warnings.warn(
            f"forcing nfreq to the even value {nfreq}", RuntimeWarning, stacklevel=2
        )

    kernel = get_window(length, name, parameter, dft_window=in_frequency_domain)
    if len(specification) > 4 and str(specification[4]).lower() == "y":
        kernel = kernel / np.sum(kernel)
    return np.asarray(kernel), length, nfreq


def quadratic_tfd(
    signal: ArrayLike,
    doppler_kernel: Sequence[object],
    lag_kernel: Sequence[object],
    ntime: int | None = None,
    nfreq: int | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.generic], NDArray[np.generic]]:
    """Compute a quadratic TFD with a separable Doppler-lag kernel."""
    analytic, n2, signal_length, _ = get_analytic_signal(signal)
    lag, lag_length, nfreq = generate_lag_kernel(lag_kernel, signal_length, nfreq)
    doppler, doppler_length, ntime = generate_doppler_kernel(
        doppler_kernel, signal_length, ntime
    )

    half_freq = int(np.ceil(nfreq / 2))
    half_lag = int(np.ceil(lag_length / 2))
    half_doppler = int(np.ceil(doppler_length / 2))
    tfd = np.zeros((ntime, nfreq), dtype=np.float64)

    sample_indices = np.arange(signal_length)
    positive_doppler = np.arange(half_doppler)
    negative_doppler = np.arange(1, half_doppler)
    for lag_index in range(half_lag):
        forward = (sample_indices + lag_index) % n2
        backward = (sample_indices - lag_index) % n2
        lag_slice = lag[lag_index] * analytic[forward] * np.conj(analytic[backward])
        lag_spectrum = np.fft.fft(lag_slice)

        smoothed = np.zeros(ntime, dtype=np.complex128)
        smoothed[positive_doppler] = (
            lag_spectrum[positive_doppler] * doppler[positive_doppler]
        )
        if negative_doppler.size:
            smoothed[-negative_doppler] = (
                lag_spectrum[-negative_doppler]
                * doppler[doppler_length - negative_doppler]
            )
        smoothed = np.fft.ifft(smoothed)
        tfd[:, lag_index] = np.real(smoothed)
        tfd[:, lag_index + half_freq] = np.imag(smoothed)

    positive_lag = np.arange(half_lag)
    negative_lag = np.arange(1, half_lag)
    for time_index in range(0, ntime, 2):
        even_half = tfd[time_index, :half_freq] + 1j * tfd[time_index, half_freq:]
        odd_half = tfd[time_index + 1, :half_freq] + 1j * tfd[time_index + 1, half_freq:]
        even_slice = np.zeros(nfreq, dtype=np.complex128)
        odd_slice = np.zeros(nfreq, dtype=np.complex128)
        even_slice[positive_lag] = even_half[positive_lag]
        odd_slice[positive_lag] = odd_half[positive_lag]
        if negative_lag.size:
            even_slice[-negative_lag] = np.conj(even_half[negative_lag])
            odd_slice[-negative_lag] = np.conj(odd_half[negative_lag])
        transformed = np.fft.fft(even_slice + 1j * odd_slice)
        tfd[time_index, :] = np.real(transformed)
        tfd[time_index + 1, :] = np.imag(transformed)

    return tfd / nfreq, lag, doppler


def cross_tfd(
    x: ArrayLike,
    y: ArrayLike,
    doppler_kernel: Sequence[object],
    lag_kernel: Sequence[object],
    ntime: int | None = None,
    nfreq: int | None = None,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Compute the separable-kernel cross TFD."""
    x_values = np.asarray(x).reshape(-1)
    y_values = np.asarray(y).reshape(-1)
    if x_values.size != y_values.size:
        raise ValueError("x and y must have the same length")
    if ntime is None:
        ntime = x_values.size
    if nfreq is None:
        nfreq = x_values.size // 2

    x_analytic = get_analytic_signal(x_values)[0] if np.isrealobj(x_values) else x_values
    y_analytic = get_analytic_signal(y_values)[0] if np.isrealobj(y_values) else y_values
    x_analytic = np.asarray(x_analytic, dtype=np.complex128)
    y_analytic = np.asarray(y_analytic, dtype=np.complex128)
    n2 = x_analytic.size
    signal_length = n2 // 2

    lag, _, nfreq = generate_lag_kernel(lag_kernel, signal_length, nfreq)
    doppler, _, ntime = generate_doppler_kernel(doppler_kernel, signal_length, ntime)
    lag = pad_window(lag, nfreq)
    doppler = pad_window(doppler, ntime)

    rows = max(ntime, signal_length)
    time_lag = np.zeros((rows, nfreq), dtype=np.complex128)
    half_freq = int(np.ceil(nfreq / 2))
    lag_indices = np.arange(-(half_freq - 1), half_freq + 1)
    wrapped_lags = lag_indices % nfreq
    keep = lag[wrapped_lags] != 0
    wrapped_lags = wrapped_lags[keep]
    kept_lags = lag_indices[keep]

    for time_index in range(signal_length):
        forward = (time_index + kept_lags) % n2
        backward = (time_index - kept_lags) % n2
        time_lag[time_index, wrapped_lags] = (
            x_analytic[forward]
            * np.conj(y_analytic[backward])
            * lag[wrapped_lags]
        )

    doppler_lag = np.fft.fft(time_lag, axis=0) * doppler[:, None]
    result = np.fft.ifft(doppler_lag, axis=0)
    if ntime > signal_length:
        result = result[:signal_length, :]
    xtfd = np.fft.fft(result, axis=1) / nfreq
    return xtfd, result
