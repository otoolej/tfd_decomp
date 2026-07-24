"""Cross-TFD phase-reconstruction decomposition method."""

from __future__ import annotations

import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid

from ._tfd import cross_tfd, quadratic_tfd
from ._tracks import find_tracks
from ._utils import get_bandwidth, scale_tfd
from .params import DecompParams


def phases_from_cross_tfd(
    tracks: list[NDArray[np.int64]], cross_distribution: ArrayLike
) -> list[NDArray[np.float64]]:
    """Sample instantaneous phase from a cross TFD along each track."""
    xtfd = np.asarray(cross_distribution)
    result: list[NDArray[np.float64]] = []
    for track in tracks:
        phase_track = np.zeros((track.shape[0], 2), dtype=np.float64)
        phase_track[:, 0] = track[:, 0]
        valid = track[:, 0] < xtfd.shape[0]
        phase_track[valid, 1] = np.angle(xtfd[track[valid, 0], track[valid, 1]])
        result.append(phase_track)
    return result


def synthesize_sinusoidal_model(
    tracks: list[NDArray[np.int64]],
    phases: list[NDArray[np.float64]] | None,
    distribution: ArrayLike,
    sampling_frequency: float = 1.0,
    correct_amplitude_bandwidth: bool = False,
) -> list[NDArray[np.float64]]:
    """Synthesize sinusoidal components from IF tracks and a TFD."""
    tfd = np.asarray(distribution)
    magnitude = np.abs(tfd) if np.iscomplexobj(tfd) else np.asarray(tfd, dtype=np.float64)
    ntime, nfreq = magnitude.shape
    frequency_scale = sampling_frequency / (2 * nfreq)
    components: list[NDArray[np.float64]] = []

    for track_index, track in enumerate(tracks):
        length = track.shape[0]
        amplitude = magnitude[track[:, 0], track[:, 1]].astype(np.float64, copy=True)
        # Keep MATLAB's one-based frequency-bin convention in the synthesis formula.
        frequency = (track[:, 1] + 1) * frequency_scale
        phase = np.zeros(length, dtype=np.float64)
        if phases is not None:
            phase = np.asarray(phases[track_index][:, 1], dtype=np.float64)

        if correct_amplitude_bandwidth and length:
            bandwidths = np.zeros(length, dtype=np.float64)
            usable: list[int] = []
            for sample_index, (time_index, frequency_index) in enumerate(track):
                bandwidth = get_bandwidth(
                    magnitude[time_index, :],
                    int(frequency_index),
                    10 ** (-3 / 10),
                )
                if bandwidth is not None:
                    bandwidths[sample_index] = bandwidth
                    if bandwidth > 1:
                        usable.append(sample_index)
            max_amplitude = float(np.max(amplitude))
            maximum_index = int(np.argmax(amplitude))
            reference_bandwidth = bandwidths[maximum_index]
            if reference_bandwidth > 0 and usable:
                normalized = bandwidths / reference_bandwidth
                amplitude[usable] *= normalized[usable]
                amplitude[amplitude > max_amplitude] = max_amplitude

        amplitude[amplitude < 0] = 0
        integrated_phase = (
            2
            * np.pi
            * cumulative_trapezoid(frequency, dx=1.0, initial=0.0)
            / sampling_frequency
            + phase
        )
        component = np.zeros(ntime, dtype=np.float64)
        component[track[:, 0]] = np.sqrt(amplitude) * np.cos(integrated_phase)
        components.append(component)
    return components


def decompose_cross_tfd(
    signal: ArrayLike, params: DecompParams, n_components: int | None
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    list[NDArray[np.float64]],
    NDArray[np.float64],
    NDArray[np.complex128] | None,
]:
    """Run the complete xTFD method."""
    x = np.asarray(signal, dtype=np.float64).reshape(-1)
    signal_length = x.size
    raw_tfd, lag_window, _ = quadratic_tfd(
        x,
        params.doppler_kernel,
        params.lag_kernel,
        ntime=signal_length,
        nfreq=params.nfreq,
    )
    tfd = scale_tfd(raw_tfd, lag_window, signal_length)
    tracks = find_tracks(tfd, 1.0, signal_length, params)[0]
    if n_components is not None:
        tracks = tracks[:n_components]
    if not tracks:
        if params.warn_no_tracks:
            warnings.warn("no tracks found in the TFD", RuntimeWarning, stacklevel=2)
        empty = np.empty(0, dtype=np.float64)
        return empty, empty.copy(), [], tfd, None

    without_phase = synthesize_sinusoidal_model(
        tracks,
        None,
        tfd,
        sampling_frequency=1.0,
        correct_amplitude_bandwidth=False,
    )
    combined = np.sum(np.column_stack(without_phase), axis=1)

    phase_distribution: NDArray[np.complex128] | None = None
    if params.phase_correction:
        phase_distribution, _ = cross_tfd(
            x,
            combined,
            params.doppler_kernel,
            params.lag_kernel,
            ntime=signal_length,
            nfreq=params.nfreq,
        )
        phases = phases_from_cross_tfd(tracks, phase_distribution)
        components = synthesize_sinusoidal_model(
            tracks,
            phases,
            tfd,
            sampling_frequency=1.0,
            correct_amplitude_bandwidth=params.correct_amplitude_bw,
        )
    else:
        components = without_phase

    output = np.nansum(np.vstack(components), axis=0)
    residual = x - output
    return output, residual, components, tfd, phase_distribution
