"""Time-varying-filter decomposition method."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import windows

from ._tfd import quadratic_tfd
from ._tracks import tracks_for_tv_filtering
from ._utils import scale_tfd
from .params import DecompParams


def time_varying_filter(
    tracks: list[NDArray[np.int64]],
    track_mask: ArrayLike,
    signal: ArrayLike,
    filter_length: int,
    n_components: int | None,
) -> NDArray[np.float64]:
    """Extract ridge-defined components using lossless time-varying FIR bands."""
    selected = tracks if n_components is None else tracks[:n_components]
    x = np.asarray(signal, dtype=np.float64).reshape(-1)
    x = np.concatenate((x, np.zeros(filter_length + 1, dtype=np.float64)))
    total_length = x.size
    mask = np.asarray(track_mask)
    nfreq = mask.shape[0]
    components = np.zeros((len(selected), total_length), dtype=np.float64)
    taper = windows.hamming(filter_length, sym=True)
    half = filter_length // 2
    tap_numbers = np.arange(filter_length)
    alpha = (filter_length - 1) / 2
    offsets = np.arange(-half, half + 1)

    for component_index, track in enumerate(selected):
        start = int(np.min(track[:, 0]))
        end = min(int(np.max(track[:, 0])), total_length - 1)
        center_frequencies = np.zeros(total_length, dtype=np.float64)
        lower = np.zeros(total_length, dtype=np.float64)
        upper = np.zeros(total_length, dtype=np.float64)

        valid_track = track[track[:, 0] < total_length]
        center_frequencies[valid_track[:, 0]] = (valid_track[:, 1] + 1) / (2 * nfreq)
        for time_index, frequency_index in valid_track:
            ridge_frequencies = np.flatnonzero(mask[:, time_index] == 1)
            offsets_from_center = ridge_frequencies - frequency_index
            below = offsets_from_center[offsets_from_center < 0]
            above = offsets_from_center[offsets_from_center > 0]
            if below.size:
                lower[time_index] = center_frequencies[time_index] + (
                    np.floor(below[-1] / 2) / (2 * nfreq)
                )
            else:
                lower[time_index] = 0.0
            if above.size:
                upper[time_index] = center_frequencies[time_index] + (
                    np.floor(above[0] / 2) / (2 * nfreq)
                )
            else:
                upper[time_index] = 0.5

        for time_index in range(start, end + 1):
            relative = tap_numbers - alpha
            coefficients = np.empty(filter_length, dtype=np.float64)
            center = relative == 0
            coefficients[center] = 2 * (upper[time_index] - lower[time_index])
            noncenter = ~center
            coefficients[noncenter] = (
                np.sin(2 * np.pi * upper[time_index] * relative[noncenter])
                - np.sin(2 * np.pi * lower[time_index] * relative[noncenter])
            ) / (np.pi * relative[noncenter])
            coefficients *= taper

            source_indices = time_index - offsets
            valid = (source_indices >= 0) & (source_indices < total_length)
            components[component_index, time_index] = np.sum(
                coefficients[valid] * x[source_indices[valid]]
            )
    return components


def decompose_tv_filter(
    signal: ArrayLike, params: DecompParams, n_components: int | None
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Run the complete TV-filter method and return components plus its TFD."""
    x = np.asarray(signal, dtype=np.float64).reshape(-1)
    signal_length = x.size
    raw_tfd, lag_window, _ = quadratic_tfd(
        x,
        params.doppler_kernel,
        params.lag_kernel,
        ntime=signal_length,
        nfreq=signal_length,
    )
    tfd = scale_tfd(raw_tfd, lag_window, signal_length)
    tracks, mask = tracks_for_tv_filtering(
        tfd,
        signal_length,
        params.filter_length + 1,
        1.0,
        params.delta_freq_samples,
        params.min_if_length,
        params.max_no_peaks,
        params.qtfd_max_threshold,
    )
    if not tracks:
        return np.empty((0, signal_length), dtype=np.float64), tfd
    components = time_varying_filter(
        tracks, mask, x, params.filter_length, n_components
    )
    return components[:, :signal_length], tfd
