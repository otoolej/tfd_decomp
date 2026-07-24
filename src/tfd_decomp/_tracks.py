"""Instantaneous-frequency track extraction."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._utils import peak_pick


IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class TrackParameters:
    """Small internal parameter view used by the McAulay-Quatieri linker."""

    delta_freq_samples: int
    min_if_length: int
    max_no_peaks: int


def extract_if_tracks(
    distribution: ArrayLike, parameters: TrackParameters
) -> tuple[NDArray[np.float64], list[IntArray]]:
    """Link time-slice peaks with the McAulay-Quatieri birth/death method.

    Returned tracks use zero-based ``[time, frequency_bin]`` coordinates.
    """
    values = np.asarray(distribution)
    if values.ndim != 2:
        raise ValueError("the time-frequency distribution must be two-dimensional")
    ntime, nfreq = values.shape
    states = np.zeros((ntime, nfreq), dtype=np.float64)
    track_ids = np.full((ntime, nfreq), -1, dtype=np.int64)
    individual: list[list[list[int]]] = []

    previous_peaks = np.empty(0, dtype=np.int64)
    for time_index in range(ntime):
        time_slice = values[time_index]
        magnitude = np.abs(time_slice) ** 2 if np.iscomplexobj(values) else time_slice
        peaks = np.flatnonzero(peak_pick(magnitude)).astype(np.int64)
        if peaks.size > parameters.max_no_peaks:
            strongest = np.argsort(-magnitude[peaks], kind="stable")[
                : parameters.max_no_peaks
            ]
            peaks = np.sort(peaks[strongest])

        if time_index == 0:
            for frequency in peaks:
                track_id = len(individual)
                individual.append([[0, int(frequency)]])
                track_ids[0, frequency] = track_id
            previous_peaks = peaks.copy()
            continue

        previous = previous_peaks.copy()
        current = peaks.copy()

        for previous_index in range(previous.size):
            if previous[previous_index] < 0:
                continue
            available = np.flatnonzero(current >= 0)
            if available.size == 0:
                _kill_track(
                    states,
                    track_ids,
                    individual,
                    time_index,
                    previous,
                    previous_index,
                )
                continue

            distances = np.abs(previous[previous_index] - current[available])
            relative = int(np.argmin(distances))
            current_index = int(available[relative])
            smallest_distance = int(distances[relative])

            if smallest_distance > parameters.delta_freq_samples:
                _kill_track(
                    states,
                    track_ids,
                    individual,
                    time_index,
                    previous,
                    previous_index,
                )
            else:
                candidate = int(current[current_index])
                future = np.flatnonzero(previous[previous_index + 1 :] >= 0)
                if future.size:
                    future_indices = future + previous_index + 1
                    future_distance = np.min(np.abs(previous[future_indices] - candidate))
                else:
                    future_distance = np.inf

                if future_distance >= smallest_distance:
                    _match_track(
                        states,
                        track_ids,
                        individual,
                        time_index,
                        previous,
                        current,
                        previous_index,
                        current_index,
                    )
                else:
                    # Preserve the original algorithm's immediate-next-track choice.
                    _match_track(
                        states,
                        track_ids,
                        individual,
                        time_index,
                        previous,
                        current,
                        previous_index + 1,
                        current_index,
                    )

            if previous[previous_index] >= 0 and current_index > 0:
                alternate_index = current_index - 1
                if current[alternate_index] >= 0:
                    alternate_distance = abs(
                        int(current[alternate_index]) - int(previous[previous_index])
                    )
                    if alternate_distance > parameters.delta_freq_samples:
                        _kill_track(
                            states,
                            track_ids,
                            individual,
                            time_index,
                            previous,
                            previous_index,
                        )
                    else:
                        _match_track(
                            states,
                            track_ids,
                            individual,
                            time_index,
                            previous,
                            current,
                            previous_index,
                            alternate_index,
                        )

        for previous_index in np.flatnonzero(previous >= 0):
            _kill_track(
                states,
                track_ids,
                individual,
                time_index,
                previous,
                int(previous_index),
            )

        for current_index in np.flatnonzero(current >= 0):
            _birth_track(
                states,
                track_ids,
                individual,
                time_index,
                current,
                int(current_index),
            )

        previous_peaks = peaks.copy()

    states[0, states[0] == 1] = 0.5
    tracks = [np.asarray(track, dtype=np.int64) for track in individual]
    return states, tracks


def _match_track(
    states: NDArray[np.float64],
    track_ids: IntArray,
    individual: list[list[list[int]]],
    time_index: int,
    previous: IntArray,
    current: IntArray,
    previous_index: int,
    current_index: int,
) -> None:
    if not (0 <= previous_index < previous.size and 0 <= current_index < current.size):
        return
    previous_frequency = int(previous[previous_index])
    current_frequency = int(current[current_index])
    if previous_frequency < 0 or current_frequency < 0:
        return

    states[time_index - 1, previous_frequency] = 1
    states[time_index, current_frequency] = 1
    track_id = int(track_ids[time_index - 1, previous_frequency])
    if track_id < 0:
        return
    track_ids[time_index, current_frequency] = track_id
    individual[track_id].append([time_index, current_frequency])
    previous[previous_index] = -1
    current[current_index] = -1


def _kill_track(
    states: NDArray[np.float64],
    track_ids: IntArray,
    individual: list[list[list[int]]],
    time_index: int,
    previous: IntArray,
    previous_index: int,
) -> None:
    if not 0 <= previous_index < previous.size:
        return
    frequency = int(previous[previous_index])
    if frequency < 0:
        return
    states[time_index - 1, frequency] = 1
    states[time_index, frequency] = -1
    track_id = int(track_ids[time_index - 1, frequency])
    if track_id >= 0:
        track_ids[time_index, frequency] = track_id
        individual[track_id].append([time_index, frequency])
    previous[previous_index] = -1


def _birth_track(
    states: NDArray[np.float64],
    track_ids: IntArray,
    individual: list[list[list[int]]],
    time_index: int,
    current: IntArray,
    current_index: int,
) -> None:
    frequency = int(current[current_index])
    if frequency < 0:
        return
    states[time_index - 1, frequency] = 0.5
    states[time_index, frequency] = 1
    track_id = len(individual)
    track_ids[time_index - 1, frequency] = track_id
    track_ids[time_index, frequency] = track_id
    individual.append(
        [[time_index - 1, frequency], [time_index, frequency]]
    )
    current[current_index] = -1


def find_tracks(
    distribution: ArrayLike,
    sampling_frequency: float,
    signal_length: int,
    parameters: object,
    frequency_limits: tuple[float, float] | None = None,
) -> tuple[list[IntArray], float, float, list[IntArray], NDArray[np.float64]]:
    """Extract, filter, and energy-rank instantaneous-frequency tracks."""
    tfd = np.asarray(distribution)
    if tfd.ndim != 2:
        raise ValueError("the time-frequency distribution must be two-dimensional")
    ntime, nfreq = tfd.shape
    duration = signal_length / sampling_frequency
    time_scale = duration / ntime
    frequency_scale = sampling_frequency / (2 * nfreq)
    delta = int(np.floor(parameters.delta_freq_samples * (nfreq / ntime)))

    if frequency_limits is not None:
        lower, upper = frequency_limits
        if lower < 0 or upper <= lower or upper > sampling_frequency / 2:
            raise ValueError("frequency limits must lie within [0, sampling_frequency / 2]")
        first = max(0, int(np.ceil(lower / frequency_scale)))
        last = min(nfreq - 1, int(np.floor(upper / frequency_scale)))
        narrowed = np.zeros_like(tfd)
        narrowed[:, first : last + 1] = tfd[:, first : last + 1]
        used = narrowed
    else:
        used = tfd

    link_parameters = TrackParameters(
        delta_freq_samples=delta,
        min_if_length=int(parameters.min_if_length),
        max_no_peaks=int(parameters.max_no_peaks),
    )
    _, raw_tracks = extract_if_tracks(used, link_parameters)
    if not raw_tracks:
        return [], frequency_scale, time_scale, [], np.empty(0, dtype=np.float64)

    retained: list[IntArray] = []
    energies: list[float] = []
    for track in raw_tracks:
        if track.shape[0] < parameters.min_if_length:
            continue
        energy = float(np.sum(tfd[track[:, 0], track[:, 1]]))
        retained.append(track.copy())
        energies.append(energy)

    if not energies:
        return [], frequency_scale, time_scale, raw_tracks, np.empty(0, dtype=np.float64)

    order = np.argsort(-np.asarray(energies), kind="stable")
    ordered_tracks: list[IntArray] = []
    ordered_energies: list[float] = []
    for index in order:
        energy = energies[int(index)]
        if energy > 0:
            ordered_tracks.append(retained[int(index)])
            ordered_energies.append(energy)
    return (
        ordered_tracks,
        frequency_scale,
        time_scale,
        raw_tracks,
        np.asarray(ordered_energies, dtype=np.float64),
    )


def tracks_for_tv_filtering(
    distribution: ArrayLike,
    signal_length: int,
    track_padding: int,
    sampling_frequency: float,
    bandwidth: int,
    min_if_length: int,
    max_peaks: int,
    max_threshold: float | None,
    frequency_limits: tuple[float, float] | None = None,
) -> tuple[list[IntArray], NDArray[np.float64]]:
    """Threshold, extract, extend, and mask tracks for TV filtering."""
    tfd = np.asarray(distribution, dtype=np.float64).copy()
    if max_threshold is not None and tfd.size:
        threshold = max_threshold * np.max(tfd)
        tfd[tfd < threshold] = 0

    parameters = TrackParameters(
        delta_freq_samples=bandwidth,
        min_if_length=min_if_length,
        max_no_peaks=max_peaks,
    )
    tracks = find_tracks(
        tfd, sampling_frequency, signal_length, parameters, frequency_limits
    )[0]
    full_length = signal_length + track_padding
    if not tracks:
        return [], np.zeros((tfd.shape[1], full_length), dtype=np.float64)

    tracks = [track.copy() for track in tracks]
    starts = np.asarray([track[0, 0] for track in tracks])
    earliest_index = int(np.argmin(starts))
    earliest = tracks[earliest_index]
    if earliest[0, 0] > 0:
        prefix_times = np.arange(earliest[0, 0], dtype=np.int64)
        prefix = np.column_stack(
            (prefix_times, np.full(prefix_times.size, earliest[0, 1], dtype=np.int64))
        )
        tracks[earliest_index] = np.vstack((prefix, earliest))

    for index, track in enumerate(tracks):
        end = int(track[-1, 0])
        if end + track_padding >= signal_length:
            extra_times = np.arange(end + 1, end + track_padding + 1, dtype=np.int64)
            extra = np.column_stack(
                (extra_times, np.full(extra_times.size, track[-1, 1], dtype=np.int64))
            )
            tracks[index] = np.vstack((track, extra))

    ends = np.asarray([track[-1, 0] for track in tracks])
    latest_index = int(np.argmax(ends))
    target_end = full_length - 1
    latest = tracks[latest_index]
    if latest[-1, 0] < target_end:
        extra_times = np.arange(latest[-1, 0] + 1, target_end + 1, dtype=np.int64)
        extra = np.column_stack(
            (extra_times, np.full(extra_times.size, latest[-1, 1], dtype=np.int64))
        )
        tracks[latest_index] = np.vstack((latest, extra))

    mask = np.zeros((tfd.shape[1], full_length), dtype=np.float64)
    for track in tracks:
        valid = track[:, 0] < full_length
        mask[track[valid, 1], track[valid, 0]] = 1
    return tracks, mask
