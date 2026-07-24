import numpy as np

from tfd_decomp._tracks import TrackParameters, extract_if_tracks, find_tracks


def test_constant_ridge_becomes_one_contiguous_track():
    distribution = np.zeros((12, 16))
    distribution[:, 5] = 2.0
    parameters = TrackParameters(
        delta_freq_samples=1, min_if_length=4, max_no_peaks=3
    )

    states, raw = extract_if_tracks(distribution, parameters)
    tracks, _, _, _, energies = find_tracks(distribution, 1.0, 12, parameters)

    assert states.shape == distribution.shape
    assert len(raw) == 1
    assert len(tracks) == 1
    np.testing.assert_array_equal(tracks[0][:, 0], np.arange(12))
    np.testing.assert_array_equal(tracks[0][:, 1], 5)
    np.testing.assert_allclose(energies, [24.0])


def test_tracks_are_ranked_by_accumulated_tfd_energy():
    distribution = np.zeros((10, 12))
    distribution[:, 3] = 1.0
    distribution[:, 8] = 3.0
    parameters = TrackParameters(
        delta_freq_samples=1, min_if_length=4, max_no_peaks=4
    )

    tracks, _, _, _, energies = find_tracks(distribution, 1.0, 10, parameters)

    assert len(tracks) == 2
    assert np.all(tracks[0][:, 1] == 8)
    assert np.all(tracks[1][:, 1] == 3)
    np.testing.assert_allclose(energies, [30.0, 10.0])


def test_flat_distribution_has_no_tracks():
    distribution = np.zeros((10, 12))
    parameters = TrackParameters(
        delta_freq_samples=1, min_if_length=4, max_no_peaks=4
    )

    tracks = find_tracks(distribution, 1.0, 10, parameters)[0]

    assert tracks == []
