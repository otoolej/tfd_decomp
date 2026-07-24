from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from tfd_decomp import DecompParams, tfd_decomposition
from tfd_decomp._tfd import cross_tfd, quadratic_tfd
from tfd_decomp._tracks import find_tracks
from tfd_decomp._utils import get_analytic_signal, scale_tfd
from tfd_decomp._xtfd import phases_from_cross_tfd, synthesize_sinusoidal_model


FIXTURE = Path(__file__).parent / "fixtures" / "matlab_reference.mat"
DATASET_FIXTURE = Path(__file__).parent / "fixtures" / "matlab_dataset_reference.mat"


@pytest.fixture(scope="module")
def reference():
    if not FIXTURE.exists():
        pytest.skip(
            "MATLAB golden fixture is absent; run tests/matlab/generate_reference_fixtures.m"
        )
    return loadmat(FIXTURE, squeeze_me=False)


@pytest.fixture(scope="module")
def dataset_reference():
    if not DATASET_FIXTURE.exists():
        pytest.skip(
            "MATLAB dataset fixture is absent; run the MATLAB fixture generator"
        )
    return loadmat(DATASET_FIXTURE, simplify_cells=True)["dataset_references"]


def _rows(cell_array):
    return [np.asarray(value, dtype=np.float64).reshape(-1) for value in cell_array.ravel()]


def test_shared_numerical_stages_match_matlab(reference):
    source = reference["x"].reshape(-1)
    params = DecompParams(source.size, "xtfd")
    params.nfreq = 128
    params.min_if_length = 8
    raw, lag, doppler = quadratic_tfd(
        source, params.doppler_kernel, params.lag_kernel, source.size, params.nfreq
    )
    scaled = scale_tfd(raw, lag, source.size)
    analytic = get_analytic_signal(source)[0]

    np.testing.assert_allclose(
        analytic, reference["analytic_z"].reshape(-1), rtol=1e-10, atol=1e-12
    )
    np.testing.assert_allclose(lag, reference["g2"].reshape(-1), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(doppler, reference["G1"].reshape(-1), rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(raw, reference["qtfd_raw"], rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(scaled, reference["qtfd_scaled"], rtol=1e-10, atol=1e-12)


def test_tracks_match_matlab_after_index_conversion(reference):
    source = reference["x"].reshape(-1)
    params = DecompParams(source.size, "xtfd")
    params.nfreq = 128
    params.min_if_length = 8
    raw, lag, _ = quadratic_tfd(
        source, params.doppler_kernel, params.lag_kernel, source.size, params.nfreq
    )
    tracks = find_tracks(scale_tfd(raw, lag, source.size), 1.0, source.size, params)[0]
    matlab_tracks = [np.asarray(item, dtype=np.int64) - 1 for item in reference["tracks"].ravel()]

    assert len(tracks) == len(matlab_tracks)
    for actual, expected in zip(tracks, matlab_tracks, strict=True):
        np.testing.assert_array_equal(actual, expected)


def test_cross_tfd_phase_and_synthesis_match_matlab(reference):
    source = reference["x"].reshape(-1)
    params = DecompParams(source.size, "xtfd")
    params.nfreq = 128
    params.min_if_length = 8
    raw, lag, _ = quadratic_tfd(
        source, params.doppler_kernel, params.lag_kernel, source.size, params.nfreq
    )
    distribution = scale_tfd(raw, lag, source.size)
    tracks = find_tracks(distribution, 1.0, source.size, params)[0]
    no_phase = synthesize_sinusoidal_model(tracks, None, distribution)
    signal_no_phase = np.sum(np.column_stack(no_phase), axis=1)
    cross_distribution = cross_tfd(
        source,
        signal_no_phase,
        params.doppler_kernel,
        params.lag_kernel,
        source.size,
        params.nfreq,
    )[0]
    phases = phases_from_cross_tfd(tracks, cross_distribution)
    corrected = synthesize_sinusoidal_model(
        tracks, phases, distribution, correct_amplitude_bandwidth=True
    )

    np.testing.assert_allclose(
        signal_no_phase, reference["signal_no_phase"].reshape(-1), rtol=1e-8, atol=1e-10
    )
    np.testing.assert_allclose(
        cross_distribution, reference["cross_distribution"], rtol=1e-10, atol=1e-12
    )
    expected_phases = reference["phase_tracks"].ravel()
    expected_components = reference["components_with_phase"].ravel()
    for actual, expected in zip(phases, expected_phases, strict=True):
        np.testing.assert_allclose(actual[:, 1], expected[:, 1], rtol=1e-8, atol=1e-10)
    for actual, expected in zip(corrected, expected_components, strict=True):
        np.testing.assert_allclose(actual, expected.reshape(-1), rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("method", ["tvfilt", "xtfd"])
def test_final_decompositions_match_matlab(reference, method):
    source = reference["x"].reshape(-1)
    params = DecompParams(source.size, method)
    params.min_if_length = 8
    if method == "xtfd":
        params.nfreq = 128
    output, components = tfd_decomposition(source, method, 3, params)
    expected_output = reference[f"y_{method}"].reshape(-1)
    expected_components = _rows(reference[f"components_{method}"])

    np.testing.assert_allclose(output, expected_output, rtol=1e-8, atol=1e-10)
    assert len(components) == len(expected_components)
    for actual, expected in zip(components, expected_components, strict=True):
        np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize(
    "dataset_name", ["white_noise", "ffgn", "eeg", "eeg1", "bat", "whale"]
)
@pytest.mark.parametrize("method", ["tvfilt", "xtfd"])
def test_bundled_signal_outputs_match_matlab(dataset_reference, dataset_name, method):
    expected = dataset_reference[dataset_name]
    source = np.asarray(expected["x"]).reshape(-1)
    params = DecompParams(source.size, method)
    if method == "xtfd":
        params.nfreq = 1024
    output, components = tfd_decomposition(source, method, 3, params)
    expected_output = np.asarray(expected[f"y_{method}"]).reshape(-1)
    expected_components = expected[f"components_{method}"]
    if isinstance(expected_components, np.ndarray) and expected_components.dtype != object:
        expected_components = [expected_components]
    expected_components = [np.asarray(item).reshape(-1) for item in expected_components]

    np.testing.assert_allclose(output, expected_output, rtol=1e-8, atol=1e-10)
    assert len(components) == len(expected_components)
    for actual, wanted in zip(components, expected_components, strict=True):
        np.testing.assert_allclose(actual, wanted, rtol=1e-8, atol=1e-10)
