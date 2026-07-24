import numpy as np
import pytest

from tfd_decomp import DecompParams, tfd_decomposition
from tfd_decomp._xtfd import decompose_cross_tfd


def chirp(length=64):
    samples = np.arange(length)
    return np.cos(
        2 * np.pi * (0.06 * samples + ((0.24 - 0.06) / (2 * length)) * samples**2)
    )


def test_tv_filter_returns_1d_components_whose_sum_is_y():
    source = chirp()
    params = DecompParams(source.size, "tvfilt")
    params.min_if_length = 8

    output, components = tfd_decomposition(source, "tvfilt", 3, params)

    assert output.shape == source.shape
    assert 1 <= len(components) <= 3
    assert all(component.shape == source.shape for component in components)
    np.testing.assert_allclose(output, np.nansum(np.vstack(components), axis=0))
    assert np.all(np.isfinite(output))


def test_xtfd_returns_1d_components_whose_sum_is_y():
    source = chirp()
    params = DecompParams(source.size, "xtfd")
    params.nfreq = 128
    params.min_if_length = 8

    output, components = tfd_decomposition(source, "xtfd", 3, params)

    assert output.shape == source.shape
    assert 1 <= len(components) <= 3
    assert all(component.shape == source.shape for component in components)
    np.testing.assert_allclose(output, np.nansum(np.vstack(components), axis=0))
    assert np.all(np.isfinite(output))


def test_xtfd_residual_identity_and_determinism():
    source = chirp()
    params = DecompParams(source.size, "xtfd")
    params.nfreq = 128
    params.min_if_length = 8

    first = decompose_cross_tfd(source, params, 2)
    second = decompose_cross_tfd(source, params, 2)

    np.testing.assert_allclose(first[1], source - first[0])
    np.testing.assert_array_equal(first[0], second[0])


@pytest.mark.parametrize("method", ["tvfilt", "xtfd"])
def test_zero_signal_returns_empty_result(method):
    source = np.zeros(32)
    params = DecompParams(source.size, method)
    if method == "xtfd":
        params.nfreq = 64
    params.warn_no_tracks = False

    output, components = tfd_decomposition(source, method, params=params)

    assert output.shape == (0,)
    assert components == []


def test_padding_is_removed_once_from_result():
    source = chirp()
    params = DecompParams(source.size, "tvfilt")
    params.min_if_length = 8
    params.pad_signal = True

    output, components = tfd_decomposition(source, params=params)

    assert output.shape == source.shape
    assert all(component.shape == source.shape for component in components)


def test_xtfd_padding_and_optional_corrections_preserve_contract():
    source = chirp()
    params = DecompParams(source.size, "xtfd")
    params.nfreq = 128
    params.min_if_length = 8
    params.pad_signal = True
    params.phase_correction = False
    params.correct_amplitude_bw = False

    output, components = tfd_decomposition(source, "xtfd", 2, params)

    assert output.shape == source.shape
    assert 1 <= len(components) <= 2
    assert all(component.shape == source.shape for component in components)


def test_component_limit_is_applied_before_tv_filtering():
    length = 96
    samples = np.arange(length)
    source = np.cos(2 * np.pi * 0.08 * samples) + np.cos(2 * np.pi * 0.28 * samples)
    params = DecompParams(length, "tvfilt")
    params.min_if_length = 8

    one_output, one_component = tfd_decomposition(source, "tvfilt", 1, params)
    all_output, all_components = tfd_decomposition(source, "tvfilt", None, params)

    assert len(one_component) == 1
    assert len(all_components) >= len(one_component)
    np.testing.assert_allclose(one_output, one_component[0])
    np.testing.assert_allclose(one_component[0], all_components[0])
    np.testing.assert_allclose(all_output, np.nansum(np.vstack(all_components), axis=0))


def test_odd_input_warns_and_drops_last_sample():
    source = chirp(65)
    params = DecompParams(source.size, "tvfilt")
    params.min_if_length = 8

    with pytest.warns(RuntimeWarning, match="dropping the final sample"):
        output, components = tfd_decomposition(source, params=params)

    assert output.shape == (64,)
    assert all(component.shape == (64,) for component in components)


def test_plotting_uses_noninteractive_matplotlib_backend(monkeypatch):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    source = chirp(32)
    params = DecompParams(source.size, "tvfilt")
    params.min_if_length = 6
    monkeypatch.setattr(plt, "show", lambda: None)

    output, _ = tfd_decomposition(source, params=params, plot=True)

    assert output.shape == source.shape
    assert len(plt.get_fignums()) >= 2
    plt.close("all")


@pytest.mark.parametrize(
    "source, error, message",
    [
        (np.zeros((2, 4)), ValueError, "one-dimensional"),
        (np.array([0.0, np.nan]), ValueError, "finite"),
        (np.array([1 + 1j, 2 + 0j]), ValueError, "real-valued"),
    ],
)
def test_invalid_input_is_rejected(source, error, message):
    with pytest.raises(error, match=message):
        tfd_decomposition(source)


def test_parameter_method_and_signal_length_must_match_call():
    source = chirp()

    with pytest.raises(ValueError, match="method must match"):
        tfd_decomposition(source, "xtfd", params=DecompParams(source.size, "tvfilt"))
    with pytest.raises(ValueError, match="signal_length"):
        tfd_decomposition(source, params=DecompParams(source.size + 2, "tvfilt"))
