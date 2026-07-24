import numpy as np

from tfd_decomp import DecompParams
from tfd_decomp._tfd import quadratic_tfd
from tfd_decomp._utils import (
    get_analytic_signal,
    get_window,
    pad_window,
    peak_pick,
    preprocess_signal,
    scale_tfd,
)


def test_analytic_signal_has_matlab_zero_padding():
    source = np.arange(8, dtype=float)
    analytic, doubled_length, signal_length, half = get_analytic_signal(source)

    assert analytic.shape == (16,)
    assert doubled_length == 16
    assert signal_length == 8
    assert half == 4
    np.testing.assert_allclose(analytic.real[:8], source, atol=1e-12)
    np.testing.assert_array_equal(analytic[8:], 0)


def test_shifted_window_and_odd_padding_layout():
    raw_hamming = np.array([0.08, 0.54, 1.0, 0.54, 0.08])
    shifted = get_window(5, "hamm")
    np.testing.assert_allclose(shifted, np.roll(raw_hamming, 3), atol=1e-15)

    padded = pad_window(np.array([1.0, 2.0, 3.0]), 7)
    np.testing.assert_array_equal(padded, [1.0, 2.0, 0, 0, 0, 0, 3.0])


def test_peak_picker_uses_strict_interior_maxima():
    peaks = peak_pick([3, 0, 2, 1, 1, 4, 0, 8])
    np.testing.assert_array_equal(np.flatnonzero(peaks), [2, 5])


def test_quadratic_tfd_dimensions_and_scaling_formula():
    source = np.cos(2 * np.pi * 0.1 * np.arange(32))
    params = DecompParams(source.size, "tvfilt")
    raw, lag, doppler = quadratic_tfd(
        source, params.doppler_kernel, params.lag_kernel, 32, 64
    )
    scaled = scale_tfd(raw, lag, source.size)

    assert raw.shape == (32, 64)
    assert lag.shape == (11,)
    assert doppler.shape == (23,)
    expected_factor = (64 / 32) * 32 / np.sum(lag)
    np.testing.assert_allclose(scaled, raw * expected_factor)


def test_positive_windowed_padding_adds_one_sixteenth_per_side():
    source = np.arange(32, dtype=float)
    padded = preprocess_signal(source, True)

    assert padded.shape == (36,)
    np.testing.assert_array_equal(padded[2:-2], source)
