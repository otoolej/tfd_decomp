import pytest

from tfd_decomp import DecompParams


def test_tv_filter_defaults_match_matlab_formulas():
    params = DecompParams(400, "tvfilt")

    assert params.delta_freq_samples == 10
    assert params.min_if_length == 80
    assert params.doppler_kernel == (79, "hamm")
    assert params.lag_kernel == (39, "dolph", 100.0)
    assert params.filter_length == 39
    assert params.qtfd_max_threshold == 0.01


def test_xtfd_aliases_update_python_fields():
    params = DecompParams(256, "xtfd")
    params.Nfreq = 32768
    params.L_filt = 17
    params.qtfd_max_thres = 0.02
    params.db_warn = False

    assert params.nfreq == 32768
    assert params.filter_length == 17
    assert params.qtfd_max_threshold == 0.02
    assert not params.warn_no_tracks
    assert params.N == 256


def test_low_pass_option_fails_instead_of_silently_doing_nothing():
    params = DecompParams(64, "tvfilt")
    params.low_pass_filter = True

    with pytest.raises(NotImplementedError, match="no filter definition"):
        params.validate_for(64, "tvfilt")


@pytest.mark.parametrize("method", ["bad", "stft", ""])
def test_invalid_method_is_rejected(method):
    with pytest.raises(ValueError, match="tvfilt.*xtfd"):
        DecompParams(64, method)
