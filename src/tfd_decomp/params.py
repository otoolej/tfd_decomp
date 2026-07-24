"""Parameter model for the TFD decomposition methods."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, floor, sqrt
from typing import Any


KernelSpec = tuple[Any, ...]


def make_odd(value: float | int) -> int:
    """Match the MATLAB helper: floor a value and lower even values by one."""
    result = floor(value)
    if result % 2 == 0:
        result -= 1
    return result


@dataclass
class DecompParams:
    """Configuration for :func:`tfd_decomposition`.

    Defaults intentionally follow ``decomp_params.m``.  Parameters may be
    changed after construction; they are validated when decomposition starts.
    """

    signal_length: int
    method: str = "tvfilt"

    max_no_peaks: int = 8
    delta_freq_samples: int | None = None
    min_if_length: int | None = None

    doppler_kernel: KernelSpec | None = None
    lag_kernel: KernelSpec | None = None
    doppler_kernel_type: KernelSpec = ("hamm",)
    lag_kernel_type: KernelSpec = ("dolph", 100.0)

    pad_signal: bool = False
    low_pass_filter: bool = False
    warn_no_tracks: bool = True

    correct_amplitude_bw: bool = True
    phase_correction: bool = True
    nfreq: int = 256 * 32

    filter_length: int | None = None
    qtfd_max_threshold: float | None = None

    def __post_init__(self) -> None:
        self.method = str(self.method).lower()
        if not isinstance(self.signal_length, int) or self.signal_length < 2:
            raise ValueError("signal_length must be an integer of at least 2")
        if self.method not in {"tvfilt", "xtfd"}:
            raise ValueError("method must be either 'tvfilt' or 'xtfd'")

        root_n = sqrt(self.signal_length)
        if self.delta_freq_samples is None:
            self.delta_freq_samples = floor(root_n / 2)
        if self.min_if_length is None:
            self.min_if_length = floor(4 * root_n)
        if self.doppler_kernel is None:
            self.set_doppler_kernel(make_odd(ceil(4 * root_n)))
        else:
            self.doppler_kernel = tuple(self.doppler_kernel)
        if self.lag_kernel is None:
            self.set_lag_kernel(make_odd(ceil(2 * root_n)))
        else:
            self.lag_kernel = tuple(self.lag_kernel)

        if self.method == "tvfilt":
            if self.filter_length is None:
                self.filter_length = make_odd(ceil(2 * root_n))
            if self.qtfd_max_threshold is None:
                self.qtfd_max_threshold = 0.01

    def set_doppler_kernel(
        self, length: int, parameters: KernelSpec | None = None
    ) -> "DecompParams":
        """Set the Doppler kernel and return this object for convenient chaining."""
        values = self.doppler_kernel_type if parameters is None else tuple(parameters)
        self.doppler_kernel = (int(length), *values)
        return self

    def set_lag_kernel(
        self, length: int, parameters: KernelSpec | None = None
    ) -> "DecompParams":
        """Set the lag kernel and return this object for convenient chaining."""
        values = self.lag_kernel_type if parameters is None else tuple(parameters)
        self.lag_kernel = (int(length), *values)
        return self

    # MATLAB-compatible method spelling.
    set_dopp_kernel = set_doppler_kernel

    def validate_for(self, signal_length: int, method: str) -> None:
        """Validate mutable fields immediately before computation."""
        if self.signal_length != signal_length:
            raise ValueError(
                "params.signal_length must equal the unprocessed input signal length"
            )
        if self.method != method:
            raise ValueError("params.method must match the requested decomposition method")
        if self.low_pass_filter:
            raise NotImplementedError(
                "low_pass_filter has no filter definition in the MATLAB implementation"
            )
        if not isinstance(self.max_no_peaks, int) or self.max_no_peaks < 1:
            raise ValueError("max_no_peaks must be a positive integer")
        if not isinstance(self.delta_freq_samples, int) or self.delta_freq_samples < 0:
            raise ValueError("delta_freq_samples must be a non-negative integer")
        if not isinstance(self.min_if_length, int) or self.min_if_length < 1:
            raise ValueError("min_if_length must be a positive integer")

        self._validate_kernel("doppler_kernel", self.doppler_kernel)
        lag_length = self._validate_kernel("lag_kernel", self.lag_kernel)

        if method == "tvfilt":
            if not isinstance(self.filter_length, int) or self.filter_length < 1:
                raise ValueError("filter_length must be a positive integer")
            if self.filter_length % 2 == 0:
                raise ValueError("filter_length must be odd")
            if self.qtfd_max_threshold is None or not 0 <= self.qtfd_max_threshold <= 1:
                raise ValueError("qtfd_max_threshold must be between 0 and 1")
        else:
            if not isinstance(self.nfreq, int) or self.nfreq < lag_length + 1:
                raise ValueError("nfreq must be at least lag-kernel length + 1")
            if self.nfreq % 2:
                raise ValueError("nfreq must be even")

    @staticmethod
    def _validate_kernel(name: str, kernel: KernelSpec | None) -> int:
        if kernel is None or len(kernel) < 2:
            raise ValueError(f"{name} must contain at least a length and window name")
        length = kernel[0]
        if not isinstance(length, int) or length < 1:
            raise ValueError(f"{name} length must be a positive integer")
        if not isinstance(kernel[1], str):
            raise ValueError(f"{name} window name must be a string")
        return length

    # Familiar aliases for users translating existing MATLAB scripts.
    @property
    def N(self) -> int:
        return self.signal_length

    @property
    def Nfreq(self) -> int:
        return self.nfreq

    @Nfreq.setter
    def Nfreq(self, value: int) -> None:
        self.nfreq = value

    @property
    def L_filt(self) -> int | None:
        return self.filter_length

    @L_filt.setter
    def L_filt(self, value: int) -> None:
        self.filter_length = value

    @property
    def qtfd_max_thres(self) -> float | None:
        return self.qtfd_max_threshold

    @qtfd_max_thres.setter
    def qtfd_max_thres(self, value: float) -> None:
        self.qtfd_max_threshold = value

    @property
    def db_warn(self) -> bool:
        return self.warn_no_tracks

    @db_warn.setter
    def db_warn(self, value: bool) -> None:
        self.warn_no_tracks = value
