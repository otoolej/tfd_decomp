"""Public decomposition API."""

from __future__ import annotations

import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._plotting import plot_decomposition, plot_tfd
from ._tvfilt import decompose_tv_filter
from ._utils import preprocess_signal
from ._xtfd import decompose_cross_tfd
from .params import DecompParams


def tfd_decomposition(
    x: ArrayLike,
    method: str = "tvfilt",
    n_components: int | None = 10,
    params: DecompParams | None = None,
    plot: bool = False,
) -> tuple[NDArray[np.float64], list[NDArray[np.float64]]]:
    """Decompose a real signal using TV filtering or cross-TFD phase synthesis.

    Parameters
    ----------
    x:
        One-dimensional real-valued input signal.
    method:
        ``"tvfilt"`` (default) or ``"xtfd"``; matching is case-insensitive.
    n_components:
        Maximum number of components, or ``None`` to retain all detected tracks.
    params:
        Optional :class:`DecompParams` built for this signal and method.
    plot:
        Display diagnostic Matplotlib figures when true.

    Returns
    -------
    y, components:
        The component sum and a list of component signals. All arrays have shape
        ``(N,)``. If no tracks are detected, ``y`` and the list are empty.
    """
    method_name = str(method).lower()
    if method_name not in {"tvfilt", "xtfd"}:
        raise ValueError("method must be either 'tvfilt' or 'xtfd'")
    if n_components is not None:
        if isinstance(n_components, bool) or not isinstance(n_components, int):
            raise TypeError("n_components must be a positive integer or None")
        if n_components < 1:
            raise ValueError("n_components must be a positive integer or None")

    values = np.asarray(x)
    if values.ndim != 1:
        raise ValueError("x must be a one-dimensional signal")
    if np.iscomplexobj(values):
        raise ValueError("the public API accepts real-valued signals only")
    try:
        values = values.astype(np.float64, copy=False)
    except (TypeError, ValueError) as error:
        raise TypeError("x must contain numeric values") from error
    if values.size < 2:
        raise ValueError("x must contain at least two samples")
    if not np.all(np.isfinite(values)):
        raise ValueError("x must contain only finite values")

    original_length = values.size
    if params is None:
        params = DecompParams(original_length, method_name)
    elif not isinstance(params, DecompParams):
        raise TypeError("params must be a DecompParams instance")
    params.validate_for(original_length, method_name)

    if values.size % 2:
        warnings.warn(
            "odd-length signal: dropping the final sample",
            RuntimeWarning,
            stacklevel=2,
        )
        values = values[:-1]
    source = values.copy()
    processed = preprocess_signal(source, params.pad_signal)
    trim = (processed.size - source.size) // 2

    phase_distribution = None
    if method_name == "tvfilt":
        component_matrix, distribution = decompose_tv_filter(
            processed, params, n_components
        )
        components = [row.copy() for row in component_matrix]
        if components:
            output = np.nansum(component_matrix, axis=0)
        else:
            output = np.empty(0, dtype=np.float64)
            if params.warn_no_tracks:
                warnings.warn("no tracks found in the TFD", RuntimeWarning, stacklevel=2)
    else:
        output, _, components, distribution, phase_distribution = decompose_cross_tfd(
            processed, params, n_components
        )

    if trim and output.size:
        output = output[trim:-trim]
        components = [component[trim:-trim] for component in components]

    output = np.asarray(output, dtype=np.float64).reshape(-1)
    components = [np.asarray(component, dtype=np.float64).reshape(-1) for component in components]

    if plot and output.size:
        plot_tfd(distribution, f"{method_name} time-frequency distribution")
        if phase_distribution is not None:
            plot_tfd(np.angle(phase_distribution), "cross-TFD phase")
        plot_decomposition(source, output, components, method_name)
    return output, components
