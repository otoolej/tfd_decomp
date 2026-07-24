"""Matplotlib diagnostics for decomposition results."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def plot_tfd(distribution: ArrayLike, title: str) -> None:
    """Display a time-frequency distribution."""
    import matplotlib.pyplot as plt

    values = np.asarray(distribution)
    figure, axes = plt.subplots()
    axes.imshow(
        values,
        origin="lower",
        aspect="auto",
        extent=(0.0, 0.5, 0, values.shape[0] - 1),
    )
    axes.set_xlabel("normalized frequency")
    axes.set_ylabel("sample")
    axes.set_title(title)
    figure.tight_layout()


def plot_decomposition(
    source: ArrayLike,
    output: ArrayLike,
    components: list[np.ndarray],
    method: str,
) -> None:
    """Plot reconstructed signal, residual, and individual components."""
    import matplotlib.pyplot as plt

    x = np.asarray(source).reshape(-1)
    y = np.asarray(output).reshape(-1)
    figure, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(x, label="signal")
    axes[0].plot(y, label="reconstruction")
    axes[0].set_title(method)
    axes[0].legend()
    axes[1].plot(x - y)
    axes[1].set_title("error")
    axes[1].set_xlabel("sample")
    figure.tight_layout()

    if components:
        component_figure, component_axes = plt.subplots(
            len(components) + 1, 1, sharex=True, squeeze=False
        )
        component_axes[0, 0].plot(x, color="black")
        component_axes[0, 0].set_ylabel("signal")
        for index, component in enumerate(components, start=1):
            component_axes[index, 0].plot(component)
            component_axes[index, 0].set_ylabel(f"comp. {index}")
        component_axes[-1, 0].set_xlabel("sample")
        component_figure.tight_layout()
