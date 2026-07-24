"""Quadratic time-frequency-distribution signal decomposition."""

from .api import tfd_decomposition
from .params import DecompParams

__all__ = ["DecompParams", "tfd_decomposition"]
__version__ = "0.1.0"
