"""Finite-field symplectic vectors, codes and coordinate operations (Secs. 1-2)."""

from galois import GF

from .symplectic import SymplecticVector, VElement, e, f, omega, zero
from .linear_code import LinearCode, linear_code
from .anticodes import anticode, anticodes
from .css import css_construct, reverse_css

__all__ = ["GF", "VElement", "SymplecticVector", "zero", "e", "f", "omega",
           "LinearCode", "linear_code", "anticode", "anticodes",
           "css_construct", "reverse_css"]

