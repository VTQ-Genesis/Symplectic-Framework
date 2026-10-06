"""Finite-field coordinate realization of Section 1, Definitions 1.1 and 1.2.

V = F e + F f, with omega(e, f) = 1. V^n consists of n mutually
orthogonal copies of V. The default field is F_2. This coordinate model
implements vectors and their form, not the paper's subspace theory.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field as dataclass_field
from typing import overload

from galois import GF, FieldArray


_BINARY_FIELD = GF(2)


def _require_field(field: type[FieldArray]) -> None:
    """Require a concrete finite-field class constructed by galois.GF."""
    if (not isinstance(field, type) or not issubclass(field, FieldArray)
            or field is FieldArray):
        raise TypeError("field must be a finite-field class returned by GF(q)")


def _coefficient(value: int | FieldArray, field: type[FieldArray], name: str) -> int:
    """Copy a field scalar into its immutable integer encoding."""
    if isinstance(value, FieldArray):
        if type(value) is not field:
            raise ValueError(f"{name} must belong to the selected field")
        if value.ndim != 0:
            raise TypeError(f"{name} must be a scalar, not a field array")
        return int(value)
    if type(value) is not int:
        raise TypeError(f"{name} must be a Python integer encoding or a field scalar")
    if not 0 <= value < field.order:
        raise ValueError(f"{name} must be in range 0 <= value < {field.order}")
    return value


def _require_same_field(left: type[FieldArray], right: type[FieldArray]) -> None:
    if left is not right:
        raise ValueError("operands must use the same finite-field class")


def _field_repr(field: type[FieldArray]) -> str:
    """Retain the exact field presentation when representing nonbinary values."""
    arguments = str(field.order)
    if field.degree > 1:
        arguments += f", irreducible_poly={int(field.irreducible_poly)}"
    arguments += f", primitive_element={int(field.primitive_element)}"
    return f"GF({arguments})"


@dataclass(frozen=True, slots=True)
class VElement:
    """The element a*e + b*f of V over the selected finite field.

    Use VElement(a, b, field=GF(q)); omitting field retains F_2. Coefficients
    are stored as integer encodings in [0, q), as used by galois, not as
    residues modulo q in extension fields. Matching galois scalars are also
    accepted and copied. Booleans, floats and implicit reduction are rejected.
    V is a two-dimensional vector space over F, not the field of order q^2.
    """

    a: int
    b: int
    field: type[FieldArray] = dataclass_field(default=_BINARY_FIELD, kw_only=True)

    # Let a NumPy/galois scalar on the left defer to our __rmul__ method.
    __array_ufunc__ = None

    def __post_init__(self) -> None:
        _require_field(self.field)
        object.__setattr__(self, "a", _coefficient(self.a, self.field, "a"))
        object.__setattr__(self, "b", _coefficient(self.b, self.field, "b"))

    def __repr__(self) -> str:
        suffix = "" if self.field is _BINARY_FIELD else f", field={_field_repr(self.field)}"
        return f"VElement(a={self.a}, b={self.b}{suffix})"

    def __str__(self) -> str:
        terms = []
        for coefficient, symbol in ((self.a, "e"), (self.b, "f")):
            if coefficient:
                terms.append(symbol if coefficient == 1 else f"{coefficient}*{symbol}")
        return "+".join(terms) or "0"

    def __add__(self, other: VElement) -> VElement:
        if not isinstance(other, VElement):
            return NotImplemented
        _require_same_field(self.field, other.field)
        if self.field is _BINARY_FIELD:
            return VElement(self.a ^ other.a, self.b ^ other.b)
        field = self.field
        return VElement(field(self.a) + field(other.a),
                        field(self.b) + field(other.b), field=field)

    def __neg__(self) -> VElement:
        if self.field.characteristic == 2:
            return self
        return VElement(-self.field(self.a), -self.field(self.b), field=self.field)

    def __sub__(self, other: VElement) -> VElement:
        if not isinstance(other, VElement):
            return NotImplemented
        return self + (-other)

    def __mul__(self, scalar: int | FieldArray) -> VElement:
        scalar = _coefficient(scalar, self.field, "scalar")
        if self.field is _BINARY_FIELD:
            return VElement(scalar * self.a, scalar * self.b)
        field = self.field
        return VElement(field(scalar) * field(self.a),
                        field(scalar) * field(self.b), field=field)

    def __rmul__(self, scalar: int | FieldArray) -> VElement:
        return self * scalar


# These constants deliberately keep their original F_2 meaning.
zero = VElement(0, 0)
e = VElement(1, 0)
f = VElement(0, 1)


@dataclass(frozen=True, slots=True, init=False)
class SymplecticVector:
    """An immutable tuple (v_1, ..., v_n) in the orthogonal direct sum V^n.

    Infer the field from VElement coordinates, or supply field=GF(q).
    Literal integer 0 denotes the zero coordinate in that field. Empty and
    all-literal-zero inputs default to F_2 unless a field is supplied.
    n is the coordinate count; the ordinary F dimension of V^n is 2*n.
    """

    coordinates: tuple[VElement, ...]
    field: type[FieldArray]

    
    __array_ufunc__ = None

    def __init__(self, coordinates: Iterable[VElement | int], *,
                 field: type[FieldArray] | None = None) -> None:
        entries = tuple(coordinates)
        if field is None:
            field = next((entry.field for entry in entries if isinstance(entry, VElement)),
                         _BINARY_FIELD)
        _require_field(field)
        normalized = []
        for coordinate in entries:
            if isinstance(coordinate, VElement):
                _require_same_field(field, coordinate.field)
                normalized.append(coordinate)
            elif type(coordinate) is int:
                if coordinate != 0:
                    raise ValueError("only the integer 0 denotes a V coordinate")
                normalized.append(VElement(0, 0, field=field))
            else:
                raise TypeError("coordinates must be VElement objects or integer 0")
        object.__setattr__(self, "coordinates", tuple(normalized))
        object.__setattr__(self, "field", field)

    @classmethod
    def zero(cls, n: int, *, field: type[FieldArray] = _BINARY_FIELD) -> SymplecticVector:
        """Return the zero vector of V^n over field, including () when n = 0."""
        if type(n) is not int:
            raise TypeError("n must be a Python integer")
        if n < 0:
            raise ValueError("n must be non-negative")
        return cls((0,) * n, field=field)

    @property
    def n(self) -> int:
        """Number of V coordinates (not the ordinary vector-space dimension)."""
        return len(self.coordinates)

    def __len__(self) -> int:
        return self.n

    def __iter__(self) -> Iterator[VElement]:
        return iter(self.coordinates)

    @overload
    def __getitem__(self, index: int) -> VElement: ...

    @overload
    def __getitem__(self, index: slice) -> tuple[VElement, ...]: ...

    def __getitem__(self, index: int | slice) -> VElement | tuple[VElement, ...]:
        return self.coordinates[index]

    def __repr__(self) -> str:
        suffix = "" if self.field is _BINARY_FIELD else f", field={_field_repr(self.field)}"
        return f"SymplecticVector(coordinates={self.coordinates!r}{suffix})"

    def __str__(self) -> str:
        entries = ", ".join(str(coordinate) for coordinate in self)
        if self.n == 1:
            entries += ","
        return f"({entries})"

    def _require_same_length(self, other: SymplecticVector) -> None:
        if self.n != other.n:
            raise ValueError(f"vectors must have the same length: {self.n} != {other.n}")

    def _require_compatible(self, other: SymplecticVector) -> None:
        self._require_same_length(other)
        _require_same_field(self.field, other.field)

    def __add__(self, other: SymplecticVector) -> SymplecticVector:
        if not isinstance(other, SymplecticVector):
            return NotImplemented
        self._require_compatible(other)
        return SymplecticVector((u + v for u, v in zip(self, other)), field=self.field)

    def __neg__(self) -> SymplecticVector:
        if self.field.characteristic == 2:
            return self
        return SymplecticVector((-coordinate for coordinate in self), field=self.field)

    def __sub__(self, other: SymplecticVector) -> SymplecticVector:
        if not isinstance(other, SymplecticVector):
            return NotImplemented
        return self + (-other)

    def __mul__(self, scalar: int | FieldArray) -> SymplecticVector:
        scalar = _coefficient(scalar, self.field, "scalar")
        return SymplecticVector((coordinate * scalar for coordinate in self), field=self.field)

    def __rmul__(self, scalar: int | FieldArray) -> SymplecticVector:
        return self * scalar


@overload
def omega(u: VElement, v: VElement) -> int | FieldArray: ...


@overload
def omega(u: SymplecticVector, v: SymplecticVector) -> int | FieldArray: ...


def omega(u: VElement | SymplecticVector,
          v: VElement | SymplecticVector) -> int | FieldArray:
    """Evaluate the field-valued form of Definitions 1.1-1.2.

    omega(ae+bf, ce+df) = ad-bc. For tuples, sum the local pairings in F,
    making distinct coordinate planes orthogonal. Operands must have matching
    types, fields and lengths. The empty sum on V^0 is the field's zero.
    Return a galois field scalar, except over F_2 where the original Python
    integer 0/1 return type is preserved. Use u.field(omega(u, v)) when a
    field scalar is needed uniformly across all fields.
    """
    if isinstance(u, VElement) and isinstance(v, VElement):
        _require_same_field(u.field, v.field)
        if u.field is _BINARY_FIELD:
            return (u.a * v.b + u.b * v.a) % 2
        field = u.field
        return field(u.a) * field(v.b) - field(u.b) * field(v.a)
    if isinstance(u, SymplecticVector) and isinstance(v, SymplecticVector):
        u._require_compatible(v)
        if u.field is _BINARY_FIELD:
            return sum(omega(left, right) for left, right in zip(u, v)) % 2
        return sum((omega(left, right) for left, right in zip(u, v)), u.field(0))
    raise TypeError("omega requires two VElement objects or two SymplecticVector objects")