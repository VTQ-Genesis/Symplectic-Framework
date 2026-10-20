from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from itertools import combinations

from galois import FieldArray

from ._coordinates import _coordinate_indices, _require_length
from .symplectic import (
    SymplecticVector,
    VElement,
    _require_field,
    _require_same_field,
    omega,
)


def _require_code_vector(vector: SymplecticVector, field: type[FieldArray],
                         n: int) -> None:
    """Check that a vector belongs to the selected ambient space V^n."""
    if not isinstance(vector, SymplecticVector):
        raise TypeError("code vectors must be SymplecticVector objects")
    _require_same_field(field, vector.field)
    if vector.n != n:
        raise ValueError(f"code vectors must have length {n}, got {vector.n}")


def _coefficient_matrix(vectors: tuple[SymplecticVector, ...],
                        field: type[FieldArray], n: int) -> FieldArray:
    """Copy framework vectors into a field matrix with 2*n columns."""
    matrix = field.Zeros((len(vectors), 2 * n))
    for i, vector in enumerate(vectors):
        matrix[i, :] = [coefficient for coordinate in vector
                        for coefficient in (coordinate.a, coordinate.b)]
    return matrix


def _vectors_from_rows(rows: FieldArray, field: type[FieldArray],
                       n: int) -> tuple[SymplecticVector, ...]:
    """Rebuild vectors from interleaved (a_1, b_1, ..., a_n, b_n) rows."""
    return tuple(
        SymplecticVector(
            (VElement(row[2 * j], row[2 * j + 1], field=field)
             for j in range(n)),
            field=field,
        )
        for row in rows
    )


@dataclass(frozen=True, slots=True, init=False)
class LinearCode:
    """The F-linear span of generators in the existing symplectic V^n.
    """

    basis: tuple[SymplecticVector, ...]
    field: type[FieldArray]
    n: int

    def __init__(self, generators: Iterable[SymplecticVector],
                 field: type[FieldArray], n: int) -> None:
        _require_field(field)
        _require_length(n)
        vectors = tuple(generators)
        for vector in vectors:
            _require_code_vector(vector, field, n)

        reduced = _coefficient_matrix(vectors, field, n).row_space()
        basis = _vectors_from_rows(reduced, field, n)
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "field", field)
        object.__setattr__(self, "n", n)

    @property
    def dimension_over_field(self) -> int:
        """Ordinary vector-space dimension dim_F(C), the number of basis rows."""
        return len(self.basis)

    def _radical_coefficients(self) -> FieldArray:
        """RREF basis of radical coefficients relative to self.basis.
        """
        size = self.dimension_over_field
        gram = self.field.Zeros((size, size))
        for i, left in enumerate(self.basis):
            for j, right in enumerate(self.basis):
                gram[i, j] = omega(left, right)
        return gram.left_null_space().row_reduce()

    def symplectic_decomposition(self) -> tuple[LinearCode, LinearCode]:
        """Return (rad(C), K) with C = rad(C) direct-sum K (symplectic)
        """
        coefficients = self._radical_coefficients()
        code_rows = _coefficient_matrix(self.basis, self.field, self.n)
        radical_rows = coefficients @ code_rows
        radical = LinearCode(
            _vectors_from_rows(radical_rows, self.field, self.n),
            self.field, self.n,
        )
        pivot_columns = {int(row.nonzero()[0][0]) for row in coefficients}
        complement = LinearCode(
            (vector for j, vector in enumerate(self.basis)
             if j not in pivot_columns),
            self.field, self.n,
        )
        return radical, complement

    def radical(self) -> LinearCode:
        """Return rad(C) = C intersect C-perp as a subspace
        """
        return self.symplectic_decomposition()[0]

    def symplectic_part(self) -> LinearCode:
        """Return the chosen nondegenerate complement K to rad(C).
        """
        return self.symplectic_decomposition()[1]

    @property
    def dimension(self) -> int:
        """Paper dimension dim(C): number of symplectic pairs
        """
        radical_dimension = self._radical_coefficients().shape[0]
        return (self.dimension_over_field - radical_dimension) // 2

    @property
    def isorank(self) -> int:
        """Maximum ordinary F-dimension of an isotropic subspace of C.
        """
        return self.dimension_over_field - self.dimension

    @property
    def cardinality(self) -> int:
        """Number of codewords, |F|**dim_F(C), without enumerating them."""
        return self.field.order ** self.dimension_over_field

    def puncture(self, indices: Iterable[int]) -> LinearCode:
        """Return pi_J(C), keeping the coordinates in J, then getting rid of those not.
        """
        # ensure the indices are in the appropriate python format and are valid for the length
        positions = _coordinate_indices(indices, self.n)
        # take in the generator and puncutre each basis vector and store the basis
        projected = (
            SymplecticVector((vector[i] for i in positions), field=self.field)
            for vector in self.basis
        )
        return LinearCode(projected, self.field, len(positions))

    
    def restrict(self, indices: Iterable[int]) -> LinearCode:
        """Return C intersect A_J at the original length n. Thus its essentially
        """
        # keep the code words supported on the J, then return that as new code
        # first ensure that the indices are valid 
        positions = _coordinate_indices(indices, self.n)
        # convert the positions to a set for efficiency of checking inclusion
        selected = set(positions)
        # we store vectors according to their coefficients ((a1,b1),(a2,b2)...)
        # select the coefficient columns that are not in the J subset, ie those that will be omitted
        outside_columns = [2 * i + offset for i in range(self.n)
                           if i not in selected for offset in (0, 1)]
        # constructs the matrix where each row is the coefficient representation for each basis vec,
        # represents the generaor matrix for coefficients in F xG
        matrix = _coefficient_matrix(self.basis, self.field, self.n)

        # gives the coefficients that give the combinations of the outside columns to get 0.
        coefficients = matrix[:, outside_columns].left_null_space()

        # construct the new code that will have 0s not in J, as a LC of the matrix,
        # will be encoded with the 0's at the back as we have the coefficients as 
        # the null of the columns outside
        restricted_rows = coefficients @ matrix
        return LinearCode(
            _vectors_from_rows(restricted_rows, self.field, self.n),
            self.field, self.n,
        )

    def restrictions(self) -> Iterator[tuple[tuple[int, ...], LinearCode]]:
        """Yield (J, C intersect A_J) for every coordinate subset J.
        """
        # return the iterator of ((subset), words supported)
        return (
            (indices, self.restrict(indices))
            for size in range(self.n + 1)
            for indices in combinations(range(1, self.n + 1), size)
        )

    def shorten(self, indices: Iterable[int]) -> LinearCode:
        """Return pi_J(C intersect A_J), keeping J 
        """
        # Preserve a one-use iterable for both operations.
        indices = tuple(indices)
        # we select the codewords that are supported on J, then keep only the elements on those indices
        return self.restrict(indices).puncture(indices)

    def __contains__(self, vector: SymplecticVector) -> bool:
        """Test membership; reject wrong types, fields, or ambient lengths."""
        _require_code_vector(vector, self.field, self.n)
        augmented = _coefficient_matrix(self.basis + (vector,), self.field, self.n)
        return augmented.row_space().shape[0] == self.dimension_over_field

    def codewords(self) -> Iterator[SymplecticVector]:
        """Yield each codeword once, starting with zero.
        """
        # Base-q digits enumerate coefficient tuples without storing all q
        # scalars. The last basis coefficient changes fastest.
        for index in range(self.cardinality):
            word = SymplecticVector.zero(self.n, field=self.field)
            for vector in reversed(self.basis):
                index, scalar = divmod(index, self.field.order)
                word = word + vector * scalar
            yield word


def linear_code(generators: Iterable[SymplecticVector], field: type[FieldArray],
                n: int) -> LinearCode:
    """Construct C = span_F(generators) <= V^n with an independent basis.
    """
    return LinearCode(generators, field, n)
