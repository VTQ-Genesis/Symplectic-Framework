"""F-linear codes in V^n, following Section 1 and Definition 2.2.

V = F e + F f is the existing coordinate plane. Row reduction uses the
ordinary F coordinates (a_1, b_1, ..., a_n, b_n), not scalars from V.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from galois import FieldArray

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

    Supply a finite iterable of SymplecticVector objects, a field class from
    GF(q), and the nonnegative code length n. All vectors must already use
    that exact field and length. Empty input defines the zero code.

    basis stores the nonzero RREF rows as immutable framework vectors, so
    reordering or adding dependent generators gives the same code. Row
    reduction may replace the supplied vectors by linear combinations.
    No isotropy or nondegeneracy condition is imposed on the code.

    dimension_over_field means dim_F(C); dimension and isorank use
    Definition 1.7. radical() and symplectic_part() return the subspaces
    in an orthogonal splitting from Theorem 1.6. Distance is not computed.
    """

    basis: tuple[SymplecticVector, ...]
    field: type[FieldArray]
    n: int

    def __init__(self, generators: Iterable[SymplecticVector],
                 field: type[FieldArray], n: int) -> None:
        _require_field(field)
        if type(n) is not int:
            raise TypeError("n must be a Python integer")
        if n < 0:
            raise ValueError("n must be non-negative")
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

        For G[i,j] = omega(b_i,b_j), xG = 0 is equivalent to
        sum_i x_i*b_i being orthogonal to every vector of C (Def. 1.5).
        galois solves this homogeneous system over the code's field.
        Explicit row reduction makes the pivot convention independent
        of the basis ordering returned by the null-space routine.
        """
        size = self.dimension_over_field
        gram = self.field.Zeros((size, size))
        for i, left in enumerate(self.basis):
            for j, right in enumerate(self.basis):
                gram[i, j] = omega(left, right)
        return gram.left_null_space().row_reduce()

    def symplectic_decomposition(self) -> tuple[LinearCode, LinearCode]:
        """Return (rad(C), K) with C = rad(C) direct-sum K (Thm. 1.6).

        Both results are LinearCode subspaces with the same field and n.
        The splitting is orthogonal and K is nondegenerate. The radical
        is unique, but K is not: choose original RREF code basis vectors
        at the nonpivot columns of the RREF radical coefficient matrix.
        This makes the choice deterministic for the stored coordinates
        and field presentation, not intrinsic to the abstract space.

        A radical coefficient row has its own pivot coordinate; unit rows
        at the remaining coordinates complete it to a basis of F^size.
        Hence their corresponding code vectors complement the radical.
        Any such complement is nondegenerate: a vector of K orthogonal
        to K is also orthogonal to the radical, so lies in rad(C) and K,
        whose intersection is zero.

        The returned bases are ordinary RREF bases, not ordered normalized
        symplectic pairs. No codewords are enumerated.
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
        """Return rad(C) = C intersect C-perp as a subspace (Def. 1.5).

        This is the space of codewords orthogonal to all of C, not generally
        a largest isotropic subspace. It keeps the original field and n.
        """
        return self.symplectic_decomposition()[0]

    def symplectic_part(self) -> LinearCode:
        """Return the chosen nondegenerate complement K to rad(C).

        See symplectic_decomposition() for the deterministic choice and
        nonuniqueness. For an isotropic code K is the zero subspace.
        """
        return self.symplectic_decomposition()[1]

    @property
    def dimension(self) -> int:
        """Paper dimension dim(C): number of symplectic pairs (Def. 1.7).

        This equals (dim_F(C) - dim_F(rad(C)))/2, half the rank of the
        restricted alternating form. It is not dimension_over_field.
        """
        radical_dimension = self._radical_coefficients().shape[0]
        return (self.dimension_over_field - radical_dimension) // 2

    @property
    def isorank(self) -> int:
        """Maximum ordinary F-dimension of an isotropic subspace of C.

        Definition 1.7 gives irk(C) = dim(C) + dim_F(rad(C)), equivalently
        dim_F(C) - dim(C). This returns an integer, not a subspace.
        """
        return self.dimension_over_field - self.dimension

    @property
    def cardinality(self) -> int:
        """Number of codewords, |F|**dim_F(C), without enumerating them."""
        return self.field.order ** self.dimension_over_field

    def __contains__(self, vector: SymplecticVector) -> bool:
        """Test membership; reject wrong types, fields, or ambient lengths."""
        _require_code_vector(vector, self.field, self.n)
        augmented = _coefficient_matrix(self.basis + (vector,), self.field, self.n)
        return augmented.row_space().shape[0] == self.dimension_over_field

    def codewords(self) -> Iterator[SymplecticVector]:
        """Yield each codeword once, starting with zero.

        This is optional enumeration: construction only stores a basis.
        Exhausting this iterator takes |F|**dim_F(C) iterations. Scalars
        follow galois's integer encoding order, including extension fields.
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

    Example: linear_code([u, v, u + v], field=GF(2), n=2).
    The field and length are explicit even for an empty generating set.
    """
    return LinearCode(generators, field, n)