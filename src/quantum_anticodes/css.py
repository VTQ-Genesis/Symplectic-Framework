"""Binary CSS stabilizer construction in the paper's symplectic coordinates.

The output S = e*C1 + f*C2 represents the stabilizer, not the normalizer
code C = S-perp used in Definition 2.3 and Example 2.5. The paper supplies
the symplectic framework; it does not explicitly define these CSS routines.
"""
from galois import GF, FieldArray

from .linear_code import LinearCode, _coefficient_matrix
from .symplectic import SymplecticVector, VElement


_BINARY_FIELD = GF(2)


# take in input as field arrays of binary field F([[1,0,1],[0,1,1]])
def css_construct(G1: FieldArray, G2: FieldArray) -> LinearCode:
    """Return the CSS stabilizer subspace S = e*C1 direct-sum f*C2.

    G1 and G2 must be two-dimensional galois.GF(2) arrays whose rows
    generate classical binary codes of the same length n. Require
    G1 @ G2.T = 0 over GF(2), equivalently C1 <= C2-perp for the ordinary
    binary dot product. Lists, ordinary NumPy arrays, SymplecticVector
    objects, and arrays over other fields are not implicitly converted.

    Zero, duplicate, and dependent rows are accepted. To specify a zero
    classical code use GF(2).Zeros((0, n)); n=0 is supported. Inputs are
    unchanged. The returned LinearCode has length n and ordinary dimension
    rank(G1) + rank(G2). It is isotropic, so its paper dimension is zero;
    this is not the number of encoded logical qubits of the quantum code.
    """
    # ensure we have the matrix is only of binary and correct dimension (matrix)
    for matrix in (G1, G2):
        if type(matrix) is not _BINARY_FIELD or matrix.ndim != 2:
            raise TypeError("generators must be two-dimensional GF(2) arrays")
    # ensure that the codes are the same length so same number of columns
    if G1.shape[1] != G2.shape[1]:
        raise ValueError("both classical codes must have the same length")
    # ensure that G1 contained in G2.T, ie every vector of G1 orthog to G2.T, .any returns
    # true if any entry is nonzero
    if (G1 @ G2.T != 0).any():
        raise ValueError("CSS requires C1 to be contained in C2-perp")

    # grab the code length
    n = G1.shape[1]

    # Now go through all of the rows of G1 and convert each element to e -> (1,0)
    e_generators = [
        SymplecticVector(
            (VElement(int(a), 0, field=_BINARY_FIELD) for a in row),
            field=_BINARY_FIELD,
        )
        for row in G1
    ]
    # go through rows of G2 and do the same turning into fs
    f_generators = [
        SymplecticVector(
            (VElement(0, int(b), field=_BINARY_FIELD) for b in row),
            field=_BINARY_FIELD,
        )
        for row in G2
    ]
    # return the direct sum output, so just concatenate the generators together, as they are lists
    return LinearCode(e_generators + f_generators, _BINARY_FIELD, n)


# take in a LinearCode in symplectic formalism, return two classical codes such that
# L = eC1 + fC2 such that C1 in C2.T
def reverse_css(code: LinearCode) -> tuple[FieldArray, FieldArray]:
    """Recover binary generator matrices (G1, G2) from a CSS stabilizer S.

    The input must be a LinearCode over GF(2), isotropic and equal to the
    sum of its pure-e and pure-f subspaces in the fixed coordinates.
    Reject noncommuting codes and codes that do not split this way with
    ValueError. No change of local e,f basis is attempted.

    Return independent RREF generators for C1 = {a : e*a in S} and
    C2 = {b : f*b in S}. Thus css_construct(*reverse_css(S)) == S, but
    original generator matrices need not be recovered row-for-row.
    A zero sector returns shape (0, S.n). The input is unchanged; returned
    arrays are independently allocated and may be edited by the caller.
    Codewords are not enumerated.
    """
    # ensure we have appropriate input, so much be over binary
    if not isinstance(code, LinearCode):
        raise TypeError("expected a LinearCode")
    if code.field is not _BINARY_FIELD:
        raise ValueError("expected a code over GF(2)")

    # store the inputs by there coefficents
    matrix = _coefficient_matrix(code.basis, code.field, code.n)

    # A takes e's we start at 0 and take every other componet, B takes f components
    A = matrix[:, 0::2]
    B = matrix[:, 1::2]

    # Code word has form (xA, xB), to get the pure e components, we need xB to be 0
    # B.leftnull gives all such x with xB = 0
    # multiply this by A to get the pure e coefficients, then return the independent row space
    G1 = (B.left_null_space() @ A).row_space()
    G2 = (A.left_null_space() @ B).row_space()

    # Validate orthogonality as well as spanning the entire input subspace.
    reconstructed = css_construct(G1, G2)
    if reconstructed != code:
        raise ValueError("the code does not split into pure-e and pure-f parts")
    return G1, G2