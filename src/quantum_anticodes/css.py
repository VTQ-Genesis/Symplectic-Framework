from galois import GF, FieldArray
import numpy as np
from .linear_code import LinearCode, _coefficient_matrix
from .symplectic import SymplecticVector, VElement


_BINARY_FIELD = GF(2)


# take in input as field arrays of binary field F([[1,0,1],[0,1,1]])
def css_construct(G1: FieldArray, G2: FieldArray) -> LinearCode:
    """Return the CSS stabilizer subspace S = e*C1 direct-sum f*C2.
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


# take in GF(2)([[]]) matrix objects and return the tuple of Lx and Lz new codes
# take in the generator matrices
# take in classical normalizer generators and return the Choi X/Z generators
def choi(Gx: FieldArray, Gz: FieldArray) -> tuple[FieldArray, FieldArray]:
    """Return generator matrices for a binary CSS encoder's Choi state.
    """
    # First we must ensure that the matrices are over a binary field and are actually
    # matrices
    for matrix in (Gx, Gz):
        if type(matrix) is not _BINARY_FIELD or matrix.ndim != 2:
            raise TypeError("generators must be two-dimensional GF(2) arrays")
    # ensure the length is the same, num of cols
    if Gx.shape[1] != Gz.shape[1]:
        raise ValueError("both classical codes must have the same length")


    # Compute the row space of the generators, and then grab the length
    cx = Gx.row_space()
    cz = Gz.row_space()
    n = cx.shape[1]

    # We set hx to cz perp, and it follows that we need cz perp to be contained in
    # cx, if not its unvalid.
    # Thus check if cz perp adds any additional directions, if yes not contained
    hx = cz.null_space()
    combined = np.vstack((cx, hx)).row_space()
    if combined.shape[0] != cx.shape[0]:
        raise ValueError("Choi construction requires Cz-perp to be contained in Cx")

    # Keep each row in cx that adds a independent direction to span which is contained
    # of hx = cz perp, and the other independent directions cx
    ax_rows = []
    span = hx
    for v in cx:
        extended = np.vstack((span, v)).row_space()
        if extended.shape[0] > span.shape[0]:
            ax_rows.append(v)
            span = extended


    # create ax, then if empty, just make the empty code
    ax = (np.vstack(ax_rows) if ax_rows
          else _BINARY_FIELD.Zeros((0, n)))

    # grab the number of rows of the ax matrix
    k = ax.shape[0]

    # Compute the Lx maually as the block matrix
    Lx = _BINARY_FIELD(np.block([
        [hx, _BINARY_FIELD.Zeros((hx.shape[0], k))],
        [ax, _BINARY_FIELD.Identity(k)],
    ])).row_space()

    # Lz must be Lx perp
    Lz = Lx.null_space()
    
    return Lx, Lz
