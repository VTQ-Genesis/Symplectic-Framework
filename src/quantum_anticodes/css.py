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
def choice(Gx: FieldArray, Gz: FieldArray) -> tuple[FieldArray, FieldArray]:
    """Choose classical complements and append identity-coordinate blocks.
    """
    # Ensure that the given is of the correct type, a matrix of F = 2
    for matrix in (Gx, Gz):
        if type(matrix) is not _BINARY_FIELD or matrix.ndim != 2:
            raise TypeError("generators must be two-dimensional GF(2) arrays")
    # ensure the length (number of columns is the same)
    if Gx.shape[1] != Gz.shape[1]:
        raise ValueError("both classical codes must have the same length")


    # Grab the codes cz and cx
    cz = Gz.row_space()
    cx = Gx.row_space()



    # grab dx = cz perp and same for dz
    dx = cz.null_space()
    dz = cx.null_space()
   
    # compute the basis for these 
    hx = dx.row_space()
    hz = dz.row_space()
   

    
    # Choose rows in Cz that add new independent directions to Dz
    bz_rows = []
    # span starts as all of dz
    span = dz
    # loop through basis in cz looking if it adds independent direction
    for v in cz:
        # stack v in cz under span dz and generate basis
        extended = np.vstack((span, v)).row_space()
        # check if adding this vector increased the # rows, ie new direction
        if extended.shape[0] > span.shape[0]:
            bz_rows.append(v)
            span = extended
    # if non empty rows store the matrix of the rows, else make it 0 code, so no error      
    bz = (np.vstack(bz_rows) if bz_rows
          else _BINARY_FIELD.Zeros((0, cz.shape[1])))

    # Choose rows in Cx that add new independent directions to Dx
    bx_rows = []
    span = dx
    for v in cx:
        extended = np.vstack((span, v)).row_space()
        if extended.shape[0] > span.shape[0]:
            bx_rows.append(v)
            span = extended
    bx = (np.vstack(bx_rows) if bx_rows
          else _BINARY_FIELD.Zeros((0, cx.shape[1])))

    # create final matrices for the block
    ax = bx.row_space()
    az = bz.row_space()
  
    # grab num rows for ax and az
    kx = ax.shape[0]
    kz = az.shape[0]

    # Construct the Lx and Lz matrices
    Lx = _BINARY_FIELD(np.block([
        [hx, _BINARY_FIELD.Zeros((hx.shape[0], kx))],
        [ax, _BINARY_FIELD.Identity(kx)],
    ]))
    Lz = _BINARY_FIELD(np.block([
        [hz, _BINARY_FIELD.Zeros((hz.shape[0], kz))],
        [az, _BINARY_FIELD.Identity(kz)],
    ]))
    return Lx.row_space(), Lz.row_space()
