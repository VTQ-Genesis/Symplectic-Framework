from galois import FieldArray

from .symplectic import _require_same_field


def schur_product_vectors(u: FieldArray, v: FieldArray) -> FieldArray:
    # Ensure vectors are of valid types as well as a vector so dim is 1
    for vector in (u, v):
        if not isinstance(vector, FieldArray) or vector.ndim != 1:
            raise TypeError("vectors must be one-dimensional finite-field arrays")
    # require the same field to be used
    _require_same_field(type(u), type(v))
    if u.shape != v.shape:
        raise ValueError("both vectors must have the same length")

    # galois coord mult
    return u * v

# take in generators for the codes and return the span of their schur product
def schur_product_codes(G1: FieldArray, G2: FieldArray) -> FieldArray:
    # Ensure the matrices are approproiate type and are in fact matrices
    for matrix in (G1, G2):
        if not isinstance(matrix, FieldArray) or matrix.ndim != 2:
            raise TypeError("generators must be two-dimensional finite-field arrays")
    # require the same field 
    _require_same_field(type(G1), type(G2))
    
    # check they both have same lenght (# columns)
    if G1.shape[1] != G2.shape[1]:
        raise ValueError("both classical codes must have the same length")

    # compute basis for both, so remove the redundant rows
    basis1 = G1.row_space()
    basis2 = G2.row_space()

    # store the field
    field = type(G1)
    # grab the length
    n = G1.shape[1]
    # create the 0 matrix in the appropriate field of product of basis vectors and appropriate length
    # each pairwise schur product becomes a row, then take basis
    products = field.Zeros((basis1.shape[0] * basis2.shape[0], n))

    # current row
    row = 0
    for u in basis1:
        for v in basis2:
            # fill up row with schur product
            products[row] = schur_product_vectors(u, v)
            row += 1

    # Take span and return
    return products.row_space()