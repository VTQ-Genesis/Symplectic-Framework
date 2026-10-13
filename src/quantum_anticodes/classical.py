from collections.abc import Iterable

from galois import FieldArray

from ._coordinates import _coordinate_indices


def puncture_classical(G: FieldArray, indices: Iterable[int]) -> FieldArray:
    # ensure that our array is valid
    if not isinstance(G, FieldArray) or G.ndim != 2:
        raise TypeError("generators must be a two-dimensional finite-field array")
    # ensure the cooridnates are valid
    positions = _coordinate_indices(indices, G.shape[1])

    # Select all of the rows and the supported cols then get basis
    return G[:, positions].row_space()


def shorten_classical(G: FieldArray, indices: Iterable[int]) -> FieldArray:
    # ensure its a valid matrix 
    if not isinstance(G, FieldArray) or G.ndim != 2:
        raise TypeError("generators must be a two-dimensional finite-field array")

    # ensure the postions are appropriate for pythion indexing
    positions = _coordinate_indices(indices, G.shape[1])

    # form a set for effcient comparison
    selected = set(positions)
    #  loop through the columns, then choose the is to be the indices not in j
    outside = [i for i in range(G.shape[1]) if i not in selected]

    # compute the basis
    basis = G.row_space()
    # select all of the rows,the the matrix is the columns not selected
    # left null here indicates the the combos such that the outside rows combine to 0
    coefficients = basis[:, outside].left_null_space()
    return (coefficients @ basis[:, positions]).row_space()
