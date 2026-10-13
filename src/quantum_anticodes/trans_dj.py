from collections.abc import Iterable

import numpy as np
from galois import FieldArray

from .classical import puncture_classical, shorten_classical
from .schur import schur_product_codes

# take in classical code, subset to puncture/shorten on then a int to iterate and return bool for inclusion param
def trans_DJ(C1: FieldArray, J: Iterable[int], l: int) -> bool:
    # ensure type checks are valid
    if type(l) is not int:
        raise TypeError("l must be a Python integer")
    if l < 1:
        raise ValueError("l must be positive")

    # return tuple for the coords, _coordinate_indices called through the puncture
    J = tuple(J)
    punctured = puncture_classical(C1, J)

    # basis for the code, then 
    basis = C1.row_space()
    D1 = basis
    # no loop variable required, then compute the l schur product of C1
    for _ in range(1, l):
        D1 = schur_product_codes(D1, basis)

    # compute the dual for D1 and shorten
    dual = D1.null_space()
    shortened = shorten_classical(dual, J)

    # Containment condition is simply if punctured <= shortened
    # if punctured adds addtional dimensions cant be contained, so check the basis
    # for the matrix
    combined = np.vstack((shortened, punctured)).row_space()
    return bool(combined.shape[0] == shortened.shape[0])
