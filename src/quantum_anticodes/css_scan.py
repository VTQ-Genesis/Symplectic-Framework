"""Test the existing trans_DJ condition on every CSS Choi subset."""

from itertools import combinations

from galois import FieldArray

from .css import choi
from .trans_dj import trans_DJ


# take in the field array of GF(2) of the x type and z type generator matrices
def scan_css_code(
    Gx: FieldArray,
    Gz: FieldArray,
    l: int = 2,
    *,
    parameters: dict | None = None,
) -> dict:

    # create a dict of parameters to pass in
    parameters = {} if parameters is None else dict(parameters)
    # grab the choi state
    Lx, Lz = choi(Gx, Gz)
    # n is number qubits described
    n = Gx.shape[1]
    # grab cols for choi state, larger by some k
    N = Lx.shape[1]

    passing_subsets = []
    # test range size from 0 to N (how large the subset will be)
    for size in range(N + 1):
        # for each size, get every subset of that size
        for J in combinations(range(1, N + 1), size):
            # test condition, append if passed
            if trans_DJ(Lx, J, l):
                passing_subsets.append(J)

    return {
        "n": n,
        "k": N - n,
        "l": l,
        "parameters": parameters,
        "passing_subsets": passing_subsets,
        "has_nonempty_pass": any(passing_subsets),
    }