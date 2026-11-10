import numpy as np

from quantum_anticodes import choi, schur_product_codes


def scan_css_code_exact(Gx, Gz, l=2, *, parameters=None):

    # basic type checking
    if type(l) is not int:
        raise TypeError("l must be a Python integer")
    if l < 1:
        raise ValueError("l must be positive")

    # grab the parameters and compute choi state for the code
    parameters = {} if parameters is None else dict(parameters)
    Lx, Lz = choi(Gx, Gz)
    n = Gx.shape[1]
    N = Lx.shape[1]

    # Keep in mind the condition is punc_j(C)<= short_j((C^*l).T)

    # The extra factor comes from the code being punctured. Thus l=1
    # needs the square and l=2 needs the cube, including repeated factors.
    power = Lx
    for _ in range(l):
        power = schur_product_codes(power, Lx)
    indicator_basis = power.null_space()

    # Every binary linear combination, including the zero indicator.
    # At each step, keep existing words and add their sums with the new row.
    indicators = type(Lx).Zeros((1, N))
    for row in indicator_basis:
        indicators = np.vstack((indicators, indicators + row))

    passing_subsets = [
        tuple((np.flatnonzero(indicator) + 1).tolist())
        for indicator in indicators
    ]
    passing_subsets.sort(key=lambda J: (len(J), J))

    return {
        "n": n,
        "k": N - n,
        "l": l,
        "parameters": parameters,
        "passing_subsets": passing_subsets,
        "has_nonempty_pass": indicator_basis.shape[0] > 0,
    }
