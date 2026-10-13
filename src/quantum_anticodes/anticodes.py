from collections.abc import Iterable, Iterator
from itertools import combinations

from galois import FieldArray

from ._coordinates import _coordinate_indices, _require_length
from .linear_code import LinearCode
from .symplectic import SymplecticVector, VElement, _require_field


def anticode(indices: Iterable[int], field: type[FieldArray], n: int) -> LinearCode:
    """Return A_J = {v in V^n : v_j = 0 for j outside J}
    """
    # ensure we have the appropriate field type, make sure its a valid subclass
    _require_field(field)
    # ensure the length is valid for the indices
    _require_length(n)
    # compute the python indices that we need for the anticode
    positions = _coordinate_indices(indices, n)

    # Generate the local basis (e and f) for the given field, the symplectic pair generators for each coord block
    local_basis = (VElement(1, 0, field=field), VElement(0, 1, field=field))

    # generate the basis vectors supported at the J components
    generators = (
        SymplecticVector(
            (element if i == position else 0 for i in range(n)), field=field,
        )
        # for each position, for each element
        for position in positions
        for element in local_basis
    )
    return LinearCode(generators, field, n)


# take in the field subclass and the lenght of the code then return a iterable structure of all possible anticodes for any subset
# return type represents a tuple of any length of integers (1,) or (1,2), etc
# we return the iterator of the J subset and the code ie ((1,2), code)
# Thus we can loop through these, and more memory efficient
def anticodes(field: type[FieldArray], n: int
              ) -> Iterator[tuple[tuple[int, ...], LinearCode]]:
    """Iterate over every (J, A_J), including empty and full J
    """
    # require the field to be an appropriate field, ie subsclass
    _require_field(field)
    # require lenght to be valid
    _require_length(n)
    # we first iterate through all possible subset size 0,...,n
    # Then combinations returns all possible subsets of size size and we return the indicy it self as well as the code generated.
    # for size = 3 we do subsets of (1,2,3),(1,3),(1,2),(2,3), (1),(2),(3),()
    return (
        (indices, anticode(indices, field, n))
        for size in range(n + 1)
        for indices in combinations(range(1, n + 1), size)
    )
