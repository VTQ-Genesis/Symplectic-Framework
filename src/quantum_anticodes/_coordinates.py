"""Validation shared by the coordinate operations of Section 2."""

from collections.abc import Iterable


def _require_length(n: int) -> None:
    if type(n) is not int:
        raise TypeError("n must be a Python integer")
    if n < 0:
        raise ValueError("n must be non-negative")


def _coordinate_indices(indices: Iterable[int], n: int) -> tuple[int, ...]:
    """Validate paper labels 1,...,n and return sorted zero-based positions.

    A subset may be supplied in any order, but repetitions are rejected
    rather than silently discarded. Consume an input iterator only once.
    """
    labels = tuple(indices)
    for label in labels:
        if type(label) is not int:
            raise TypeError("coordinate indices must be Python integers")
        if not 1 <= label <= n:
            raise ValueError(f"coordinate indices must be in 1,...,{n}")
    if len(set(labels)) != len(labels):
        raise ValueError("coordinate indices must be distinct")
        # returns the 1 label position, so important when using to do labeling in terms of starting at 1
    return tuple(label - 1 for label in sorted(labels))
