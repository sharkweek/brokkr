"""Validation tools for ``brokkr``."""

from brokkr.config import DEFAULT_USYS
from brokkr.core.exceptions import (
    UnitDimensionError,
    BoundedValueError
)

__all__ = ['out_of_bounds']

def out_of_bounds(val, mn, mx, condition):
    """Check if value satisfies boundary conditions.

    Parameters
    ----------
    val : float or int
        Value to evaluate against boundary conditions
    mn, mx: float or int
        Minimum and maximum boundaries
    condition: {'g', 'ge', 'g-l', 'ge-l', 'g-le', 'ge-le', 'l', 'le'}
        Boundary condition to evaluate. Conditions are defined:

        =========== ===================
        Value       Boundary Condition
        =========== ===================
        ``'g'``     ``val > mn``
        ``'ge'``    ``val >= mn``
        ``'g-l'``   ``mn < val < mx``
        ``'ge-l'``  ``mn <= val < mx``
        ``'g-le'``  ``mn < val <= mx``
        ``'ge-le'`` ``mn <= val <= mx``
        ``'l'``     ``val < mx``
        ``'le'``    ``val <= mx``
        =========== ===================

    Returns
    -------
    bool
        True if ``val`` satisfies ``condition``; false otherwise

    """

    return not {
        'g': lambda: val > mn and mx is None,
        'ge': lambda: val >= mn and mx is None,
        'g-l': lambda: mn < val < mx,
        'ge-l': lambda: mn <= val < mx,
        'g-le': lambda: mn < val <= mx,
        'ge-le': lambda: mn <= val <= mx,
        'l': lambda: mn is None and val < mx,
        'le': lambda: mn is None and val <= mx
    }.get(condition)()
