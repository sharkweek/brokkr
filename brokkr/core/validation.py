"""Validation tools for ``brokkr``."""

from brokkr.config import USYS
from brokkr.core.exceptions import (
    UnitDimensionError,
    BoundedValueError
)

def check_dimension(key, obj, dims, usys=USYS):
    """Validate dimensionality of an object."""

    # set dimensions for required objibutes
    if key in dims:
        correct_dim = dims.get(key)
        # check if object has units
        if not hasattr(obj, 'units'):
            obj *= usys[correct_dim[0]]  # assign first in tuple

        # if units assigned, check units for dimensionality
        else:
            if obj.units.dimensions not in correct_dim:
                raise UnitDimensionError(
                    key, ' OR '.join([str(i) for i in correct_dim])
                    )
            else:
                # units are converted to supplied unit system
                obj.convert_to_base(usys)

    else:
        raise KeyError("No dimensions found for `" + key + "` in `dims`.")

    return obj


def check_limit(key, obj, limits, usys=USYS):
    """Check that a value is within certain limits."""

    if key in limits:
        if out_of_bounds(obj, limits.get(key)):
            raise BoundedvalueError(key, limits.get(key))

    else:
        raise KeyError("No limits found for `" + key + "` in `limits`.")

def out_of_bounds(val, mn, mx, condition):
    """Check if value satisfies boundary conditions.

    Parameters
    ----------
    val : float or int
        value to evaluate against boundary conditions
    mn, mx: float or int
        the minimum and maximum boundaries
    condition: {'g', 'ge', 'g-l', 'ge-l', 'g-le', 'ge-le', 'l', 'le'}
        the boundary condition to evaluate

    Returns
    -------
    bool
        True if ``val`` is outside boundaries, False if within boundaries

    Notes
    -----
    ``condition`` should may be any of the values defined for each of the
    boundary conditions described in the table below:

    =========== ==================
    Value       Boundary Condition
    =========== ==================
    ``'g'``     ``val > mn``
    ``'ge'``    ``val >= mn``
    ``'g-l'``   ``mn < val < mx``
    ``'ge-l'``  ``mn <= val < mx``
    ``'g-le'``  ``mn < val <= mx``
    ``'ge-le'`` ``mn <= val <= mx``
    ``'l'``     ``val < mx``
    ``'le'``    ``val <= mx``
    =========== ==================

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
    }.get(condition)