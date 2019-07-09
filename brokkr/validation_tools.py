"""Validation tools for ``brokkr``."""

from brokkr.config import USYS
from brokkr.exceptions import (
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

