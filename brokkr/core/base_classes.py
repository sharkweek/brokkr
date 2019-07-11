"""Base classes for ``brokkr``"""

from abc import ABCMeta, abstractmethod
from brokkr.config import USYS
from brokkr.core.exceptions import (
    UnitDimensionError,
    BoundedValueError
)
from brokkr.core.validation import out_of_bounds
from numpy import array, empty
from unyt.array import unyt_array


def abstract_attribute(obj=None):
    """Decorator for abstract attributes in ABCs."""

    if obj is None:
        obj = object()
    obj.__is_abstract_attribute__ = True
    return obj


class ExtendedABCMeta(ABCMeta):
    """Extended ABCMeta class that includes abstract attributes."""

    def __call__(cls, *args, **kwargs):
        instance = ABCMeta.__call__(cls, *args, **kwargs)
        abstract_attributes = {
            name
            for name in dir(instance)
            if getattr(getattr(instance, name),
                       '__is_abstract_attribute__', False)
        }

        if abstract_attributes:
            raise NotImplementedError(
                "Can't instantiate abstract class {} with"
                " abstract attributes: {}".format(
                    cls.__name__,
                    ', '.join(abstract_attributes)
                )
            )
        return instance


class DimensionedABC(metaclass=ExtendedABCMeta):
    """Abstract base class requiring units to be applied."""

    __slots__ = ['usys']
    _dimensions = {}
    _limits = {}

    @abstract_attribute
    def usys(self):
        pass

    def __setattr__(self, name, value):
        # set dimensions for required attributes
        if name in self._dimensions:
            dim = self._dimensions.get(name)

        # set dimensions for required attributes
        if name in dims:
            dim = self._dimensions.get(name)
            # check if value has units
            if not hasattr(value, 'units'):
                value *= usys[dim[0]]  # assign first in tuple

            # if units assigned, check units for dimensionality
            else:
                if value.units.dimensions not in dim:
                    raise UnitDimensionError(
                        name, ' OR '.join([str(i) for i in dim])
                        )
                else:
                    # units are converted to supplied unit system
                    value.convert_to_base(usys)

        # make sure value is within limits
        if name in self._limits:
            if out_of_bounds(value.value, **self._limits.get(name)):
                raise BoundedValueError(name, **self._limits.get(name))

        super().__setattr__(name, value)


class BaseTensor(unyt_array):
    """Base 3x3 symmetric tensor.

    Parameters
    ----------
    matrix : numpy.array_like
        an array-like object with six values
    units : string or unyt.Units
        the physical units for the tensor

    """

    def __new__(cls, matrix, units):
        """Create `BaseTensor` instance."""
        if array(matrix).size != 6:
            raise TypeError("`matrix` must have six values")
        else:
            new = unyt_array(empty((3, 3)), units).view(cls)
            return new

    def __init__(self, matrix):
        """Initialize `BaseTensor` instance."""
        for i, j in enumerate(array(matrix).ravel()):
            self.set_v_item(i, j)

    @property
    def voigt(self):
        """The Voigt representation of the array."""
        # create the indexing tuple for returning specific items in self
        voigt = [0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]

        return self[voigt].reshape((6, 1))

    def set_v_item(self, vindex, new_val):
        """Set an item in-place using Voigt index.

        Parameters
        ----------
        vindex : int
            the Voigt index of the item to set
        new_val : float
            the new value for the specified ``vindex``

        """

        i = [((0, 0),),
             ((1, 1),),
             ((2, 2),),
             ([1, 2], [2, 1]),
             ([2, 0], [0, 2]),
             ([0, 1], [1, 0])]

        for j, k in i[vindex]:
            self[j, k] = new_val

    def has_dimension(self, dim):
        """Check if dimensions of vector match `dim`.

        Parameters
        ----------
        dim : ``unyt.dimensions`` dimension
            the dimension to compare the instance against

        Returns
        -------
        bool
            True if instance's dimension matches ``dim``; False otherwise.

        """

        return self.units.dimensions == dim


class BaseVector(unyt_array):
    """Base 3x1 vector array.

    Parameters
    ----------
    matrix : numpy.array_like
        an array-like object with three items
    units : string or unyt.Unit object
        the physical units for the vector

    Notes
    -----
    If ``units`` is supplied as a string, it must conform to the same
    formatting rules as the ``unyt.array.unyt_array.to()`` method.

    """

    def __new__(cls, matrix, units):
        """Create `BaseVector` instance."""
        new = array(matrix)
        if new.size != 3:
            raise TypeError("`matrix` must have six values")
        else:
            vector = unyt_array(new, units).view(cls)
            vector.resize((3, 1))
            return vector

    @property
    def resultant(self):
        """The resultant magnitude."""

        return (self[0]**2 + self[1]**2 + self[2]**2)**(0.5)

    def has_dimension(self, dim):
        """Check if dimensions of vector match `dim`.

        Parameters
        ----------
        dim : ``unyt.dimensions`` dimension
            the dimension to compare the instance against

        Returns
        -------
        bool
            True if instance's dimension matches ``dim``; False otherwise.

        """

        return self.units.dimensions == dim