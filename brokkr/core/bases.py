"""Base classes for ``brokkr``"""

# standard imports
from abc import ABCMeta, abstractmethod

# third party imports
from numpy import array, empty
from unyt import UnitSystem, unyt_array

# local imports
from brokkr.config import DEFAULT_USYS
from brokkr.core.exceptions import (
    UnitDimensionError,
    BoundedValueError
)
from brokkr.core.validation import out_of_bounds

__all__ = ['abstract_attribute', 'DimensionedABC', 'BaseTensor', 'BaseVector']

def abstract_attribute(obj=None):
    """Decorator for abstract attributes.

    It creates abstract attributes that must be defined in all classes that are
    subclassed from ABCs that use the ``BrokkrABCMeta`` metaclass.

    Parameters
    ----------
    obj : optional

    Example
    -------
    To use this decorator, subclass attributes should be defined as functions
    in the ABC. For example::

        class MyABC(metaclass=BrokkrABCMeta):
            @abstract_attribute
            def hello(self):
                pass

        class GoodFoo(MyABC):
            def __init__(self, hello, world):
                self.hello = hello
                self.world = world

        class BadFoo(MyABC):
            def __init__(self, world):
                self.world = world

        >>> x = GoodFoo('hello', 'world')
        >>> x.hello
        'hello'

        >>> y = BadFoo('world')
        NotImplementedError: Can't instantiate abstract class BadFoo with
        abstract attributes: hello

    Alternatively, attributes can be declared as class variables::

        class MyABC(metaclass=BrokkrABCMeta):
            hello = abstract_attribute()

    This should yield the same results as using the decorator for an empty
    class method as shown above.

    """

    # assign __is_abstract_attribute__ flag
    if obj is None:
        class Obj: pass  # dummy class for attribute assignment
        obj = Obj()

    obj.__is_abstract_attribute__ = True

    return obj


class BrokkrABCMeta(ABCMeta):
    """Extended ABCMeta class including provisions for abstract attributes."""

    def __call__(cls, *args, **kwargs):
        instance = ABCMeta.__call__(cls, *args, **kwargs)

        # create dictionary of abstract attributes not declared in subclass
        abstract_attributes = {
            name
            for name in dir(instance)
            if getattr(getattr(instance, name),
                       '__is_abstract_attribute__', False)
        }

        # throw error if any abstract attributes are not declared
        if abstract_attributes:
            raise NotImplementedError(
                "Can't instantiate abstract class {} with"
                " abstract attributes: {}".format(
                    cls.__name__,
                    ', '.join(abstract_attributes)
                )
            )

        return instance


class DimensionedABC(metaclass=BrokkrABCMeta):
    """ABC that requires dimensions and limits to specified attributes.

    Attributes
    ----------
    _dimensions : (class) dict of {str: (unyt.dimension,)}
        dimensions for each dimensioned attributed
    _limits : (class) dict of {str: {'mn': float or int, 'mx': float or int,
        'condition': str}}
        limits for values that are bounded

        Each key in ``_limits`` represents an attribute bounded by limits. The
        corresponding dictionary defines the types of limits. See documentation
        for ``out_of_bounds`` in the ``validation`` module for definitions of
        ``mn``, ``mx``, and ``condition``.
    usys : unyt.UnitSystem
        instance's unit system

    """

    _dimensions = {}
    _limits = {}

    usys = abstract_attribute()

    def __setattr__(self, name, value):
        """Extended to validate dimensioned and limited attributes."""

        # set dimensions for required attributes
        if name in self._dimensions:
            dim = self._dimensions.get(name)

            # check if value has units
            if not hasattr(value, 'units'):
                value *= self.usys[dim[0]]  # assign first in tuple

            # if units are assigned, check dimensionality
            else:
                if value.units.dimensions not in dim:
                    raise UnitDimensionError(
                        name, ' OR '.join([str(i) for i in dim])
                        )
                else:
                    # units are converted to supplied unit system
                    value.convert_to_base(self.usys)

        # make sure value is within limits
        if name in self._limits:
            if out_of_bounds(value.value, **self._limits.get(name)):
                raise BoundedValueError(name, **self._limits.get(name))

        # catch unit system change and convert all attributes with units
        if name == 'usys':
            # validate unit system
            if type(value) != UnitSystem:
                raise AttributeError(f"`{name}` must be a `UnitSystem`")
            elif value['temperature'].base_offset != 0:
                raise AttributeError(
                    "Unit system must have an absolute temperature unit",
                    "(e.g. R or K)"
                )

            # convert only for attributes that exist in the instance
            for each in {x: self._dimensions[x]
                         for x in self._dimensions
                         if hasattr(self, x)}:
                getattr(self, each).convert_to_base(value)

        super().__setattr__(name, value)


class BaseTensor(unyt_array):
    """Base 3x3 symmetric tensor.

    Parameters
    ----------
    matrix : array_like
        an ``array_like`` object with six values
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
        dim : dimension from ``unyt.dimensions``
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
    matrix : array_like
        an ``array_like`` object with three items
    units : str or unyt.Unit
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
            raise TypeError("`matrix` must have 3 values")
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
        dim : dimension from ``unyt.dimensions``
            the dimension to compare the instance against

        Returns
        -------
        bool
            True if instance's dimension matches ``dim``; False otherwise.

        """

        return self.units.dimensions == dim
