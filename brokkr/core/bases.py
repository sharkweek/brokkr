"""Base classes for ``brokkr``"""

# standard imports
from abc import ABCMeta, abstractmethod

# third party imports
from numpy import array, empty
from unyt import UnitSystem, unyt_array

# local imports
from brokkr.config import DEFAULT_USYS
from brokkr.core.decorators import abstract_attribute
from brokkr.core.exceptions import (
    BoundedValueError,
    DerivedAttributeError,
    UnitDimensionError,
)
from brokkr.core.validation import out_of_bounds

__all__ = ['abstract_attribute', 'DimensionedABC', 'CalculatedABC',
           'BaseTensor', 'BaseVector']


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
    """ABC that requires dimensions and limits for specified attributes.

    Dimensions for subclass attributes are limited by the dimensionality (i.e.
    length, time, force, etc.) definitions provided by the :attr:`_dimensions`
    attribute. Value limits (e.g. >0 or >=1, etc.) can also be defined by use of
    the :attr:`_limits` attribute.

    Parameters
    ----------
    usys : unyt.Unit_System
        the unit system in which to evaluate dimensions

        .. note::
            If an attribute is given a value without units, the default unit
            for the dimensions of that attribute defined in :attr:`_limits` is
            assigned automatically.

    Attributes
    ----------
    _dimensions : (class) dict of {str: (unyt.dimension, )}
        dimensions for each dimensioned attributed. (For dimensions types see
        :py:mod:`unyt.dimensions`.)
    _limits : (class) dict of {str: {'mn': float or int, 'mx': float or int,
        'condition': str}}
        limits for values that are bounded

        Each key in ``_limits`` represents an attribute bounded by limits. The
        corresponding dictionary defines the types of limits. See documentation
        for ``out_of_bounds`` in the ``validation`` module for definitions of
        ``mn``, ``mx``, and ``condition``.

    Example
    -------
    ::

        class Displacement(DimensionedABC):
            _dimensions = {'x': (unyt.dimensions.length, ),
                           't': (unyt.dimensions.time, )}
            _limits = {'x': {'mn': 0, 'mx': None, 'condition': 'g'},
                       't': {'mn': 10, 'mx': 20, 'condition': 'ge-le'}}

            def __init__(self, x, t, usys=DEFAULT_USYS):
                self.usys = usys
                self.x = x
                self.t = t

        >>> u = Displacement(5, 15)
        >>> u.x
        unyt_quantity(5, 'inch')
        >>> u.t
        unyt_quantity(15, 's')

        >>> u.x = 1 * unyt.Unit('psi')
        UnitDimensionError: `x` must have units with dimensions (length)

        >>> u.t = 1
        BoundedValueError: `t` must be a number >=10 and <=20

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
            if out_of_bounds(val=value.value, **self._limits.get(name)):
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


class CalculatedABC(metaclass=BrokkrABCMeta):
    """ABC for classes with derived attributes.

    For objects possessing attributes that are interdependent, this ABC provides
    functionality to 'lock' certain derived attributes. This is useful for
    situations in which adding many ``@property`` decorators may be cumbersome.
    It also may provide quicker processing time due to values being stored in
    instance attributes rather than having to be recalculated each time the
    property is called.

    Attributes
    ----------
    _base_attr : list or tuple of str
        base attributes
    _calc_attr : list or tuple of str
        derived or 'calculated' attributes
    _free_attr : list or tuple of str
        free attributes that are independent of any other attribute

    Methods
    -------
    _update()
        updates all derived attributes based on values of base attributes

    Notes
    -----
    The :func:`~brokkr.core.decorators.unlock` decorator must be used with
    with the :meth:`~brokkr.core.CalculatedABC._update` method.

    It is generally recommended that subclasses created from
    :class:`~brokkr.core.bases.CalculatedABC` contain a :attr:`__slots__`
    attribute, given that the class is designed for computing speed.

    Example
    -------
    An object defining a circle would have a base attribute for ``radius`` and
    a derived attribute for ``area``::

        class Circle(CalculatedABC):
            _base_attr = ('radius')
            _calc_attr = ('area')
            _free_attr = ('__locked')

            __slots__ = ('radius', 'area', '__locked')

            def __init__(self, radius):
                self.radius = radius

                self._update()

            @unlock
            def _update(self):
                self.area = 3.14 * self.radius ** 2

        >>> x = Circle(5)
        >>> x.radius
        5
        >>> x.area
        25

        >>> x.area = 4
        DerivedAttributeError: `area` is a derived attribute and cannot be
        set manually.

    """

    _base_attr = (,)
    _calc_attr = (,)
    _free_attr = (,)

    __locked = abstract_attribute()

    def __setattr__(self, name, value):
        """Extended to protect calculated attributes."""

        if self.__locked:
            # udpate ply and laminate after updated properties are set
            if name in self._base_attr:
                super().__setattr__(name, value)
                self._update()

            # don't set protected values
            elif name in self._calc_attr:
                raise DerivedAttributeError(name)

        else:
            super().__setattr__(name, value)

    @abstractmethod
    def _update(self):
        """Updates all derived attributes based on the base attributes."""
        pass


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
