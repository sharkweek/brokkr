"""Decorators for ``brokkr``."""

from unyt import UnitSystem

from brokkr.core.exceptions import (
    BoundedValueError,
    UnitDimensionError,
)
from brokkr.core.validation import out_of_bounds

def abstract_attribute(obj=None):
    """Decorator for abstract attributes.

    It creates abstract attributes that must be defined in all classes that are
    subclassed from ABCs that use the ``BrokkrABCMeta`` metaclass.

    Parameters
    ----------
    obj : var or function, optional
        the variable to assign an abstract attribute to or the function being
        decorated to behave as an abstract attribute.

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
        ... 'hello'

        >>> y = BadFoo('world')
        ... NotImplementedError: Can't instantiate abstract class BadFoo with
            abstract attributes: hello

    Alternatively, attributes can be declared as class variables::

        class MyABC(metaclass=BrokkrABCMeta):
            hello = abstract_attribute()

    This should yield the same results as using the decorator for an empty
    class method as shown above.

    """

    # assign __is_abstract_attribute__ flag
    if obj is None:

        class Obj:
            pass  # dummy class for attribute assignment

        obj = Obj()

    obj.__is_abstract_attribute__ = True

    return obj

def unlock(locked_method):
    """Decorate class methods to unlock derived attributes.

    Unlocks derived attirbutes of an object to allow ``locked_method`` to
    operate freely by setting the ``__locked`` attribute to ``False``. Once
    ``locked_method`` is finished executing, ``__locked`` is set to True,
    re-locking the derived attributes.

    Parameters
    ----------
    locked_method : method
        The function to unlock.

    Returns
    -------
    method
        An unlocked function.

    """
    def unlocked_method(self, *args, **kwargs):
        cls = self.__class__
        cls_name = cls.__name__

        super(cls, self).__setattr__('_' + cls_name + '__locked', False)

        locked_method(self, *args, **kwargs)

        super(cls, self).__setattr__('_' + cls_name + '__locked', True)

    return unlocked_method


def validate_attr(method):
    """Attribute validation decorator.

    Validates an attribute before executing the decorated method.

    """

    def validated_setattr(self, name, value):
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
                        name, ' OR '.join([str(i) for i in dim]))
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
                    "(e.g. R or K)")

            # convert only for attributes that exist in the instance
            for each in {x: self._dimensions[x]
                         for x in self._dimensions if hasattr(self, x)}:
                getattr(self, each).convert_to_base(value)

        method(self, name, value)

    return validated_setattr