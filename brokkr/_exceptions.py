"""Exceptions for ``brokkr``"""

from .mech_math import out_of_bounds


class UnitDimensionError(Exception):
    """Error raised if wrong unit type is assigned.

    Parameters
    ----------
    name : str
        name of entity to print in error message
    correct_dim : unyt.dimensions.dimension
        the unit that should have been assigned

    """

    def __init__(self, name, correct_dim):
        self.name = name
        self.correct_dim = correct_dim

    def __str__(self):
        return (
            f"`{self.name}` must have units with dimensions {self.correct_dim}"
        )


class BoundedValueError(Exception):
    """Error for values that are out of specified boundaries.

    Parameters
    ----------
    name : string
        name of object checked for boundary
    mn : float
        the minimum value checked against
    mx : float
        the maximum value checked against
    condition: {'g', 'ge', 'g-l', 'ge-l', 'g-le', 'ge-le', 'l', 'le'}
        the condition ``name`` was evaluated against ;see 'Notes' for
        acceptable values

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

    def __init__(self, name, mn, mx, condition):
        self.name = name
        self.both = True

        if mn is not None:
            self.mn = str(mn) + ' '
        else:
            self.mn = ''
            self.both = False

        if mx is not None:
            self.mx = str(mx) + ' '
        else:
            self.mx = ''
            self.both = False

        args = {
            'g': {
                'mn': str(mn) + ' ',
                'mn_eq': '>',
                'mx': '',
                'mx_eq': '',
                'nd': ''
            },
            'ge': {
                'mn': str(mn) + ' ',
                'mn_eq': '>=',
                'mx': '',
                'mx_eq': '',
                'nd': ''
            },
            'g-l': {
                'mn': str(mn) + ' ',
                'mn_eq': '>',
                'mx': str(mx) + ' ',
                'mx_eq': '<',
                'nd': 'and '
            },
            'ge-l': {
                'mn': str(mn) + ' ',
                'mn_eq': '>=',
                'mx': str(mx) + ' ',
                'mx_eq': '<',
                'nd': 'and '
            },
            'g-le': {
                'mn': str(mn) + ' ',
                'mn_eq': '>',
                'mx': str(mx) + ' ',
                'mx_eq': '<=',
                'nd': 'and '
            },
            'ge-le': {
                'mn': str(mn) + ' ',
                'mn_eq': '>=',
                'mx': str(mx) + ' ',
                'mx_eq': '<= ',
                'nd': 'and '
            },
            'l': {
                'mn': '',
                'mn_eq': '',
                'mx': str(mx) + ' ',
                'mx_eq': '<',
                'nd': ''
            },
            'le': {
                'mn': '',
                'mn_eq': '',
                'mx': str(mx) + ' ',
                'mx_eq': '<=',
                'nd': ''
            }
        }.get(condition)

        for each in args:
            self.__setattr__(each, args[each])


    def __str__(self):
        return (
            f"`{self.name}` must be a number {self.mn_eq}{self.mn}{self.nd}"
            + f"{self.mx_eq}{self.mx}"
        )


class CoefficientError(BoundedValueError):
    """Error for values that must be between 0 and 1.

    .. warning:: ``CoefficientError`` will accept all ``condition`` values
        that ``BoundedValueError`` accepts, but should limited to ``3``
        through ``6``. Use ``BoundedValueError`` for any other conditions.

    """

    def __init__(cls, name, condition):
        super().__init__(name, 0, 1, condition)
