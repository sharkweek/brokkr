"""Exceptions for ``brokkr``"""


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
    condition: int
        the condition ``name`` was evaluated against ;see 'Notes' for
        acceptable values

    Notes
    -----
    ``condition`` should may be any of the values defined for each of the
    boundary conditions described in the table below:

    ===== ==================
    Value Boundary Condition
    ===== ==================
    ``1`` ``val > mn``
    ``2`` ``val >= mn``
    ``3`` ``mn < val < mx``
    ``4`` ``mn <= val < mx``
    ``5`` ``mn < val <= mx``
    ``6`` ``mn <= val <= mx``
    ``7`` ``val < mx``
    ``8`` ``val <= mx``
    ===== ==================

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
            1: {'mn': str(mn) + ' ',
                'mn_eq': '> ',
                'mx': '',
                'mx_eq': '',
                'nd': ''},
            2: {'mn': str(mn) + ' ',
                'mn_eq': '>= ',
                'mx': '',
                'mx_eq': '',
                'nd': ''},
            3: {'mn': str(mn) + ' ',
                'mn_eq': '< ',
                'mx': str(mx) + ' ',
                'mx_eq': '< ',
                'nd': 'and '},
            4: {'mn': str(mn) + ' ',
                'mn_eq': '<= ',
                'mx': str(mx) + ' ',
                'mx_eq': '< ',
                'nd': 'and '},
            5: {'mn': str(mn) + ' ',
                'mn_eq': '< ',
                'mx': str(mx) + ' ',
                'mx_eq': '<= ',
                'nd': 'and '},
            6: {'mn': str(mn) + ' ',
                'mn_eq': '<= ',
                'mx': str(mx) + ' ',
                'mx_eq': '<= ',
                'nd': 'and '},
            7: {'mn': '',
                'mn_eq': '',
                'mx': str(mx) + ' ',
                'mx_eq': '< ',
                'nd': ''},
            8: {'mn': '',
                'mn_eq': '',
                'mx': str(mx) + ' ',
                'mx_eq': '<= ',
                'nd': ''}
        }.get(condition)

        for each in args:
            self.__setattr__(each, args[each])


    def __str__(self):
        return (
            f"`{self.name}` must be a number {self.mn_eq}{self.mn}{self.nd}"
            + f"{self.mx_eq}{self.mx}"
        )


class CoefficientError(BoundedValueError):
    """Error for values that must be between 0 and 1."""

    def __init__(cls, name, condition):
        super().__init__(name, 0, 1, condition)


def check_bounds(val, mn, mx, condition):
    """Check if value satisfies boundary conditions.

    Parameters
    ----------
    val : float or int
        value to evaluate against boundary conditions
    mn, mx: float or int
        the minimum and maximum boundaries
    condition: int {1 through 8}
        the boundary condition to evaluate

    Notes
    -----
    ``condition`` should may be any of the values defined for each of the
    boundary conditions described in the table below:

    ===== ==================
    Value Boundary Condition
    ===== ==================
    ``1`` ``val > mn``
    ``2`` ``val >= mn``
    ``3`` ``mn < val < mx``
    ``4`` ``mn <= val < mx``
    ``5`` ``mn < val <= mx``
    ``6`` ``mn <= val <= mx``
    ``7`` ``val < mx``
    ``8`` ``val <= mx``
    ===== ==================

    """

    return {
        1: lambda: val > mn and mx is None,
        2: lambda: val >= mn and mx is None,
        3: lambda: mn < val < mx,
        4: lambda: mn <= val < mx,
        5: lambda: mn < val <= mx,
        6: lambda: mn <= val <= mx,
        7: lambda: mn is None and val < mx,
        8: lambda: mn is None and val <= mx
    }.get(condition)()
