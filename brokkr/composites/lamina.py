"""Defines the Lamina class for use with brokkr.

TODO:
-----
* [ ] add a property for the compliance matrix Sk in Lamina
* [ ] move stiffness matrix Q to Lamina
* [ ] add failure indices

"""

from numpy import zeros, cos, sin, ndarray
from numpy.linalg import inv
from brokkr.mech_math import ms, out_of_bounds
from brokkr.config import USYS
from brokkr.core.exceptions import (
    UnitDimensionError,
    BoundedValueError,
    CalculatedAttributeError
)
from unyt import unyt_array, UnitSystem
from unyt.dimensions import (
    length,
    pressure,
    dimensionless,
    temperature,
    angle
)


__all__ = ['Lamina', 'Ply']


class Lamina:
    """An individual lamina.

    The ``Lamina`` class exists for material property assignment. To consider
    loading, thermal effects, or deflection, see the ``Ply`` and ``Laminate``
    classes.

    Parameters
    ----------
    t : float
        lamina thickness (``dim: (length)``)
    E1, E2, G12 : float
        elastic moduli in the lamina 1-, 2-, and 12-directions (``dim:
         (mass)/((length)*(time)**2)``)
    nu12 : float
        Poisson's ratio in the 12-plane (``dim: (dimensionless)``)
    a11, a22 : float
        coefficients of thermal expansion (CTE) in the lamina 1- and 2-
        directions (``dim: (dimensionless)``)
    b11, b22 : float
        coefficients of hygroscopic expansion (CTE) in the lamina 1- and 3-
        directions (``dim: (dimensionless)``)
    F1, F2, F12 : float
        lamina strengths in each direction; may be in strain or stress (``dim:
         (dimensionless) or (mass)/((length)*(time)**2)``)
    usys : unyt.UnitSystem
        the unit system for all dimensioned parameters

    Raises
    ------
    BoundedValueError
        if assigned attribute value is not within required boundaries
    UnitDimensionError
        if assigned attribute units have incorrect dimensionality

    Notes
    -----
    Parameters may be entered including or not including units. If units are
    not included, they will be assigned from the default unit system
    (``usys``). Dimensionless units such as strain and percent moisture gain
    are assigned a ``dimensionless`` unit. Calculated attributes produce the
    appropriate units based on the base attributes.

    Parameters that are entered with units, are checked for the correct
    dimensionality and converted to the default unit system before being stored
    in the instance for consistency.

    """

    _dims = {
        't': (length,),
        'E1': (pressure,),
        'E2': (pressure,),
        'nu12': (dimensionless,),
        'G12': (pressure,),
        'a11': (1 / temperature,),
        'a22': (1 / temperature,),
        'b11': (dimensionless,),  # percent moisture change
        'b22': (dimensionless,),
        'F1': (dimensionless, pressure),  # dimensionless for strain
        'F2': (dimensionless, pressure),
        'F12': (dimensionless, pressure)
    }
    # define attribute limits for use with `out_of_bounds()`
    _limits = {
        't': {'mn': 0, 'mx': None, 'condition': 'g'},
        'E1': {'mn': 0, 'mx': None, 'condition': 'g'},
        'E2': {'mn': 0, 'mx': None, 'condition': 'g'},
        'nu12': {'mn': 0, 'mx': 1, 'condition': 'g-l'},
        'G12': {'mn': 0, 'mx': None, 'condition': 'g'}
    }
    __slots__ = ('usys', *_dims)

    def __init__(self, t, E1, E2, nu12, G12, a11, a22, b11, b22, F1=0, F2=0,
                 F12=0, usys=USYS):

        super().__setattr__('usys', usys)
        self.t = t
        self.E1 = E1
        self.E2 = E2
        self.nu12 = nu12
        self.G12 = G12
        self.a11 = a11
        self.a22 = a22
        self.b11 = b11
        self.b22 = b22
        self.F1 = F1
        self.F2 = F2
        self.F12 = F12

    def __setattr__(self, name, value):
        """Extend __setattr__() to validate units."""

        # set dimensions for required attributes
        if name in self._dims:
            correct_dim = self._dims.get(name)
            # check if object has units
            if not hasattr(value, 'units'):
                value *= self.usys[correct_dim[0]]  # assign first in tuple

            # if units assigned, check units for dimensionality
            else:
                if value.units.dimensions not in correct_dim:
                    raise UnitDimensionError(
                        name, ' OR '.join([str(i) for i in correct_dim])
                        )
                else:
                    # units are converted to master unit system for
                    # consistency when creating unit_arrays in `Ply` class
                    value.convert_to_base(self.usys)

            # make sure value is within limits
            if name in self._limits:
                if out_of_bounds(value.value, **self._limits.get(name)):
                    raise BoundedValueError(
                        name, **self._limits.get(name)
                        )

        # catch unit system change and convert all attributes with units
        if name == 'usys':
            # validate unit system
            if type(value) != UnitSystem:
                raise AttributeError(f"`{name}` must be a `UnitSystem`")
            elif value['temperature'].base_offset != 0:
                raise AttributeError(
                    "Unit system must have an absolute temperature unit "
                    + "(e.g. R or K)"
                )

            for each in self._dims:
                getattr(self, each).convert_to_base(value)

        super().__setattr__(name, value)

    def __repr__(self):
        r = f"{type(self).__name__}:"
        for each in self.__slots__:
            r += f"\n    {each + ':':<6}{getattr(self, each)}"
        return r


class Ply(Lamina):
    """A Ply for use in a Laminate.

    Extended ``Lamina``. While the ``Lamina`` class exists for defining
    material properties, the ``Ply`` class is intended to extend its
    functionality further for considering loading and thermal effects. ``Ply``
    instances may exist on their own, but they are intended to function as
    constituent items of a ``Laminate``.

    | **Attribute Types**
    | Not all attributes of are able to be directly modified. Attributes are
      divided into 'base' and 'calculated' values categories and are
      prescribed by the class attributes ``_base_attr`` and ``_calc_attr``,
      respectively. Base attributes may be set freely, while calculated
      attributes are 'locked' and updated based on the values of base
      attributes.

    | **Assumptions**
    | The following assumptions apply to all Ply objects:

      * Ply z, zk, and zk1 are all measured assuming that positive is upward,
        TOWARD the top surface of the laminate.
      * Theta is in degrees, measured from the laminate x-axis to the lamina
        1- axis.

    Parameters
    ----------
    laminate : Laminate
        the Laminate object the Ply belongs to
    theta : float
        the angle the Ply is oriented in w.r.t. the Laminate coordinate system
        (``dim: (angle)``)
    failure_theory : {'strain', 'stress', 'Tsai-Hill'}
        The failure theory for calculating the failure index

    Attributes
    ----------
    Q : 3x1 unyt.unyt_array
        Ply stiffness matrix in the Ply coordinate system (``dim: (mass)
        /((length)*(time)**2)``)
    Qbar : 3x1 unyt.unyt_array
        Ply stiffness matrix in the Laminate coordinate system (``dim:
         (mass)/((length)*(time)**2)``)
    T : 3x1 unyt.unyt_array
        Ply transformation matrix (``dim: (dimensionless)``)
    Tinv : 3x1 unyt.unyt_array
        Inverse of the Ply transformation matrix (``dim: (dimensionless)``)
    e_m, e_t, e_h : 3x1 unyt.unyt_array
        Ply strains due to mechanical, thermal, and hygroscopic loading (``dim:
         (dimensionless)``)
    z : float
        Vertical location of the ply midplane in the laminate (``dim:
        (length)``)
    failure_index : float
        the failure index (``dim: dimensionless``)

    """

    _dims = {
        **Lamina._dims,
        'theta': (angle,),
        'z': (length,),
        'e_m': (dimensionless,),
        's': (pressure,),
        'Q': (pressure,),
        'Qbar': (pressure,),
        'T': (dimensionless,),
        'e_t': (dimensionless,),
        'e_h': (dimensionless,),
    }

    _base_attr = ('theta', 'z', 'usys', *Lamina._dims)
    _calc_attr = ('Q', 'Qbar', 'T', 'Tinv', 'e_t', 'e_h', 's',
                  'laminate', 'failure_theory', 'failure_index')
    __slots__ = _base_attr + _calc_attr + ('e_m', '__locked',)

    def __unlock(locked_func):
        """Decorate methods to unlock attributes.

        Parameters
        ----------
        locked_func : bool
            The function to unlock.

        Returns
        -------
        function
            An unlocked function.

        """

        def unlocked_func(self, *args, **kwargs):
            super().__setattr__('_Ply__locked', False)
            locked_func(self, *args, **kwargs)
            super().__setattr__('_Ply__locked', True)
        return unlocked_func

    @__unlock
    def __init__(self, laminate, t, theta, E1, E2, nu12, G12, a11, a22, b11,
                 b22, F1=0, F2=0, F12=0, usys=USYS, failure_theory='strain'):
        """Extend ``__init__`` to account for Ply-only attributes."""

        super().__init__(t, E1, E2, nu12, G12, a11, a22, b22, b22, F1, F2, F12,
                         usys)
        self.laminate = laminate
        self.z = 0
        self.theta = theta
        self.e_m = zeros((3, 1))
        self.e_t = zeros((3, 1))
        self.e_h = zeros((3, 1))
        self.s = zeros((3, 1))
        self.failure_theory = failure_theory
        self.failure_index = 0

        self.__update()

    def __setattr__(self, name, attr):
        """Extend ``__setattr__`` to protect calculated attributes."""

        if self.__locked:
            # udpate ply and laminate after updated properties are set
            if name in self._base_attr:
                super().__setattr__(name, attr)
                self.__update()
                if self.laminate:
                    self.laminate._Laminate__update()

            # don't set protected values
            elif name in self._calc_attr:
                raise CalculatedAttributeError(name)

        else:
            super().__setattr__(name, attr)

    def __repr__(self):
        r = f'{type(self).__name__}'
        for each in sorted(self._dims):
            attr = getattr(self, each)

            if issubclass(attr.__class__, ndarray):
                lines = attr.__str__().splitlines()
                r += f"\n    {each + ':':<7}{lines[0]}"

                for i in range(1, len(lines)):
                    r += "\n{0:<11}{1:}".format(' ', lines[i])

            else:
                r += f"\n    {each + ':':<7}{attr}"

        r += (
            f"\n    {'failure theory' + ': '}{getattr(self, 'failure_theory')}"
            )
        r += f"\n    {'failure index' + ': '}{getattr(self, 'failure_index')}"

        return r

    @__unlock
    def __update(self):
        """Update calculated attributes."""

        # on-axis reduced stiffness matrix, Q
        # NASA-RP-1351, Eq (15)
        nu21 = self.nu12 * self.E2 / self.E1  # Jones, Eq (2.67)
        q11 = self.E1 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E2 / (1 - self.nu12 * nu21)
        q22 = self.E2 / (1 - self.nu12 * nu21)
        q66 = self.G12

        self.Q = unyt_array([[q11, q12, 0], [q12, q22, 0], [0, 0, q66]],
                            self.usys['pressure'])

        # the transformation matrix and its inverse
        # create intermediate trig terms
        m = cos(self.theta)
        n = sin(self.theta)

        # create transformation matrix and inverse
        self.T = unyt_array([[m**2, n**2, 2 * m * n],
                             [n**2, m**2, -2 * m * n],
                             [-m * n, m * n, m**2 - n**2]])
        self.Tinv = inv(self.T)

        # the transformed reduced stiffness matrix (laminate coordinate system)
        # Jones, Eq (2.84)
        self.Qbar = (self.Tinv @ self.Q @ self.Tinv.T).v * self.Q.units

        # thermal and hygroscopic strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.e_t = (
            unyt_array([[self.a11], [self.a22], [0]],
                       1 / self.usys['temperature'])
            ) * self.laminate.dT
        self.e_h = (
            unyt_array([[self.b11], [self.b22], [0]], dimensionless)
            ) * self.laminate.dM

        # calculate failure index
        #self.failure_index = self.calc_failure_index()

    @classmethod
    def from_lamina(cls, lamina, laminate, theta):
        """Create a new Ply object from a Lamina object.

        Parameters
        ----------
        lamina : Lamina
            ``Lamina`` from which to create ``Ply``
        laminate : Laminate
            ``Laminate`` object the ``Ply`` belongs to
        theta : float
            ``Ply`` orientation w.r.t. the ``Laminate`` coordinate system

        Returns
        -------
        ``Ply`` object

        """

        # ensure that unit system is consistent
        lamina.usys = laminate.usys

        return cls(laminate=laminate,
                   theta=theta,
                   t=lamina.t,
                   E1=lamina.E1,
                   E2=lamina.E2,
                   nu12=lamina.nu12,
                   G12=lamina.G12,
                   a11=lamina.a11,
                   a22=lamina.a22,
                   b11=lamina.b11,
                   b22=lamina.b22,
                   F1=lamina.F1,
                   F2=lamina.F2,
                   F12=lamina.F12,
                   usys=lamina.usys)

    @property
    def zk(self):
        """The vertical location of the lamina's top plane."""
        return self.z + self.t / 2

    @property
    def zk1(self):
        """The vertical location of the lamina's bottom plane."""
        return self.z - self.t / 2

    @zk1.setter
    def zk1(self, new_zk1):
        self.z = new_zk1 + self.t / 2

    @property
    def S(self):
        """The compliance matrix."""
        return inv(self.Q.v) * (1 / self.Q.units)

    @staticmethod
    def calc_failure_index(theory, s1, s2, s3, F1, F2, F12):
        r"""Calculate the failure index for a given failure theory.

        Parameters
        ----------
        theory : {'strain', 'stress', 'Tsai-Hill'}
            failure theory for which to calculate a failure index
        s1, s2, s3 : float
            applied strain or stress values (``dim: (dimensionless) OR ``dim:
            (mass)/((length)*(time)**2)``)
        F1, F2, F12 : float
            strengths of the material (``dim: (dimensionless) OR ``dim:
            (mass)/((length)*(time)**2)``)

        .. note:: Applied and strength values must have the same
           dimensionality.

        Returns
        -------
        float
            The failure index

        Notes
        -----
        A ply is considered to fail if the failure index for an applied load
        is equal to or greater than one. Failure indicies for each failure
        theory are calculated according to the following equations per the
        BJSFM User Manual [#3]_.

        Max strain (``theory='strain'``):

        .. math::
           \mathrm{FI} = \mathrm{min}
           \left( \frac{\varepsilon_1}{F_1},\quad
           \frac{\varepsilon_2}{F_2},\quad
           \frac{\gamma_{12}}{F_{12}} \right)

        Max stress (``theory='stress'``):

        .. math::
           \mathrm{FI} = \mathrm{min}
           \left( \frac{\sigma_1}{F_1} = 1,\quad
           \frac{\sigma_2}{F_2} = 1,\quad
           \frac{\tau_{12}}{F_{12}} \right)

        Tsai-Hill (``theory='Tsai-Hill'``):

        .. math::
           \mathrm{FI} = \left( \frac{\sigma_1}{F_1} \right)^{2}
           + \left( \frac{\sigma_2}{F_2} \right)^{2}
           + \left( \frac{\tau_{12}}{F_{12}} \right)^{2}
           - \frac{\sigma_1 \sigma_2}{F_1^2}

        .. note:: Modified Tsai-Wu and Hoffman criteria are not supported
           as they require separate strength values for tension and compression

        References
        ----------
        .. [#3] Ogonowski, J.M, *Effect of Variances and Manufacturing
           Tolerances on the Design Strength and Life of Mechanically Fastened
           Composite Joints, Volume 3 - Bolted Joint Stress Field Model
           (BJSFM) Computer Program User's Manual*, McDonnell Aircraft
           Company, AFWAL-TR-81-3041 VOLUME 3, pp. 8-9, 15 April 1981

        """

        if theory == 'strain':
            fi = 1

        elif theory == 'stress':
            fi = 1

        elif theory == 'Tsai-Hill':
            fi = 1

        else:
            raise ValueError(
                "`theory` must be 'strain', 'stress', or 'Tsai-Hill'."
                )

        return fi

    @property
    def margin(self):
        """The margin of safety."""
        return ms(1, self.failure_index)
