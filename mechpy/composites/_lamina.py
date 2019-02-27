"""Defines the Lamina class for use with mechpy.

TODO:
-----
* [ ] add error checking for attribute types
* [ ] add unyt units to input (may allow for removal of `ftype`)
* [ ] add a property for the compliance matrix Sk
* [ ] add failure indeces

"""

from numpy import array, zeros
from numpy.linalg import inv
from math import isnan, cos, sin, radians
from mechpy._cogs import ms


class Lamina:
    """
    An individual lamina.

    The ``Lamina`` class exists for material property assignment. To consider
    loading, thermal effects, or deflection, see the ``Ply`` and ``Laminate``
    classes.

    Attributes
    ----------
    t : float
        lamina thickness
    E1, E2, G12 : float
        elastic moduli in the lamina 1-, 2-, and 12-directions
    nu12 : float
        Poisson's ratio in the 12-plane
    a11, a22 : float
        coefficients of thermal expansion (CTE) in the lamina 1- and 2-
        directions
    b11, b22 : float
        coefficients of hygroscopic expansion (CTE) in the lamina 1- and 3-
        directions
    F1, F2, F12 : float
        lamina strengths in each direction
    ftype : {'strain', 'stress'}
        lamina strength type

    """

    __name__ = 'Lamina'
    __slots__ = ['t', 'E1', 'E2', 'nu12', 'G12', 'a11', 'a22', 'b11', 'b22',
                 'F1', 'F2', 'F12', 'ftype']

    def __init__(self, t, E1, E2, nu12, G12, a11, a22, b11, b22, F1=1, F2=1,
                 F12=1, ftype='strain'):
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
        self.ftype = ftype

    @property
    def is_fully_defined(self):
        """Check if lamina is fully defined with proper attr data types."""

        if self.t == 0 or isnan(self.t):
            raise TypeError("lamina.tk must be a non-zero number")
        elif self.E1 == 0 or isnan(self.E1):
            raise TypeError("lamina.E1 must be a non-zero number")
        elif self.E2 == 0 or isnan(self.E2):
            raise TypeError("lamina.E2 must be a non-zero number")
        elif (self.nu12 >= 1) or isnan(self.nu12):
            raise TypeError("""lamina.nu12 must be less than or equal to 1""")
        elif self.G12 == 0 or isnan(self.G12):
            raise TypeError("lamina.G12 must be a non-zero number")
        else:
            return True


class Ply(Lamina):
    """
    A Ply for use in a Laminate.

    Extended ``Lamina``. While the ``Lamina`` class exists for defining
    material properties, the ``Ply`` class is intended to extend its
    functionality further for considering loading and thermal effects. ``Ply``
    instances may exist on their own, but they are intended to function as
    constituent items of a ``Laminate``.

    | **Attribute Types**
    | Not all attributes of are able to be directly modified. Attributes are
      divided into 'base' and 'calculated' values categories and are
      prescribed by the class attributes ``__baseattr__`` and ``__calcattr__``,
      respectively. Base attributes may be set freely, while calculated
      attributes are 'locked' and updated based on the values of base
      attributes.

    | **Assumptions**
    | The following assumptions apply to all Ply objects:

      * Ply z, zk, and zk1 are all measured assuming that positive is upward,
        TOWARD the top surface of the laminate.
      * Theta is in degrees, measured from the laminate x-axis to the lamina
        1- axis.

    Attributes
    ----------
    laminate : Laminate
        the Laminate object the Ply belongs to
    t, E1, E2, nu12, G12, a11, a22, b11, F1, F2, F12, ftype
        See ``Lamina`` attribute definitions
    theta : float
        the angle the Ply is oriented in w.r.t. the Laminate coordinate system
    Q : 3x1 numpy.ndarray
        Ply stiffness matrix in the Ply coordinate system
    Qbar : 3x1 numpy.ndarray
        Ply stiffness matrix in the Laminate coordinate system
    T : 3x1 numpy.ndarray
        Ply transformation matrix
    Tinv : 3x1 numpy.ndarray
        Inverse of the Ply transformation matrix
    e_m, e_t, e_h : 3x1 numpy.ndarray
        Ply strains due to mechanical, thermal, and hygroscopic loading
    s_m, s_t, s_h : 3x1 numpy.ndarray
        Ply stresses due to mechanical, thermal, and hygroscopic loading
    z : float
        Vertical location of the ply midplane in the laminate
    failure_theory : {'strain', 'stress', 'Tsai-Hill'}
        The failure theory for calculating the failure index
    failure_index : float
        the failure index

    """

    __name__ = 'Ply'
    __baseattr__ = Lamina.__slots__ + ['theta', 'z', 'e_m', 's_m']
    __calcattr__ = ['Q', 'Qbar', 'T', 'Tinv', 'e_t', 'e_h', 's_t', 's_h',
                    'laminate', 'failure_theory', 'failure_index']
    __slots__ = __baseattr__ + __calcattr__ + ['__locked']

    def __unlock(func):
        """Decorate methods to unlock attributes.

        Parameters
        ----------
        func : bool
            The function to unlock.

        Returns
        -------
        function
            An unlocked function.

        """

        def wrapper(self, *args, **kwargs):
            super().__setattr__('_Ply__locked', False)
            func(self, *args, **kwargs)
            super().__setattr__('_Ply__locked', True)
        return wrapper

    @__unlock
    def __init__(self, laminate, t, theta, E1, E2, nu12, G12, a11, a22, b11,
                 b22, F1=1, F2=1, F12=1, ftype='strain',
                 failure_theory='strain'):
        """Extend ``__init__`` to account for Ply-only attributes."""

        self.laminate = laminate
        self.z = 0
        self.theta = theta
        self.e_m = zeros((3, 1))
        self.e_t = zeros((3, 1))
        self.e_h = zeros((3, 1))
        self.s_m = zeros((3, 1))
        self.s_t = zeros((3, 1))
        self.s_h = zeros((3, 1))
        self.failure_theory = failure_theory
        self.failure_index = 0

        super().__init__(t, E1, E2, nu12, G12, a11, a22, b22, b22, F1, F2, F12,
                         ftype)
        self.__update()

    def __setattr__(self, attr, val):
        """Extend ``__setattr__`` to protect calculated attributes."""

        if self.__locked:
            # udpate laminate after updated properties are set
            if attr in self.__baseattr__:
                super().__setattr__(attr, val)
                self.__update()

            # don't set protected values
            elif attr in self.__calcattr__:
                raise AttributeError(self.__name__ + ".%s" % attr
                                     + " is a derived value and cannot be set")

            # update the laminate
            if self.laminate:
                self.laminate._Laminate__update()

        else:
            super().__setattr__(attr, val)

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

        self.Q = array([[q11, q12, 0], [q12, q22, 0], [0, 0, q66]])

        # the transformation matrix and its inverse
        # create intermediate trig terms
        m = cos(radians(self.theta))
        n = sin(radians(self.theta))

        # create transformation matrix and inverse
        self.T = array([[m**2, n**2, 2 * m * n],
                        [n**2, m**2, -2 * m * n],
                        [-m * n, m * n, m**2 - n**2]])
        self.Tinv = inv(self.T)

        # the transformed reduced stiffness matrix (laminate coordinate system)
        # Jones, Eq (2.84)
        self.Qbar = self.Tinv @ self.Q @ self.Tinv.T

        # thermal and hygroscopic strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.e_t = array([[self.a11], [self.a22], [0]]) * self.laminate.dT
        self.e_h = array([[self.b11], [self.b22], [0]]) * self.laminate.dM

        # thermal and hygroscopic stresses in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.s_t = self.Q @ self.e_t
        self.s_h = self.Q @ self.e_h

        # calculate failure index
        self.failure_index = self.calc_failure_index(self.failure_theory,
                                                     self.F1,
                                                     self.F2,
                                                     self.F12)

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
                   ftype=lamina.ftype)

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
    def e(self):
        """Total strain."""
        return self.e_m + self.e_t + self.e_h

    @property
    def s(self):
        """Total stress."""
        return self.s_m + self.s_t + self.s_h

    @staticmethod
    def calc_failure_index(theory, F1, F2, F12):
        r"""Calculate the failure index for a given failure theory.

        Parameters
        ----------
        theory : {'strain', 'stress', 'Tsai-Hill'}
            failure theory for which to calculate a failure index
        F1, F2, F12 : float
            strengths of the material

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

    @property
    def margin(self):
        """The margin of safety."""
        return ms(1, self.failure_index)
