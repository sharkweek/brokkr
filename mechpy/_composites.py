"""Tools for composite laminates.

Notes
-----
All calculations are performed by using classical laminate theory (CLT) under
the following assumptions:

* Lamina are generally orthotropic.
* Lamina properties and values are assumed to be in the lamina coordinate
  system (i.e. 12-plane).
* Lamina angles (`Lamina.theta`) are in degrees, measured from the laminate
  x-axis to the lamina 1-axis.
* The laminate z-direction is upward, away from the bottom surface.
* Reported strain values are in engineering strain.
* All loads and moments are running values supplied as force or moment per unit
  width.
* Unit systems are consistent (i.e. SI, US, or Imperial).
* Equations and symbol conventions are per NASA-RP-1351 'Basic Mechanics of
  Laminated Composite Plates' and Jones' 'Mechanics Of Composite Materials'.

References
----------
.. [1] NASA-RP-1351, "Basic Mechanics of Laminated Composite Plates"
.. [2] Jones, Robert M. "Mechanics Of Composite Materials"

"""

from numpy import array, zeros
from numpy.linalg import inv
from math import isnan, cos, sin, radians


class Lamina:
    """Lamina(t, E1, E2, nu12, G12, a11, a22, b11, b22)

    An individual lamina.

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
        coefficients of hygral expansion (CTE) in the lamina 1- and 3-
        directions

    """

    __name__ = 'Lamina'
    __slots__ = ['t', 'E1', 'E2', 'nu12', 'G12', 'a11', 'a22', 'b11', 'b22']

    def __init__(self, t, E1, E2, nu12, G12, a11, a22, b11, b22):
        self.t = t
        self.E1 = E1
        self.E2 = E2
        self.nu12 = nu12
        self.G12 = G12
        self.a11 = a11
        self.a22 = a22
        self.b11 = b11
        self.b22 = b22

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
    """Ply(t, E1, E2, nu12, G12, a11, a22, b11, b22)

    A Ply for use in a Laminate.

    Attributes
    ----------
    laminate : Laminate
        the Laminate object the Ply belongs to
    t, E1, E2, nu12, G12, a11, a22, b11
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
        Ply strains due to mechanical, thermal, and hygral loading
    z : float
        Vertical location of the ply midplane in the laminate

    Notes
    -----
    **Attribute Types**
    Not all attributes of are able to be directly modified. Attributes are
    divided into 'base' and 'calculated' values categories and are prescribed
    by the class attributes ``__baseattr__`` and ``__calcattr__``,
    respectively. Base attributes may be set freely, while calculated
    attributes are 'locked' and updated based on the values of base attributes.

    **Assumptions**
    The following assumptions apply to all Ply objects:

    * All Lamina object assumptions apply.
    * Ply z, zk, and zk1 are all measured assuming that positive is upward,
      TOWARD the top surface of the laminate.
    * Theta is in degrees, measured from the laminate x-axis to the lamina 1-
      axis.

      """

    __name__ = 'Ply'
    __baseattr__ = Lamina.__slots__ + ['theta', 'F', 'z']
    __calcattr__ = ['Q', 'Qbar', 'T', 'Tinv', 'e_t', 'e_h', 'e_m', 'laminate']
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
                 b22):
        self.laminate = laminate
        self.z = 0
        self.theta = theta
        self.e_m = zeros((3, 1))
        self.e_t = zeros((3, 1))
        self.e_h = zeros((3, 1))

        super().__init__(t, E1, E2, nu12, G12, a11, a22, b22, b22)
        self.__update()

    def __setattr__(self, attr, val):
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

        # thermal and hygral strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.e_t = array([[self.a11], [self.a22], [0]]) * self.laminate.dT
        self.e_h = array([[self.b11], [self.b22], [0]]) * self.laminate.dM

    @staticmethod
    def new_from_lamina(lamina,
                        laminate,
                        theta):
        """Create a new Ply object from a Lamina object.

        Parameters
        ----------
        lamina : Lamina
            Lamina object to use to create Ply
        laminate : Laminate
            Laminate object the Ply belongs to
        theta : float
            Ply orientation w.r.t. the Laminate coordinate system

        Returns
        -------
        Ply object

        """

        return Ply(laminate=laminate,
                   theta=theta,
                   t=lamina.t,
                   E1=lamina.E1,
                   E2=lamina.E2,
                   nu12=lamina.nu12,
                   G12=lamina.G12,
                   a11=lamina.a11,
                   a22=lamina.a22,
                   b11=lamina.b11,
                   b22=lamina.b22)

    @property
    def zk(self):
        """The vertical location of the lamina's top plane. """
        return self.z + self.t / 2

    @property
    def zk1(self):
        """The vertical location of the lamina's bottom plane."""
        return self.z - self.t / 2

    @zk1.setter
    def zk1(self, new_zk1):
        self.z = new_zk1 + self.t / 2

    # TODO : remove e_tbar and e_hbar from Ply and perform translation in
    #        laminate
    @property
    def e_tbar(self):
        """Transformed thermal strain."""
        return self.Tinv @ self.e_t  # NASA-RP-1351 Eq (90)

    @property
    def e_hbar(self):
        """Transformed hygral strain."""
        return self.Tinv @ self.e_h  # NASA-RP-1351 Eq (95)

    @property
    def e(self):
        """Total strain."""
        return self.e_m + self.e_t + self.e_h

    def failure_index(self, theory='strain'):
        r"""Return the failure index for the given failure theory.

        A ply is considered to fail if the failure index for an applied load
        is equal to or greater than one. The margin of safety for a ply is
        calculated using the failure index and the following equation:

        .. math:: \mathrm{MS} = \frac{\mathrm{1}}{\mathrm{FI}} - 1

        where :math:`\mathrm{MS}` is the margin of safety, and :math:
        `\mathrm{FI}` is the failure index.

        Parameters
        ----------
        theory : {'strain', 'stress', 'Tsai-Hill', 'Tsai-Wu', 'Hoffman'}
            Failure theory for which to calculate a margin of safety.
        index : bool
            Flag to return failure index instead of margin of safety

        Returns
        -------
        float
            The margin of safety (default) or the failure index

        Notes
        -----
        Failure theories are defined by the following relationships:

        Max strain (``theory='strain'``):

        .. math::
           \frac{\varepsilon_1}{F_1} = 1,\quad
           \frac{\varepsilon_2}{F_2} = 1,\quad
           \frac{\gamma_{12}}{F_{12}} = 1

        Max stress (``theory='stress'``):

        .. math::
           \frac{\sigma_1}{F_1} = 1,\quad
           \frac{\sigma_2}{F_2} = 1,\quad
           \frac{\tau_{12}}{F_{12}} = 1

        Tsai-Hill (``theory='Tsai-Hill'``):

        .. math::
           (\frac{\sigma_1}{F_1})^2 + (\frac{\sigma_2}{F_2})^2
           + (\frac{\tau_{12}}{F_{12}})^2
           - \frac{\sigma_1 \sigma_2}{F_1^2} = 1

        Modified Tsai-Wu (``theory='Tsai-Wu'``):

        .. math::
           \frac{\sigma_1^2}{F_1^t F_1^c}
           + \frac{\sigma_2^2}{F_2^t F_2^c}
           + (\frac{1}{F_1^t} - \frac{1}{F_1^c}) \sigma_1
           + (\frac{1}{F_2^t} - \frac{1}{F_2^c}) \sigma_2
           + \frac{\tau_{12}^2}{F_{12}^2} = 1

        Hoffman(``theory='Hoffman'``):

        .. math::
           \frac{\sigma_1^2}{F_1^t F_1^c}
           + \frac{\sigma_2^2}{F_2^t F_2^c}
           - \frac{\sigma_1 \sigma_2}{F_1^t F_1^c}
           + (\frac{F_1^c - F_1^t}{F_1^t F_1^c}) \sigma_1
           + (\frac{F_2^c - F_2^t}{F_2^t F_2^c}) \sigma_2
           + \frac{\tau_{12}^2}{F_{12}^2} = 1

        """
