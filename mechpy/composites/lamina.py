"""Defines the Lamina class for use with mechpy.

All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions of the lamina, respectively. It is also assumed
that lamina are generally orthotropic.

Equations and symbol conventions are per 'Mechanics Of Composite Materials' by
Robert M. Jones and 'NASA-RP-1351 BASIC MECHANICS OF LAMINATED COMPOSITE
PLATES' as noted.

TODO:
-----
* [ ] add error checking for attribute types
* [ ] add a property for the compliance matrix Sk
* [ ] add stress attributes for failure criteria
"""

import numpy as np
import math


class Lamina:
    """An individual lamina or "ply".

    Each lamina is assumed to be orthotropic.
    """

    # class constants
    # the simple matrix (Jones, Eq (2.78) and (2.79).
    r = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 2]])

    r_inv = np.linalg.inv(r)

    def __init__(self,
                 t=0.010,     # lamina thickness
                 theta=0,     # orientation of lamina in degrees
                 E1=10e6,     # Young's modulus in 1-direction
                 E2=10e6,     # Young's modulus in 2-direction
                 nu12=0.1,    # Poisson's ratio in 12-plane
                 G12=725e3,   # shear modulus in 12-plane
                 a11=2.1e-6,  # coeff. of thermal expansion in 1-direction
                 a22=2.1e-6,  # coeff. of thermal expansion in 2-direction
                 b11=3e-8,    # coeff. of moisture expansion in 1-direction
                 b22=3e-8,    # coeff. of moisture expansion in 2-direction
                 dT=0,        # change in temperature
                 dm=0):       # moisture absorption
        """Initialize the Lamina instance."""

        self.ID = 0
        self.__t = t
        self.__z = 0
        self.__theta = theta
        self.__E1 = E1
        self.__E2 = E2
        self.__nu12 = nu12
        self.__G12 = G12
        self.__alpha = np.array([[a11],
                                 [a22],
                                 [0]],
                                dtype=float)
        self.__beta = np.array([[b11],
                                [b22],
                                [0]],
                               dtype=float)
        self.__dT = dT
        self.__dm = dm
        self.__e_m = np.zeros((3, 1))
        self.__e_t = np.zeros((3, 1))
        self.__e_h = np.zeros((3, 1))
        self.__e_mbar = np.zeros((3, 1))
        self.__e_tbar = np.zeros((3, 1))
        self.__e_hbar = np.zeros((3, 1))

        self.__update()

    @property
    def theta(self):
        """The lamina orientation in degrees."""

        return self.__theta

    @theta.setter
    def theta(self, new_theta):
        """The lamina orientation in degrees."""

        self.__theta = new_theta
        self.__update()

    @property
    def E1(self):
        """Young's modulus in the lamina 1-direction."""

        return self.__E1

    @E1.setter
    def E1(self, new_E1):
        """Young's modulus in the 1-direction."""

        self.__E1 = new_E1
        self.__update()

    @property
    def E2(self):
        """Young's modulus in the lamina 2-direction."""

        return self.__E2

    @E2.setter
    def E2(self, new_E2):
        """Young's modulus in the lamina 2-direction."""

        self.__E2 = new_E2
        self.__update()

    @property
    def nu12(self):
        """Poisson's ratio in the lamina 12-plane."""

        return self.__nu12

    @nu12.setter
    def nu12(self, new_nu12):
        """Poisson's ratio in the lamina 12-plane."""

        self.__nu12 = new_nu12
        self.__update()

    @property
    def G12(self):
        """Shear modulus in the lamina 12-plane."""

        return self.__G12

    @G12.setter
    def G12(self, new_G12):
        """Shear modulus in the lamina 12-plane."""

        self.__G12 = new_G12
        self.__update()

    @property
    def alpha(self):
        """The thermal expansion matrix for the lamina.

        Lamina are assumed to have orthotropic expansion coefficients.

        Example: [[a11],
                  [a22],
                  [0]]
        """

        return self.__alpha

    @alpha.setter
    def alpha(self, new_alpha_k):
        """The thermal expansion matrix for the lamina.

        Lamina are assumed to have orthotropic expansion coefficients.

        Example: [[a11],
                  [a22],
                  [0]]
        """

        self.__alpha = new_alpha_k
        self.__update()

    @property
    def beta(self):
        """The hygroscopic expansion matrix for the lamina.

        Lamina are assumed to have orthotropic expansion coefficients.

        Example: [[b11],
                  [b22],
                  [0]]
        """

        return self.__beta

    @beta.setter
    def beta(self, new_beta_k):
        """The hygroscopic expansion matrix for the lamina.

        Lamina are assumed to have orthotropic expansion coefficients.

        Example: [[b11],
                  [b22],
                  [0]]
        """

        self.__beta = new_beta_k
        self.__update()

    @property
    def dT(self):
        """The relative change in temperature of the lamina."""

        return self.__dT

    @dT.setter
    def dT(self, new_dT):
        """The relative change in temperature of the lamina."""

        return self.__dT

    @property
    def dm(self):
        """The relative change in moisture of the lamina."""

        return self.__dm

    @dm.setter
    def dm(self, new_dm):
        """The relative change in moisture of the lamina."""

        self.__dm = new_dm
        self.__update()

    @property
    def t(self):
        """The thickness of the lamina."""

        return self.__t

    @t.setter
    def t(self, new_t):
        """The thickness of the lamina."""

        self.__t = new_t
        self.__update()

    @property
    def z(self):
        """The vertical location of the lamina's midplane.

        z is measured from the bottom surface of the laminate."""

        return self.__z

    @z.setter
    def z(self, new_z):
        """The vertical location of the lamina's midplane.

        z is measured from the bottom surface of the laminate."""

        self.__z = new_z
        self.__update()

    @property
    def zk(self):
        """The vertical location of the lamina's top plane.

        zk is measured from the bottom surface of the laminate."""

        return self.__zk

    @zk.setter
    def zk(self, new_zk):
        """The vertical location of the lamina's top plane.

        zk is measured from the bottom surface of the laminate."""

        self.__z = new_zk - self.__t / 2
        self.__update()

    @property
    def zk1(self):
        """The vertical location of the lamina's bottom plane.

        zk1 is measured from the bottom surface of the laminate."""

        return self.__zk1

    @zk1.setter
    def zk1(self, new_zk1):
        """The vertical location of the lamina's bottom plane.

        zk1 is measured from the bottom surface of the laminate."""

        self.__z = new_zk1 + self.__t / 2
        self.__update()

    @property
    def Qk(self):
        """The on-axis reduced stiffness matrix of the lamina."""

        return self.__Qk

    @property
    def Qk_bar(self):
        """The transformed reduced stiffness matrix of the lamina."""

        return self.__Qk_bar

    @property
    def T(self):
        """The lamina transformation matrix."""

        return self.__T

    @property
    def Tinv(self):
        """The inverse of the transformation matrix."""

        return self.__Tinv

    @property
    def e_m(self):
        """Mechanical strain in the lamina orientation."""

        return self.__e_m

    @property
    def e_t(self):
        """Thermal strain in the lamina orientation."""

        return self.__e_t

    @property
    def e_h(self):
        """Hygroscopic strain in the lamina orientation."""

        return self.__e_h

    @property
    def e_mbar(self):
        """Mechanical strain in the laminate orientation."""

        return self.__e_mbar

    @property
    def e_tbar(self):
        """Thermal strain in the laminate orientation."""

        return self.__e_tbar

    @property
    def e_hbar(self):
        """Hygroscopic strain in the laminate orientation."""

        return self.__e_hbar

    @e_m.setter
    def e_mbar(self, new_strain):
        """Mechanical strain in the lamina orientation."""

        self.__e_m = np.matmul(self.__Tinv, new_strain)
        self.__update()

    @property
    def isFullyDefined(self):
        """Check if lamina is fully defined with proper attr data types."""

        if self.__t == 0 or math.isnan(self.__t):
            raise TypeError("lamina.tk must be a non-zero number")
        elif math.isnan(self.__theta):
            raise TypeError("lamina.theta must be a number (in degrees)")
        elif self.__E1 == 0 or math.isnan(self.__E1):
            raise TypeError("lamina.E1 must be a non-zero number")
        elif self.E2 == 0 or math.isnan(self.E2):
            raise TypeError("lamina.E2 must be a non-zero number")
        elif (self.nu12 <= 0 or self.nu12 > 1) or math.isnan(self.nu12):
            raise TypeError("""lamina.nu12 must be a non-zero number between
                            zero and 1""")
        elif self.G12 == 0 or math.isnan(self.G12):
            raise TypeError("lamina.G12 must be a non-zero number")
        else:
            return True

    def __update(self):
        """Update calculated properties when new values are assigned."""

        # upper and lower boundaries of lamina based on z
        # Jones, Figure 4-8
        self.__zk = self.z + self.__t / 2
        self.__zk1 = self.z - self.__t / 2

        # on-axis reduced stiffness matrix, Qk
        # Jones, Eq (2.70) and (2.71)
        nu21 = self.nu12 * self.E2 / self.__E1  # Jones, Eq (2.67)
        q11 = self.__E1 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E2 / (1 - self.nu12 * nu21)
        q22 = self.E2 / (1 - self.nu12 * nu21)
        q66 = self.G12

        self.__Qk = np.array([[q11, q12, 0],
                              [q12, q22, 0],
                              [0, 0, q66]])

        # the transformation matrix and its inverse
        # create intermediate trig terms
        m = math.cos(math.radians(self.__theta))
        n = math.sin(math.radians(self.__theta))

        # create transformation matrix and inverse
        self.__T = np.array([[m**2, n**2, 2*m*n],
                             [n**2, m**2, -2*m*n],
                             [-m*n, m*n, m**2 - n**2]])
        self.__Tinv = np.linalg.inv(self.__T)

        # the transformed reduced stiffness matrix (laminate coordinate system)
        # Jones, Eq (2.83)
        T_m1 = np.matmul(np.matmul(self.r, self.__T), self.r_inv)
        self.__Qk_bar = np.matmul(np.matmul(self.__Tinv, self.__Qk), T_m1)

        # thermal and hygral strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.__e_t = self.__alpha * self.__dT
        self.__e_h = self.__beta * self.__dm

        # the transformed CTE and hygral matrices (laminate coordinate system)
        # NASA-RP-1351 Eq (90) and (95)
        self.__e_tbar = np.matmul(self.__T, self.__alpha)
        self.__e_hbar = np.matmul(self.__T, self.__beta)


def t_matrix(theta):
    """Return the transformation matrix for a given theta theta.

    The transformation matrix is used to transform stresses and tensor strains
    from the lamina material orientation to the laminate coordinate system. The
    inverse of the matrix may be used to transform back from the laminate
    coordinate system to that of the lamina. Theta is measured counterclockwise
    from the laminate x-axis and assumed to be in degrees.

    See 'Mechanics Of Composite Materials' by Robert M. Jones, Eq (2.76) for
    the definition of the transformation matrix.
    """

    # create intermediate trig terms
    m = math.cos(math.radians(theta))
    n = math.sin(math.radians(theta))

    # create transformation matrix and inverse
    return np.array([[m**2, n**2, 2*m*n],
                     [n**2, m**2, -2*m*n],
                     [-m*n, m*n, m**2 - n**2]])
