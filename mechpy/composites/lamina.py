"""Defines the Lamina class for use with mechpy.

Assumptions
-----------
* Lamina are generally orthotropic.
* Lamina properties and values are assumed to be in the lamina coordinate system
  (i.e. 12-plane).
* Theta is in degrees, measured from the laminate x-axis to the lamina 1-axis.
* The laminate z-direction is upward, away from the bottom surface.
* Reported strain values are in engineering strain.
* All loads and moments are running values supplied as force or moment per unit
  width.
* Unit systems are consistent (i.e. SI, US, or Imperial).
* Equations and symbol conventions are per NASA-RP-1351 'Basic Mechanics of
  Laminated Composite Plates' and Jones' 'Mechanics Of Composite Materials'.

TODO:
-----
* [ ] add functionality to register laminate property to forwad updates to
      Laminate classes
* [ ] add error checking for attribute types
* [ ] add a property for the compliance matrix Sk
* [ ] add failure indeces
* [ ] remove laminate orientation calculations and have them be returned as
      calculated properties
"""

import numpy as np
from numpy import matmul as mm
from math import isnan, sin, cos, radians


class Lamina:
    """An individual lamina or "ply".

    Each lamina is assumed to be orthotropic.
    """

    def __init__(self,
                 t=1,     # lamina thickness
                 theta=0, # orientation of lamina in degrees
                 E1=1,    # Young's modulus in 1-direction
                 E2=1,    # Young's modulus in 2-direction
                 nu12=1,  # Poisson's ratio in 12-plane
                 G12=1,   # shear modulus in 12-plane
                 a11=1,   # coeff. of thermal expansion in 1-direction
                 a22=1,   # coeff. of thermal expansion in 2-direction
                 b11=1,   # coeff. of moisture expansion in 1-direction
                 b22=1,   # coeff. of moisture expansion in 2-direction
                 dT=0,    # change in temperature
                 dM=0):   # moisture absorption
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
        self.__dM = dM
        self.__e_m = np.zeros((3, 1))
        self.__e_m = np.zeros((3, 1))
        self.__e_t = np.zeros((3, 1))
        self.__e_h = np.zeros((3, 1))
        self.__e_mbar = np.zeros((3, 1))
        self.__e_tbar = np.zeros((3, 1))
        self.__e_hbar = np.zeros((3, 1))

        self.__update()

    def __update(self):
        """Update calculated properties when new values are assigned."""

        # upper and lower boundaries of lamina based on z
        # NASA-RP-1351, Figure 9
        self.__zk = self.__z + self.__t / 2
        self.__zk1 = self.__z - self.__t / 2

        # on-axis reduced stiffness matrix, Q
        # NASA-RP-1351, Eq (15)
        nu21 = self.nu12 * self.E2 / self.__E1  # Jones, Eq (2.67)
        q11 = self.__E1 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E2 / (1 - self.nu12 * nu21)
        q22 = self.E2 / (1 - self.nu12 * nu21)
        q66 = self.G12

        self.__Q = np.array([[q11, q12, 0],
                             [q12, q22, 0],
                             [0, 0, q66]])

        # the transformation matrix and its inverse
        # create intermediate trig terms
        m = cos(radians(self.__theta))
        n = sin(radians(self.__theta))

        # create transformation matrix and inverse
        self.__T = np.array([[m**2, n**2, 2 * m * n],
                             [n**2, m**2, -2 * m * n],
                             [-m * n, m * n, m**2 - n**2]])
        T_inv = np.linalg.inv(self.__T)

        # the transformed reduced stiffness matrix (laminate coordinate system)
        # Jones, Eq (2.84)
        self.__Q_bar = mm(mm(T_inv, self.__Q), T_inv.transpose())

        # thermal and hygral strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.__e_t = self.__alpha * self.__dT
        self.__e_h = self.__beta * self.__dM

        # the transformed CTE and hygral matrices (laminate coordinate system)
        # NASA-RP-1351 Eq (90) and (95)
        self.__e_tbar = mm(T_inv, self.__alpha)
        self.__e_hbar = mm(T_inv, self.__beta)

        # calculate total strain
        self.__e = self.__e_m + self.__e_t + self.__e_h

    @property
    def is_fully_defined(self):
        """Check if lamina is fully defined with proper attr data types."""

        if self.__t == 0 or isnan(self.__t):
            raise TypeError("lamina.tk must be a non-zero number")
        elif isnan(self.__theta):
            raise TypeError("lamina.theta must be a number (in degrees)")
        elif self.__E1 == 0 or isnan(self.__E1):
            raise TypeError("lamina.E1 must be a non-zero number")
        elif self.E2 == 0 or isnan(self.E2):
            raise TypeError("lamina.E2 must be a non-zero number")
        elif (self.nu12 <= 0 or self.nu12 > 1) or isnan(self.nu12):
            raise TypeError("""lamina.nu12 must be a non-zero number between
                            zero and 1""")
        elif self.G12 == 0 or isnan(self.G12):
            raise TypeError("lamina.G12 must be a non-zero number")
        else:
            return True

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
    def alpha(self, new_alpha):
        """The thermal expansion matrix for the lamina.

        Lamina are assumed to have orthotropic expansion coefficients.

        Example: [[a11],
                  [a22],
                  [0]]
        """

        self.__alpha = new_alpha
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
    def beta(self, new_beta):
        """The hygroscopic expansion matrix for the lamina.

        Lamina are assumed to have orthotropic expansion coefficients.

        Example: [[b11],
                  [b22],
                  [0]]
        """

        self.__beta = new_beta
        self.__update()

    @property
    def dT(self):
        """The relative change in temperature of the lamina."""

        return self.__dT

    @dT.setter
    def dT(self, new_dT):
        """The relative change in temperature of the lamina."""

        self.__dT = new_dT
        self.__update()

    @property
    def dM(self):
        """The relative change in moisture of the lamina."""

        return self.__dM

    @dM.setter
    def dM(self, new_dM):
        """The relative change in moisture of the lamina."""

        self.__dM = new_dM
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

        z is measured from the midplane of the laminate with the convention
        of positive pointing upward."""

        return self.__z

    @z.setter
    def z(self, new_z):
        """The vertical location of the lamina's midplane.

        z is measured from the midplane of the laminate with the convention
        of positive pointing upward."""

        self.__z = new_z
        self.__update()

    @property
    def zk(self):
        """The vertical location of the lamina's top plane.

        zk is measured from the midplane of the laminate with the convention
        of positive pointing upward."""

        return self.__zk

    @zk.setter
    def zk(self, new_zk):
        """The vertical location of the lamina's top plane.

        zk is measured from the midplane of the laminate with the convention
        of positive pointing upward."""

        self.__z = new_zk - self.__t / 2
        self.__update()

    @property
    def zk1(self):
        """The vertical location of the lamina's bottom plane.

        zk1 is measured from the midplane of the laminate with the convention
        of positive pointing upward. Note: zk1 is z_(k-1)."""

        return self.__zk1

    @zk1.setter
    def zk1(self, new_zk1):
        """The vertical location of the lamina's bottom plane.

        zk1 is measured from the midplane of the laminate with the convention
        of positive pointing upward."""

        self.__z = new_zk1 + self.__t / 2
        self.__update()

    @property
    def Q(self):
        """The on-axis reduced stiffness matrix of the lamina."""

        return self.__Q

    @property
    def Q_bar(self):
        """The transformed reduced stiffness matrix of the lamina."""

        return self.__Q_bar

    @property
    def T(self):
        """The lamina transformation matrix."""

        return self.__T

    @property
    def e(self):
        """Total lamina strain"""

        return self.__e

    @property
    def e_m(self):
        """Mechanical strain in the lamina orientation."""

        return self.__e_m

    @e_m.setter
    def e_m(self, new_strain):
        """Mechanical strain in the lamina orientation."""

        self.__e_m = new_strain
        self.__update()

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

    @e_mbar.setter
    def e_mbar(self, new_strain):
        """Mechanical strain in the laminate orientation."""

        self.__e_m = mm(self.__T, new_strain)
        self.__update()

    @property
    def e_tbar(self):
        """Thermal strain in the laminate orientation."""

        return self.__e_tbar

    @property
    def e_hbar(self):
        """Mechanical strain in the lamina orientation."""

        return self.__e_hbar

    @e_hbar.setter
    def e_hbar(self, new_strain):
        """Mechanical strain in the lamina orientation."""

        self.__e_h = mm(self.__T, new_strain)
        self.__update()
