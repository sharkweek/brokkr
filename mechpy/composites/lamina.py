"""Defines the Lamina class for use with mechpy.

All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions of the lamina, respectively. It is also assumed
that lamina are generally orthotropic.

Equations and symbol conventions are per 'NASA-RP-1351 BASIC MECHANICS OF
LAMINATED COMPOSITE PLATES' unless otherwise noted.

https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950009349.pdf

TODO:
-----
* [ ] add a property for the compliance matrix Sk
* [ ] further implement the PlyVector class
* [ ] move moduli into
"""

import numpy as np
import math


class Lamina:
    """An individual lamina or "ply".

    Each lamina is assumed to be orthotropic.
    """

    def __init__(self,
                 tk=1,       # lamina thickness
                 theta_d=0,  # orientation of lamina
                 E1=1,      # Young's modulus in 1-direction
                 E2=1,      # Young's modulus in 2-direction
                 nu12=1,     # Poisson's ratio in 12-plane
                 G12=1,      # shear modulus in 12-plane
                 a11=1,      # coeff. of thermal expansion in 1-direction
                 a22=1,      # coeff. of thermal expansion in 2-direction
                 a12=1,      # coeff. of thermal expansion in 12-plane (shear)
                 b11=1,      # coeff. of moisture expansion in 1-direction
                 b22=1,      # coeff. of moisture expansion in 2-direction
                 b12=1,      # coeff. of moisture expansion in 12-plane (shear)
                 dT=0,       # change in temperature
                 dM=0):      # moisture absorption
        """Initialize the Lamina instance."""

        self.ID = 0
        self.__tk = tk
        self.__z = 0
        self.__theta = math.radians(theta_d)
        self.__E1 = E1
        self.__E2 = E2
        self.__nu12 = nu12
        self.__G12 = G12
        self.__alpha_k = np.array([[a11], [a22], [a12]], dtype=float)
        self.__beta_k = np.array([[b11], [b22], [b12]], dtype=float)
        self.__dT = dT
        self.__dM = dM

        self.__epsilon = {'mechanical': {'lamina': np.zeros((3, 1)),
                                         'laminate': np.zeros((3, 1))},
                          'thermal': {'lamina': np.zeros((3, 1)),
                                      'laminate': np.zeros((3, 1))},
                          'hygral': {'lamina': np.zeros((3, 1)),
                                     'laminate': np.zeros((3, 1))}
                          }

        self._update()

    @property
    def theta(self):
        """Return the lamina orientation in radians."""

        return self.__theta

    @property
    def theta_d(self):
        """Return theta in degrees."""

        return math.degrees(self.__theta)

    @property
    def E1(self):
        """Return Young's modulus in the 1-direction."""

        return self.__E1

    @property
    def E2(self):
        """Return E2 for the lamina."""

        return self.__E2

    @property
    def nu12(self):
        """Return nu12 for the lamina."""

        return self.__nu12

    @property
    def G12(self):
        """Return G12 for the lamina."""

        return self.__G12

    @property
    def alpha_k(self):
        """Return alpha_k for the lamina."""

        return self.__alpha_k

    @property
    def beta_k(self):
        """Return beta_k for the lamina."""

        return self.__beta_k

    @property
    def alpha_k_bar(self):
        """Return the transformed CTE matrix for the lamina."""

        return self.__alpha_k_bar

    @property
    def beta_k_bar(self):
        """Return transformed hygral matrix for the lamina."""

        return self.__beta_k_bar

    @property
    def dT(self):
        """Return the change in temperature of the lamina."""

        return self.__dT

    @property
    def dM(self):
        """Return the moisture absorption of the lamina."""

        return self.__dM

    @property
    def z(self):
        """Return the midplane location of the lamina."""

        return self.__z

    @property
    def zk(self):
        """Return the upper boundary of the lamina."""

        return self.__zk

    @property
    def zk1(self):
        """Return the lower boundary of the lamina."""

        return self.__zk1

    @property
    def Qk(self):
        """Return the on-axis reduced stiffness matrix."""

        return self.__Qk

    @property
    def Qk_bar(self):
        """Return the transformed reduced stiffness matrix."""

        return self.__Qk_bar

    @property
    def T(self):
        """Return the strain transformation matrix."""

        return self.__T

    @property
    def Tinv(self):
        """Return the inverse of the strain transformation matrix."""

        return self.__Tinv

    @property
    def Ak(self):
        """Return the A matrix for the lamina within it's laminate."""

        return self.__Ak

    @property
    def Bk(self):
        """Return the B matrix for the lamina within it's laminate."""

        return self.__Bk

    @property
    def Dk(self):
        """Return the D matrix for the lamina within it's laminate."""

        return self.__Dk

    @property
    def isFullyDefined(self):
        """Check if lamina is fully defined with proper attr data types."""

        if self.__tk == 0 or math.isnan(self.__tk):
            raise TypeError("lamina.tk must be a non-zero number")
        elif math.isnan(self.__theta):
            raise TypeError("lamina.__theta must be a number (in degrees)")
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

    @property
    def epsilon(self, effect, c_sys):
        """Return the lamina strain matrix in the specified c-system.

        effect: 'mechanical', 'thermal' or 'hygral'
        c_sys: 'lamina' or 'laminate'
        """




    @theta.setter
    def theta(self, new_theta):
        """Update Lamina when theta is changed."""

        self.__theta = new_theta
        self._update()

    @theta_d.setter
    def theta_d(self, new_theta_d):
        """Assign new theta if theta_d is changed."""

        self.__theta = math.radians(new_theta_d)

    @E1.setter
    def E1(self, new_E1):
        """Return Young's modulus in the 1-direction."""

        self.__E1 = new_E1
        self._update()

    @E2.setter
    def E2(self, new_E2):
        """Return E2 for the lamina."""

        self.__E2 = new_E2
        self._update()

    @nu12.setter
    def nu12(self, new_nu12):
        """Return nu12 for the lamina."""

        self.__nu12 = new_nu12
        self._update()

    @G12.setter
    def G12(self, new_G12):
        """Return G12 for the lamina."""

        self.__G12 = new_G12
        self._update()

    @alpha_k.setter
    def alpha_k(self, new_alpha_k):
        """Return alpha_k for the lamina."""

        self.__alpha_k = new_alpha_k
        self._update()

    @beta_k.setter
    def beta_k(self, new_beta_k):
        """Return beta_k for the lamina."""

        self.__beta_k = new_beta_k
        self._update()

    @dT.setter
    def dT(self, new_dT):
        """Update the change in temperature of the lamina."""

        return self.__dT

    @dM.setter
    def dM(self, new_dM):
        """Update the moisture absorption of the lamina."""

        return self.__dM

    @z.setter
    def z(self, new_z):
        """Update when a new z value is assigned."""

        self.__z = new_z
        self._update()

    @zk.setter
    def zk(self, new_zk):
        """Update when a new zk value is assigned."""

        self.__z = new_zk - self.__tk / 2
        self._update()

    @zk1.setter
    def zk1(self, new_zk1):
        """Correct z if zk1 is changed."""

        self.__z = new_zk1 + self.__tk / 2
        self._update()

    def __update(self):
        """Update calculated properties when new values are assigned."""

        # upper and lower boundaries of lamina
        self.__zk = self.z + self.__tk / 2
        self.__zk1 = self.z - self.__tk / 2

        # on-axis reduced stiffness matrix, Qk
        # NASA-RP-1351 Eq (9) and (10)
        nu21 = self.nu12 * self.E2 / self.__E1  # NASA-RP-1351, Eq (6)
        q11 = self.__E1 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E2 / (1 - self.nu12 * nu21)
        q22 = self.E2 / (1 - self.nu12 * nu21)
        q66 = self.G12

        self.__Qk = np.array([[q11, q12, 0],
                              [q12, q22, 0],
                              [0, 0, q66]])

        # the transformation matrix and its inverse
        # NASA-RP-1351 Eq (15) and (16)
        m = math.cos(self.__theta)
        n = math.sin(self.__theta)

        self.__T = np.array([[m**2, n**2, 2*m*n],
                             [n**2, m**2, -2*m*n],
                             [-m*n, m*n, m**2 - n**2]])

        self.__Tinv = np.linalg.inv(self.__T)

        # the transformed reduced stiffness matrix (laminate coordinate system)
        # NASA-RP-1351 Eq (21)
        imat = np.array([[1, 0, 0],
                         [0, 1, 0],
                         [0, 0, 2]])

        self.__Qk_bar = np.matmul(np.matmul(np.matmul(self.__Tinv, self.__Qk),
                                            imat), self.__T)

        # the transformed CTE and hygral matrices (laminate coordinate system)
        self.__alpha_k_bar = np.matmul(self.__Tinv, self.__alpha_k)
        self.__beta_k_bar = np.matmul(self.__Tinv, self.__beta_k)

        # thermal and hygrla strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.__epsilon_0T = self.__alpha_k_bar * self.__dT
        self.__epsilon_kT = np.matmul(self.__T, self.__epsilon_0T)
        self.__epsilon_0H = self.__beta_k_bar * self.__dM
        self.__epsilon_kH = np.matmul(self.__T, self.__epsilon_0H)

        # the extensional stiffness matrix (laminate coordinate system)
        # NASA-RP-1351 Eq (45)
        self.__Ak = self.__Qk_bar * (self.__zk - self.__zk1)

        # the coupling matrix (laminate coordinate system)
        # NASA-RP-1351 Eq (46)
        self.__Bk = (1 / 2) * self.__Qk_bar * (self.__zk**2 - self.__zk1**2)

        # the bending matrix (laminate coordinate system)
        # NASA-RP-1351 Eq (47)
        self.__Dk = (1 / 3) * self.__Qk_bar * (self.__zk**3 - self.__zk1**3)
