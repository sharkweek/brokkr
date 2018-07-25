"""Define the Lamina class for use with mechpy.

All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions of the lamina, respectively. It is also assumed
that lamina are generally orthotropic.

See George Staab's 'Laminar Composites' for symbol conventions and
relevant equations.

TODO:
* [ ] Create better intedependency between values and matrices
* [ ] provide error checking for values that don't make sense
* [ ] change default theta to be in degrees
* [ ] clean up property functions so that the values are statically encapsulated
"""

import numpy as np
import math


class Lamina:
    """Individual lamina."""

    def __init__(self,
                 tk=1,       # lamina thickness
                 theta_d=0,  # orientation of lamina
                 E11=1,      # Young's modulus in 1-direction
                 E22=1,      # Young's modulus in 2-direction
                 nu12=1,     # Poisson's ratio in 12-plane
                 G12=1,      # shear modulus in 12-plane
                 a11=1,      # coeff. of thermal expansion in 1-direction
                 a22=1,      # coeff. of thermal expansion in 2-direction
                 a12=1,      # coeff. of thermal expansion in 12-plane (shear)
                 b11=1,      # coeff. of moisture expansion in 1-direction
                 b22=1,      # coeff. of moisture expansion in 2-direction
                 b12=1):     # coeff. of moisture expansion in 12-plane (shear)
        """Initialize the Lamina instance."""

        self.ID = 0
        self.__tk = tk
        self.__z = 0
        self.__theta = math.radians(theta_d)  # theta is assumed to be in degrees
        self.__E11 = E11
        self.__E22 = E22
        self.__nu12 = nu12
        self.__G12 = G12
        self.__alpha_k = np.array([[a11], [a22], [a12]], dtype=float)
        self.__beta_k = np.array([[b11], [b22], [b12]], dtype=float)

        self.epsilon_k = np.array((3, 1))  # on-axis ply strains

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
    def E11(self):
        """Return Young's modulus in the 1-direction."""

        return self.__E11

    @property
    def E22(self):
        """Return E22 for the lamina."""

        return self.__E22

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
        """Return the on-axis reduced stiffness matrix. """

        return self.__Qk

    @property
    def Qk_bar(self):
        """Return the transformed reduced stiffness matrix."""

        return self.__Qk_bar

    @property
    def T_e(self):
        """Return the strain transformation matrix."""

        return self.__Te

    @property
    def T_s(self):
        """Return the stress transformation matrix."""

        return self.__Ts

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
        elif self.__E11 == 0 or math.isnan(self.__E11):
            raise TypeError("lamina.E11 must be a non-zero number")
        elif self.E22 == 0 or math.isnan(self.E22):
            raise TypeError("lamina.E22 must be a non-zero number")
        elif (self.nu12 <= 0 or self.nu12 > 1) or math.isnan(self.nu12):
            raise TypeError("""lamina.nu12 must be a non-zero number between
                            zero and 1""")
        elif self.G12 == 0 or math.isnan(self.G12):
            raise TypeError("lamina.G12 must be a non-zero number")
        else:
            return True

    @theta.setter
    def theta(self, new_theta):
        """Update Lamina when theta is changed."""

        self.__theta = new_theta
        self._update()

    @theta_d.setter
    def theta_d(self, new_theta_d):
        """Assign new theta if theta_d is changed."""

        self.__theta = math.radians(new_theta_d)

    @z.setter
    def z(self, new_z):
        """Update when a new z value is assigned."""

        self.__z = new_z
        self._update()

    @E11.setter
    def E11(self, new_E11):
        """Return Young's modulus in the 1-direction."""

        self.__E11 = new_E11
        self._update()

    @E22.setter
    def E22(self, new_E22):
        """Return E22 for the lamina."""

        self.__E22 = new_E22
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

        self.__zk = self.z + self.__tk / 2
        self.__zk1 = self.z - self.__tk / 2

        # update on-axis reduced stiffness matrix
        # See Staab section 3.2.3.1, Eq 3.9.
        nu21 = self.nu12 * self.E22 / self.__E11  # Staab Eq 3.2
        q11 = self.__E11 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E22 / (1 - self.nu12 * nu21)
        q22 = self.E22 / (1 - self.nu12 * nu21)
        q66 = self.G12

        self.__Qk = np.array([[q11, q12, 0],
                              [q12, q22, 0],
                              [0, 0, q66]])

        # update the strain transformation matrix
        # Staab section 2.2.1, Eq 2.1.
        m = math.cos(self.__theta)
        n = math.sin(self.__theta)

        self.__Te = np.array([[m**2, n**2, m*n],
                              [n**2, m**2, -m*n],
                              [-2*m*n, 2*m*n, m**2 - n**2]])

        # update the stress transformation matrix
        # See Staab section 2.3, Eq 2.3.
        # The value of Ts should always be equal to Te
        m = math.cos(self.__theta)
        n = math.sin(self.__theta)

        self.__Ts = np.array([[m**2, n**2, 2*m*n],
                         [n**2, m**2, -2*m*n],
                         [-m*n, m*n, m**2 - n**2]])

        # update the transformed reduced stiffness matrix
        # See Staab section 3.2.2.
        self.__Qk_bar = np.matmul(np.matmul(np.linalg.inv(self.T_s),
                                            self.Qk),
                                  self.T_e)

        # the extensional stiffness matrix (in a laminate)
        self.__Ak = self.Qk_bar * (self.zk - self.zk1)

        # the coupling matrix (in a lamiate)
        self.__Bk = (1 / 2) * self.Qk_bar * (self.zk**2 - self.zk1**2)

        # the bending matrix (in a laminate)
        self.__Dk = (1 / 3) * self.Qk_bar * (self.zk**3 - self.zk1**3)
