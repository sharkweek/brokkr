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
                 tk=0,     # lamina thickness
                 theta=0,  # orientation of lamina
                 E11=0,    # Young's modulus in 1-direction
                 E22=0,    # Young's modulus in 2-direction
                 nu12=0,   # Poisson's ratio in 12-plane
                 G12=0,    # shear modulus in 12-plane
                 a11=0,    # coeff. of thermal expansion in 1-direction
                 a22=0,    # coeff. of thermal expansion in 2-direction
                 a12=0,    # coeff. of thermal expansion in 12-plane (shear)
                 b11=0,    # coeff. of moisture expansion in 1-direction
                 b22=0,    # coeff. of moisture expansion in 2-direction
                 b12=0):   # coeff. of moisture expansion in 12-plane (shear)
        """Initialize the Lamina instance."""

        self.ID = 0
        self.tk = tk
        self.z = 0
        self.theta = 0  # theta is assumed to be in degrees
        self.E11 = E11
        self.E22 = E22
        self.nu12 = nu12
        self.G12 = G12
        self.alpha_k = np.array([[a11], [a22], [a12]], dtype=float)
        self.beta_k = np.array([[b11], [b22], [b12]], dtype=float)

        self.epsilon_k = np.array((3, 1))  # on-axis ply strains

    @property
    def theta_r(self):
        """Return theta in radians."""

        return math.radians(self.theta)

    @theta_r.setter
    def theta_r(self, new_theta_r):
        """Assign new theta if theta_r is changed."""

        self.theta = math.degrees(new_theta_r)

    @property
    def zk(self):
        """Return the upper boundary of the lamina."""

        return self.z + self.tk / 2

    @zk.setter
    def zk(self, new_zk):
        """Correct z if zk is changed."""

        self.z = new_zk - self.tk / 2

    @property
    def zk1(self):
        """Return the lower boundary of the lamina."""

        return self.z - self.tk / 2

    @zk1.setter
    def zk1(self, new_zk1):
        """Correct z if zk1 is changed."""

        self.z = new_zk1 + self.tk / 2

    @property
    def Qk(self):
        """Return the on-axis reduced stiffness matrix of the lamina.

        The lamina assumes plane stress only.

        See Staab section 3.2.3.1, Eq 3.9.
        """

        nu21 = self.nu12 * self.E22 / self.E11  # Staab Eq 3.2
        q11 = self.E11 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E22 / (1 - self.nu12 * nu21)
        q22 = self.E22 / (1 - self.nu12 * nu21)
        q66 = self.G12

        return np.array([[q11, q12, 0],
                         [q12, q22, 0],
                         [0, 0, q66]])

    @property
    def Qk_bar(self):
        """Return the transformed reduced stiffness matrix of the lamina.

        See Staab section 3.2.2.
        """

        return np.matmul(np.matmul(np.linalg.inv(self.T_s), self.Qk), self.T_e)

    @property
    def T_e(self):
        """Return the strain transformation matrix of the lamina.

        See Staab section 2.2.1, Eq 2.1.
        """

        m = math.cos(self.theta_r)
        n = math.sin(self.theta_r)

        return np.array([[m**2, n**2, m*n],
                         [n**2, m**2, -m*n],
                         [-2*m*n, 2*m*n, m**2 - n**2]])

    @property
    def T_s(self):  # should always be equal to T_e
        """Return the stress transformation matrix of the lamina.

        See Staab section 2.3, Eq 2.3.
        """

        m = math.cos(self.theta_r)
        n = math.sin(self.theta_r)

        return np.array([[m**2, n**2, 2*m*n],
                         [n**2, m**2, -2*m*n],
                         [-m*n, m*n, m**2 - n**2]])

    @property
    def Ak(self):
        """Calculate the A matrix for the lamina within it's laminate."""

        return self.Qk_bar * (self.zk - self.zk1)

    @property
    def Bk(self):
        """Calculate the B matrix for the lamina within it's laminate."""

        return (1 / 2) * self.Qk_bar * (self.zk**2 - self.zk1**2)

    @property
    def Dk(self):
        """Calculate the D matrix for the lamina within it's laminate."""

        return (1 / 3) * self.Qk_bar * (self.zk**3 - self.zk1**3)

    @property
    def isFullyDefined(self):
        """Check if lamina is fully defined with proper attr data types."""

        if self.tk == 0 or math.isnan(self.tk):
            raise TypeError("lamina.tk must be a non-zero number")
        elif math.isnan(self.theta):
            raise TypeError("lamina.theta must be a number (in degrees)")
        elif self.E11 == 0 or math.isnan(self.E11):
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
