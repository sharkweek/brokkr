"""All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions, respectively

See George Staab's 'Laminar Composites' for symbol conventions and
relevant equations."""

import numpy as np
import math
import csv


class Lamina:
    """Individual lamina"""

    def __init__(self,
                 tk=0,
                 theta=0,
                 E11=0,
                 E22=0,
                 nu12=0,
                 G12=0,
                 a11=0,
                 a22=0,
                 b11=0,
                 b22=0):
        self.ID = 0
        self.tk = tk
        self.z = 0
        self.theta = 0  # theta is assumed to be in degrees
        self.E11 = E11
        self.E22 = E22
        self.nu12 = nu12
        self.G12 = G12
        self.alpha_k = np.array([[a11], [a22], [0]], dtype=float)
        self.beta_k = np.array([[b11], [b22], [0]], dtype=float)

        self.epsilon_k = np.array((3, 1))  # on-axis ply strains

    @property
    def theta_r(self):
        """self.theta in radians"""

        return math.radians(self.theta)

    @theta_r.setter
    def theta_r(self, new_theta_r):
        """Corrects self.theta if self.theta_r is changed"""

        self.theta = math.degrees(new_theta_r)

    @property
    def zk(self):
        """The upper boundary of the lamina"""

        return self.z + self.tk / 2

    @zk.setter
    def zk(self, new_zk):
        """Corrects z if zk is changed"""

        self.z = new_zk - self.tk / 2

    @property
    def zk1(self):
        """The lower boundary of the lamina"""

        return self.z - self.tk / 2

    @zk1.setter
    def zk1(self, new_zk1):
        """Corrects z if zk1 is changed"""

        self.z = new_zk1 + self.tk / 2

    @property
    def Qk(self):
        """The on-axis reduced stiffness matrix of the lamina in the lamina
        assuming plane stress only.

        See Staab section 3.2.3.1, Eq 3.9"""

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
        """calculates the transformed reduced stiffness matrix of the lamina

        See Staab section 3.2.2."""

        T = self.T_e
        T_inv = np.linalg.inv(T)
        Q = self.Qk

        _ = np.matmul(T_inv, Q)

        return np.matmul(_, T)

    @property
    def T_e(self):
        """The strain transformation matrix of the lamina for a given theta

        See Staab section 2.2.1, Eq 2.1"""

        m = math.cos(self.theta_r)
        n = math.sin(self.theta_r)

        return np.array([[m**2, n**2, m*n],
                         [n**2, m**2, -m*n],
                         [-2*m*n, 2*m*n, m**2 - n**2]])

    @property
    def T_s(self):  # should always be equal to T_e
        """The stress transformation matrix of the lamina for a given theta

        See Staab section 2.3, Eq 2.3"""

        m = math.cos(self.theta_r)
        n = math.sin(self.theta_r)

        return np.array([[m**2, n**2, 2*m*n],
                         [n**2, m**2, -2*m*n],
                         [-m*n, m*n, m**2 - n**2]])

    @property
    def Ak(self):
        """Calculate the A matrix for the lamina within it's laminate"""

        return self.Qk_bar * (self.zk - self.zk1)

    @property
    def Bk(self):
        """Calculate the B matrix for the lamina within it's laminate"""

        return (1 / 2) * self.Qk_bar * (self.zk**2 - self.zk1**2)

    @property
    def Dk(self):
        """Calculate the D matrix for the lamina within it's laminate"""

        return (1 / 3) * self.Qk_bar * (self.zk**3 - self.zk1**3)

    @property
    def isFullyDefined(self):
        """Checks if lamina is fully defined with proper attribute data
        types"""

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


class Laminate(list):
    """Laminate list object made up of multiple plies of lamina

    Note that the stacking sequence is from the bottom of the laminate
    moving upward."""

    def __init__(self):
        """Initialize with a list of Lamina objects"""

        # preserve list __init__() and append with additional laminate
        # properties
        super().__init__()

        # create A, B, and D matrices
        self.A = np.zeros((3, 3))
        self.B = np.zeros((3, 3))
        self.D = np.zeros((3, 3))
        self.t = 0
        self.layup = []
        self.N_M = np.zeros((3, 1))  # mechanical running loads
        self.M_M = np.zeros((3, 1))  # mechanical running moments
        self._N_T = np.zeros((3, 1))  # thermal running loads
        self._M_T = np.zeros((3, 1))  # thermal running moments
        self._N_H = np.zeros((3, 1))  # hygral running loads
        self._M_H = np.zeros((3, 1))  # hygral running moments
        self.epsilon_0 = np.zeros((3, 1))  # midplane strains
        self.kappa_0 = np.zeros((3, 1))  # midplane curvatures
        self.dT = 0  # temperature delta
        self.M_bar = 0  # average moisture content

        # determine laminate properties for non-empty Laminate object
        if len(self) > 0:
            self.__lamUpdate_()

    def append(self, newPly):
        """Extended list.append() to update laminate properties on addition of
        new ply.

        Note: added plies are assumed to be placed on TOP SURFACE
        of existing laminate."""

        # preserve list.append() and extend to update effective laminate
        # properties when a new ply is added
        super().append(newPly)
        self.__lamUpdate_()

    def remove(self, ply):
        """Extended the list.remove() method to update the laminate when a
        ply is removed"""

        super().remove(ply)
        self.__lamUpdate_()

    def insert(self, ply):
        """Extended list.insert() to update laminate when a ply is inserted"""

        super().insert(ply)
        self.__lamUpdate_()

    def __lamUpdate_(self):
        """Updates the ply and laminate attributes based on laminate stackup"""

        # Checks to make sure all plies are fully defined Lamina objects
        for ply in self:
            if type(ply) != Lamina:
                raise TypeError("Laminates may only contain Lamina objects")
            elif ply.isFullyDefined:
                pass

        self.t = 0  # reset t
        self.layup = []  # reset layup

        for ply in self:
            ply.ID += 1  # assign ply.ID
            ply.z = self.t + ply.tk / 2  # calculate ply.z w.r.t bottom
            self.t += ply.tk  # calculate laminate thickness
            self.layup.append(ply.theta)  # add ply orientation to layup

        # recalculate z with respect to the laminate mid-plane
        for ply in self:
            ply.z -= self.t / 2

        # calculate A, B, and D matrices from individual ply matrix terms
        for ply in self:
            self.A += ply.Ak
            self.B += ply.Bk
            self.D += ply.Dk

        # calculate intermediate star and prime matrices
        A_star = np.linalg.inv(self.A)
        B_star = np.matmul(self.A_star, self.B)
        C_star = np.matmul(self.B, self.A_star)
        D_star = self.D - np.matmul(np.matmul(self.B, np.matmul(self.A_star,
                                                                self.B)))

        A_prime = A_star - np.matmul(B_star, np.matmul(np.linalg.inv(D_star),
                                                       C_star))
        B_prime = np.matmul(B_star, np.linalg.inv(D_star))
        C_prime = - np.matmul(np.linalg.inv(D_star), C_star)
        D_prime = np.linalg.inv(D_star)

        ABD_prime = np.r_[np.c_[A_prime, B_prime],
                          np.c_[C_prime, D_prime]]

        # calculate thermal and hygral loads
        for ply in self:
            self._N_T += np.matmul(ply.Qk_bar, ply.alpha_k) * self.dT * \
                ply.tk  # See Staab Eq 6.25
            self._M_T += np.matmul(ply.Qk_bar, ply.alpha_k) * self.dT * \
                ply.tk * ply.z  # See Staab Eq 6.26
            self._N_H += np.matmul(ply.Qk, ply.beta_k) * self.M_bar * \
                ply.tk  # See Staab Eq 6.29
            self._M_H += np.matmul(ply.Qk, ply.beta_k) * self.M_bar * \
                ply.tk * ply.z  # See Staab Eq 6.30

        # calculate total load (See Staab Eq 6.32)
        _N_hat = self.N_M + self._N_T + self._N_H
        _M_hat = self.M_M + self._M_T + self._M_H

        # calculate laminate midplane strains (See Staab Eq 6.35)
        self.epsilon_0, self.kappa_0 = np.vsplit(np.matmul(ABD_prime,
                                                           np.vstack(_N_hat,
                                                                     _M_hat,
                                                                     2)))

        # calculate on-axis ply strains and stresses (See Staab Eq 6.37-6.38)
        for ply in self:
            ply.epsilon_k = np.matmul(ply.T_e, (self.epsilon_0 + self.z *
                                                self.kappa_0))

            ply.sigma_k = np.matmul(ply.Qk, (ply.epsilon_k -
                                             ply.alpha_k * self.dT -
                                             ply.beta_k * self.M_bar))

    def __repr__(self):
        return "Number of plies: %s\n" \
               "t: %s\n" \
               "Layup: %s\n" \
               % (len(self),
                  self.t,
                  self.layup)

    def lamFromCSV(inputFile):
        """Determines laminate properties from input CSV file."""

        with open(inputFile, 'r') as csvLaminate:
            rawLines = csv.DictReader(csvLaminate)  # create dict from csv

            # sort by plyID - assumed order is from the bottom upward
            sortedLines = sorted(rawLines,
                                 key=lambda sortField:
                                 int(sortField['plyID']))

            L = Laminate()
            for ply in sortedLines:
                L.append(Lamina(int(ply['plyID']),
                                float(ply['thick']),
                                float(ply['theta']),
                                float(ply['E11']),
                                float(ply['E22']),
                                float(ply['nu12']),
                                float(ply['G12']),
                                float(ply['a11']),
                                float(ply['a22']),
                                float(ply['b11']),
                                float(ply['b22'])))
