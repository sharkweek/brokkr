"""All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions of the lamina, respectively. It is also assumed
that lamina are generally orthotropic.

See George Staab's 'Laminar Composites' for symbol conventions and
relevant equations."""

# TODO: create functions to import Lamina and Laminate properties from Nastran
# BDF
# TODO: create function to calculate failure indices for Laminate
# TODO: Update Laminate.pop(), reverse(), and sort()

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

        return np.matmul(np.matmul(np.linalg.inv(self.T_s), self.Qk), self.T_e)

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


class Laminate:
    """Laminate object made up of multiple plies of lamina

    Note that the stacking sequence is from the bottom of the laminate
    moving upward."""

    def __init__(self, plies=[]):
        """Initialize with a list of Lamina objects"""

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
        self.__plies = plies
        self.__counter = 0

        # determine laminate properties for non-empty Laminate object
        if len(self) > 0:
            self.__update()

    def __len__(self):
        """Return the number of plies in the laminate"""

        return len(self.__plies)

    def __iter__(self):
        """Returns the plies property as an iterable"""

        return iter(self.__plies)

    def __next__(self):
        """Iterates to the next ply"""

        self.__counter += 1
        if self.__counter < len(self):
            return self.__plies[self.__counter]
        else
            raise StopIteration

    def __getitem__(self, key):
        """Retrieves a specific item at the requested index"""

        return self.__plies[key]

    def __delitem__(self, key):
        """Extended the del  method to update the laminate when a ply is
        removed"""

        del self.__plies[key]
        self.__update()

    def __setitem__(self, key, newPly):
        """Sets a specific ply to a new Lamina at the requested index"""

        self.isLamina(newPly)
        self.__plies[key] = newPly
        self.__update()

    def __update(self):
        """Updates the ply and laminate attributes based on laminate stackup"""

        # Checks to make sure all plies are fully defined Lamina objects
        for ply in self.__plies:
            self.isLamina(ply)

        self.t = 0  # reset t
        self.layup = []  # reset layup

        for ply in self.__plies:
            ply.ID += 1  # assign ply.ID
            ply.z = self.t + ply.tk / 2  # calculate ply.z w.r.t bottom
            self.t += ply.tk  # calculate laminate thickness
            self.layup.append(ply.theta)  # add ply orientation to layup

        # recalculate z with respect to the laminate mid-plane
        for ply in self.__plies:
            ply.z -= self.t / 2

        # calculate A, B, and D matrices from individual ply matrix terms
        for ply in self.__plies:
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
        for ply in self.__plies:
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
        for ply in self.__plies:
            ply.epsilon_k = np.matmul(ply.T_e, (self.epsilon_0 + self.z *
                                                self.kappa_0))

            ply.sigma_k = np.matmul(ply.Qk, (ply.epsilon_k -
                                             ply.alpha_k * self.dT -
                                             ply.beta_k * self.M_bar))

    def __repr__(self):
        return self.layup.__repr__()

    def append(self, newPly):
        """Extended list.append() to update laminate properties on addition of
        new ply.

        Note: added plies are assumed to be placed on TOP SURFACE
        of existing laminate."""

        self.isLamina(newPly)
        self.__plies.append(newPly)
        self.__update()

    def flip(self, renumber=False):
        """Flips the stacking sequence of the Lamina. Provides option to
        renumber the plies from the bottom up"""

        self.__plies = reversed(self.__plies)
        if renumber:
            _ = 0
            for each in self:
                _ += 1
                each.ID = _
            del _
        self.__update()

    def insert(self, index, ply):
        """Extended list.insert() to update laminate when a ply is inserted"""

        self.isLamina(ply)
        self.__plies.insert(index, ply)
        self.__update()

    def isLamina(ply):
        """Checks if value is a ply and raises error if not"""

        if type(ply) == Lamina:
            raise TypeError("Laminates may only contain Lamina objects")
        else:
            return True

    def lamFromCSV(inputFile):
        """Determines laminate properties from input CSV file."""

        with open(inputFile, 'r') as csvLaminate:
            rawLines = csv.DictReader(csvLaminate)  # create dict from csv

            # sort by plyID - assumed order is from the bottom upward
            sortedLines = sorted(rawLines,
                                 key=lambda sortField:
                                 int(sortField['plyID']))

            for ply in sortedLines:
                self.append(Lamina(int(ply['plyID']),
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

    def pop(self, key):
        """Replicates pop() functionality for Laminate"""

        self.__plies.pop(key)
