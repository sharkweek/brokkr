"""Define the Laminate class for use with mechpy.

All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions of the lamina, respectively. It is also assumed
that lamina are generally orthotropic.

See George Staab's 'Laminar Composites' for symbol conventions and
relevant equations.
"""

import numpy as np
import math
from mechpy.lam_tools.lamina import Lamina

# [ ] TODO: cleanup update() function including abstraction
# [ ] TODO: create functions to import Lamina and Laminate properties from
# Nastran BDF
# [ ] TODO: create function to calculate failure indices for Laminate
# [ ] TODO: create .is_symmetric() function
# [ ] TODO: create .is_balanced() function


class Laminate:
    """Laminate made up of multiple lamina objects.

    Note that the stacking sequence is from the bottom of the laminate
    moving upward.
    """

    def __init__(self, *plies):
        """Initialize with a list of Lamina objects"""

        # create A, B, and D matrices
        self.__A = np.zeros((3, 3))
        self.__B = np.zeros((3, 3))
        self.__D = np.zeros((3, 3))
        self.__t = 0
        self.__layup = []
        self.__N_M = np.zeros((3, 1))        # mechanical running loads
        self.__M_M = np.zeros((3, 1))        # mechanical running moments
        self.__N_T = np.zeros((3, 1))        # thermal running loads
        self.__M_T = np.zeros((3, 1))        # thermal running moments
        self.__N_H = np.zeros((3, 1))        # hygral running loads
        self.__M_H = np.zeros((3, 1))        # hygral running moments
        self.__epsilon_0 = np.zeros((3, 1))  # midplane strains
        self.__kappa_0 = np.zeros((3, 1))    # midplane curvatures
        self.dT = 0                         # temperature delta
        self.M_bar = 0                      # average moisture content
        self.__plies = plies
        self.__counter = 0

        # determine laminate properties for non-empty Laminate object
        if len(self) > 0:
            self.__update()

    # Abstracted properties
    @property
    def A(self):
        """Return the A matrix."""

        return self.__A

    @property
    def B(self):
        """Return the B matrix."""

        return self.__B

    @property
    def D(self):
        """Return the D matrix."""

        return self.__D

    @property
    def t(self):
        """Return the laminate thickness."""

        return self.__t

    @property
    def layup(self):
        """Return the layup."""


        return self.__layup

    @property
    def N_M(self):
        """Return the mechanical running load."""

        return self.__N_M

    @property
    def M_M(self):
        """Return the mechanical running moment."""

        return self.__M_M

    @property
    def N_T(self):
        """Return the thermal running load."""

        return self.__N_T

    @property
    def M_T(self):
        """Return the thermal running moment."""

        return self.__M_T

    @property
    def N_H(self):
        """Return the hygral running load."""

        return self.__N_H

    @property
    def M_H(self):
        """Return the hygral running moment."""

        return self.__M_H

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
        """Retrieve a specific item at the requested index."""

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

        # make sure all plies are fully defined Lamina objects
        for ply in self.__plies:
            self.isLamina(ply)

        self._t = 0  # reset thickness
        self._layup = []  # reset layup

        for ply in self.__plies:
            ply.ID += 1  # assign ply.ID
            ply.z = self._t + ply.tk / 2  # calculate ply.z w.r.t bottom
            self._t += ply.tk  # calculate laminate thickness
            self._layup.append(ply.theta)  # add ply orientation to layup

        # recalculate z with respect to the laminate mid-plane
        for ply in self.__plies:
            ply.z -= self._t / 2

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
        return self._layup.__repr__()

    def append(self, newPly):
        """Extended list.append() to update laminate properties on addition of
        new ply.

        Note: added plies are assumed to be placed on TOP SURFACE
        of existing laminate."""

        self.isLamina(newPly)
        self.__plies.append(newPly)
        self.__update()

    def insert(self, index, ply):
        """Insert a lamina into the laminate."""

        self.isLamina(ply)
        self.__plies.insert(index, ply)
        self.__update()

    def isLamina(ply):
        """Check if value is a ply and raises error if not."""

        if type(ply) == Lamina:
            raise TypeError("Laminates may only contain Lamina objects")
        else:
            return True

    def lamFromCSV(inputFile):
        """Determine laminate properties from input CSV file."""

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
        """Pop lamina out of Laminate."""

        self.__plies.pop(key)
        self.__update()

    def sort(self):
        """Sorts the laminate plies by ID."""

        self.__plies.sort(key=lambda ply: ply.ID)
        self.__update()

    def flip(self, renumber=False):
        """Flips the stacking sequence of the Lamina.

        Provides option to renumber the plies from the bottom up.
        """

        self.__plies = reversed(self.__plies)
        if renumber:
            _ = 0
            for each in self:
                _ += 1
                each.ID = _
            del _
        self.__update()

    def is_balanced(self):
        """Return if True if laminate is balanced."""

        pass

    def is_symmetric(self):
        """Return if True if laminate is symmetric."""

        pass
