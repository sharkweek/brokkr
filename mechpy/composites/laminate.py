"""Define the Laminate class for use with mechpy.

All mechanical and thermal properties assume that x- and y-axes are in the
0- and 90-degree directions of the lamina, respectively. It is also assumed
that lamina are generally orthotropic.

See George Staab's 'Laminar Composites' for symbol conventions and
relevant equations.

Equations and symbol conventions are per 'NASA-RP-1351 BASIC MECHANICS OF
LAMINATED COMPOSITE PLATES' unless otherwise noted.

https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950009349.pdf

TODO:
-----
* [ ] Update lamFromCSV()
* [ ] Update is_balanced()
* [ ] Update is_symmetric()
"""

import numpy as np
import math
from .lamina import Lamina

# [ ] TODO: cleanup update() function including abstraction
# [ ] TODO: create functions to import Lamina and Laminate properties from
# Nastran BDF
# [ ] TODO: create function to calculate failure indices for Laminate
# [ ] TODO: create .is_symmetric() function
# [ ] TODO: create .is_balanced() function
# [ ] TODO: force insert and append functions to ensure new ply has a unique
#           ID


class Laminate:
    """Laminate made up of multiple lamina objects.

    Note that the stacking sequence is from the bottom of the laminate
    moving upward.
    """

    def __init__(self, *plies):
        """Initialize with a list of Lamina objects."""

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
        self.__dT = 0                        # temperature delta
        self.__dM = 0                        # moisture content
        self.__plies = list(plies)

        if len(self) > 0:
            # make sure all plies are Lamina objects
            for ply in self.__plies:
                self.check_lamina(ply)
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
        """Return the running load applied to the laminate."""

        return self.__N_M

    @property
    def M_M(self):
        """Return the running moment applied to the laminate."""

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

    # Setters for properties
    @N_M.setter
    def N_M(self, newLoad):
        """Update when running load is changed."""

        self.__N_M = newLoad
        self.__update()

    @M_M.setter
    def M_M(self, newLoad):
        """Update when running moment is changed."""

        self.__M_M = newLoad
        self.__update()

    # Hidden functions
    def __len__(self):
        """Return the number of plies in the laminate."""

        return len(self.__plies)

    def __iter__(self):
        """Return an iterator of self.__plies."""

        return iter(self.__plies)

    def __getitem__(self, key):
        """Retrieve a specific item at the requested index."""

        return self.__plies[key]

    def __delitem__(self, key):
        """Delete a ply."""

        del self.__plies[key]
        self.__update()

    def __setitem__(self, key, new_ply):
        """Set a specific ply to a new Lamina at the requested index."""

        self.check_lamina(new_ply)
        self.__plies[key] = new_ply
        self.__update()

    def __repr__(self):
        """Set the representation of the laminate."""

        return self.__layup.__repr__()

    # new functions
    def __update(self):
        """Update the ply and laminate attributes based on laminate stackup."""

        self.__t = 0       # reset thickness
        self.__layup = []  # reset layup

        for ply in self.__plies:
            ply.z = self.__t + ply.tk / 2   # calculate ply.z w.r.t bottom
            self.__t += ply.tk              # calculate total thickness
            self.__layup.append(ply.theta)  # add ply orientation to layup

        # recalculate z with respect to the laminate mid-plane
        for ply in self.__plies:
            ply.z -= self.__t / 2

        # calculate A, B, and D matrices from individual ply matrix terms
        # NASA-RP-1351 Eq (45) through (47)
        for ply in self.__plies:
            self.__A += ply.Ak
            self.__B += ply.Bk
            self.__D += ply.Dk

        # calculate intermediate star matrices
        # NASA-RP-1351 Eq (51)
        A_star = np.linalg.inv(self.__A)
        B_star = np.matmul(A_star, self.__B)
        C_star = np.matmul(self.__B, A_star)
        D_star = self.__D - np.matmul(np.matmul(self.__B, A_star), self.__B)

        # calculate intermediate prime matrices
        # NASA-RP-1351 Eq (52) and (52a)
        D_prime = np.linalg.inv(D_star)
        A_prime = A_star - np.matmul(np.matmul(B_star, D_prime), C_star)
        B_prime = np.matmul(B_star, D_prime)
        C_prime = - np.matmul(D_prime, C_star)

        ABD_prime = np.r_[np.c_[A_prime, B_prime],
                          np.c_[C_prime, D_prime]]

        # calculate the midplane mechanical strains
        # NASA-RP-1351 Eq (52)
        mech_load = np.vstack(self.__N_M, self.__M_M)
        mech_strn = np.matmult(ABD_prime, mech_load)

        # calculate thermal and hygral ply strains
        # NASA-RP-1351 Eq (91)



        # NASA-RP-1251, Eq (43), (44), (91) through (93), and (95)
        for ply in self.__plies:
            self.__N_T += (np.matmul(ply.Qk_bar, (self.dT * ply.alpha_k)) *
                           (ply.zk - ply.zk1))  # NASA-RP-1251 Eq (91)
            self.__M_T += (np.matmul(ply.Qk_bar, (self.dT * ply.beta_k)) *
                           (ply.zk**2 - ply.zk1**2) / 2)  # NASA-RP-1251 Eq (93)
            self.__N_H += (np.matmul(ply.Qk_bar, (self.dT * ply.alpha_k)) *
                           (ply.zk - ply.zk1))  # NASA-RP-1251 Eq (91)
            self.__M_H += np.matmul(ply.Qk, ply.beta_k) * self.M_bar * \
                ply.tk * ply.z  # See Staab Eq 6.30

        # calculate total load
        # See Staab Eq 6.32, and NASA-RP-1251 Eq (94) and (95)
        N_hat = self.__N_M + self.__N_T + self.__N_H
        M_hat = self.__M_M + self.__M_T + self.__M_H

        # calculate laminate midplane strains (See Staab Eq 6.35)
        self.epsilon_0, self.kappa_0 = np.vsplit(np.matmul(ABD_prime,
                                                           np.vstack(N_hat,
                                                                     M_hat,
                                                                     2)))

        # calculate on-axis ply strains and stresses (See Staab Eq 6.37-6.38)
        for ply in self.__plies:
            ply.epsilon_k = np.matmul(ply.T, (self.epsilon_0 + self.z *
                                              self.kappa_0))

            ply.sigma_k = np.matmul(ply.Qk, (ply.epsilon_k -
                                             ply.alpha_k * self.dT -
                                             ply.beta_k * self.M_bar))

    # Externally accessible functions
    def append(self, newPly):
        """Add a new ply to the Laminate.

        Note: added plies are assumed to be placed on TOP SURFACE
        of existing laminate.
        """

        self.check_lamina(newPly)
        self.__plies.append(newPly)
        self.__update()

    def insert(self, index, ply):
        """Insert a lamina into the laminate."""

        self.check_lamina(ply)
        self.__plies.insert(index, ply)
        self.__update()

    def check_lamina(self, ply):
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
        """Sort the laminate plies by ID."""

        self.__plies.sort(key=lambda ply: ply.ID)
        self.__update()

    def flip(self, renumber=False):
        """Flip the stacking sequence of the Lamina.

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
