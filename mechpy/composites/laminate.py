"""Define the Laminate class for use with mechpy.

All mechanical and thermal properties assume that x- and y-axes are in the 0-
and 90-degree directions of the lamina, respectively. It is also assumed that
lamina are generally orthotropic.

Equations and symbol conventions are per 'Mechanics Of Composite Materials' by
Robert M. Jones and 'NASA-RP-1351 BASIC MECHANICS OF LAMINATED COMPOSITE
PLATES' as noted.

TODO:
-----
* [ ] Finish lamFromCSV()
* [ ] Create is_balanced()
* [ ] Create is_symmetric()
* [ ] create functions to import Lamina and Laminate properties from
      Nastran BDF
* [ ] create function(s) to calculate failure indices for Laminate
* [ ] create .is_symmetric() function
* [ ] create .is_balanced() function
* [ ] force insert and append functions to ensure new ply has a unique
      ID
* [ ] Figure out a way to force .__update() when any ply attributes are updated
"""

import numpy as np
import csv
from .lamina import Lamina


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
        self.__t = 0                  # thickness
        self.__layup = []
        self._e_0 = np.zeros((3, 1))  # midplane strains
        self._k_0 = np.zeros((3, 1))  # midplane curvatures
        self.__dT = 0                 # change in temperature
        self.__dm = 0                 # moisture content
        self.__plies = list(plies)

        if len(self.__plies) > 0:
            # make sure all plies are Lamina objects
            for ply in self.__plies:
                self.check_lamina(ply)
            self.__update()

    # Encapsulated properties
    @property
    def A(self):
        """A matrix."""

        return self.__A

    @property
    def B(self):
        """B matrix."""

        return self.__B

    @property
    def D(self):
        """D matrix."""

        return self.__D

    @property
    def t(self):
        """Laminate thickness."""

        return self.__t

    @property
    def layup(self):
        """Laminate layup."""

        return self.__layup

    @property
    def N_m(self):
        """The mechanical running load applied to the laminate."""

        return self.__N_m

    @property
    def M_m(self):
        """The mechanical running moment applied to the laminate."""

        return self.__M_m

    @property
    def N_t(self):
        """Thermally induced running load."""

        return self.__N_t

    @property
    def M_t(self):
        """Thermally induced running moment."""

        return self.__M_t

    @property
    def N_h(self):
        """Hygroscopically induced running load."""

        return self.__N_h

    @property
    def M_h(self):
        """Hygroscopically induced running moment."""

        return self.__M_h

    # Setters for properties
    @N_m.setter
    def N_m(self, newLoad):
        """The mechanical running load applied to the laminate."""

        self.__N_m = newLoad
        self.__update()

    @M_m.setter
    def M_m(self, newLoad):
        """The mechanical running moment applied to the laminate."""

        self.__M_m = newLoad
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

        return "Laminate layup: " + self.__layup.__repr__()

    # new functions
    def __update(self):
        """Update the ply and laminate attributes based on laminate stackup."""

        # reset internal values
        self.__t = 0                    # thickness
        self.__layup = []               # layup
        self.__N_m = np.zeros((3, 1))   # mechanical running loads
        self.__M_m = np.zeros((3, 1))   # mechanical running moments
        self.__N_t = np.zeros((3, 1))   # thermal running loads
        self.__M_t = np.zeros((3, 1))   # thermal running moments
        self.__N_h = np.zeros((3, 1))   # hygroscopic running loads
        self.__M_h = np.zeros((3, 1))   # hygroscopic running moments

        for ply in self.__plies:
            ply.z = self.__t + ply.t / 2    # ply.z w.r.t bottom
            self.__t += ply.t               # total thickness
            self.__layup.append(ply.theta)  # add ply orientation to layup

        # recalculate z with respect to the laminate mid-plane
        for ply in self.__plies:
            ply.z -= self.__t / 2

        for ply in self.__plies:
            # A, B, and D matrices from individual ply matrix terms
            # Jones, Eq (4.22)-(4.24)
            self.__A += ply.Qk_bar * (ply.zk - ply.zk1)  # extension
            self.__B += ply.Qk_bar * (ply.zk**2 - ply.zk1**2)/2  # coupling
            self.__D += ply.Qk_bar * (ply.zk**3 - ply.zk1**3)/3  # bending

            # thermal running loads
            # NASA-RP-1351, Eq (92)
            ply.dT = self.__dT
            self.__N_t += np.matmul(ply.Qk_bar, ply.e_tbar)*(ply.zk - ply.zk1)
            self.__M_t += (np.matmul(ply.Qk_bar, ply.e_tbar) *
                           (ply.zk**2 - ply.zk1**2)/2)

            # hygroscopic running loads
            # NASA-RP-1351, Eq (93)
            ply.dm = self.__dm
            self.__N_h += np.matmul(ply.Qk_bar, ply.e_hbar)*(ply.zk - ply.zk1)
            self.__M_h += (np.matmul(ply.Qk_bar, ply.e_hbar) *
                           (ply.zk**2 - ply.zk1**2)/2)

        # intermediate star matrices
        # NASA-RP-1351 Eq (50)-(52a) provide a more straightforward
        # description of creating the 'prime' ABD matrix.
        A_star = np.linalg.inv(self.__A)
        B_star = np.matmul(A_star, self.__B)
        C_star = np.matmul(self.__B, A_star)
        D_star = self.__D - np.matmul(np.matmul(self.__B, A_star), self.__B)

        D_prime = np.linalg.inv(D_star)
        A_prime = A_star - np.matmul(np.matmul(B_star, D_prime), C_star)
        B_prime = np.matmul(B_star, D_prime)
        C_prime = - np.matmul(D_prime, C_star)

        ABD_prime = np.vstack((np.hstack((A_prime, B_prime)),
                               np.hstack((C_prime, D_prime))))

        # the midplane strains
        # Jones, (4.108) and (4.109). Also see Section 4.5.4
        mech_load = np.vstack((self.__N_m, self.__M_m))
        thrm_load = np.vstack((self.__N_t, self.__M_t))
        hygr_load = np.vstack((self.__N_h, self.__M_h))
        self.__e_0m, self.__k_0M = np.vsplit(np.matmul(ABD_prime, mech_load))
        self.__e_0t, self.__k_0T = np.vsplit(np.matmul(ABD_prime, thrm_load))
        self.__e_0h, self.__k_0H = np.vsplit(np.matmul(ABD_prime, hygr_load))

        # TODO: add stress and strain terms to Lamina class
        # add ply strains due to mechanical loading (in laminate orientation)
        # Jones, Eq (4.13)
        for ply in self.__plies:
            ply.e_mbar = np.matmul(ply.Tinv, self.__e_0m + ply.zk*self.__k_0M)

    # Externally accessible functions
    def append(self, newPly):
        """Add a new ply to the Laminate.

        Added plies are assumed to be placed on TOP SURFACE of existing
        laminate.
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
