"""Define the Laminate class for use with mechpy.

Assumptions
-----------
* Lamina are generally orthotropic.
* The laminate z-direction is upward, away from the bottom surface.
* Reported strain values are in engineering strain.
* All loads and moments are running values supplied as force or moment per unit
  width.
* Unit systems are consistent (i.e. SI, US, or Imperial).
* Equations and symbol conventions are per NASA-RP-1351 'Basic Mechanics of
  Laminated Composite Plates' and Jones' 'Mechanics Of Composite Materials'.

TODO:
-----
* [ ] Create is_balanced() method
* [ ] create functions to import Lamina and Laminate properties from Nastran
      BDF
* [ ] create methods to calculate failure indices for Laminate
* [ ] figure out a way to update Laminate when a Lamina property is changed
* [ ] fix calculations for effective laminate properties
* [ ] reorganize properties and setters to be adjacent
"""

from .lamina import Lamina, Ply
import numpy as np
from numpy import matmul as mm
import pandas as pd


class Laminate:
    """
    Laminate(*plies)

    A composite laminate made up of multiple orthotropic lamina objects.

    Parameters
    ----------
    *plies : mechpy.composites.Lamina or array_like of
            mechpy.composites.Laminate or string
        lamina the make up the laminate ordered from the bottom surface to the
        top surface. In creating the Laminate, either a sequence of Lamina
        objects or a file path to a CSV file containing all the appropriate
        data can be supplied.

    Attributes
    ----------
    t : float
        total thickness of the laminate
    layup : list of int or list of float
        laminate stacking sequence
    dT : int or float
        change in temperature
    dM : int or float
        change in moisture
    N_m : numpy.array_like
        mechanical running load
    N_t : numpy.array_like
        thermal running load
    N_h : numpy.array_like
        hygral (moisture) running load
    M_m : numpy.array_like
        mechanical running moment
    M_t : numpy.array_like
        thermal running moment
    M_h : numpy.array_like
        hygral (moisture) running moment
    E1_eff : float
        effective laminate modulus in the x-direction
    E2_eff : float
        effective laminate modulus in the y-direction
    G12_eff : float
        effective shear modulus in the xy-plane
    e_0 : numpy.array_like
        total laminate midplane strain
    e_0m : numpy.array_like
        laminate midplane strain due to mechanical loading
    e_0t : numpy.array_like
        laminate midplane strain due to thermal effects
    e_0h : numpy.array_like
        laminate midplane strain due to hygral effects
    k_0 : numpy.array_like
        total laminate midplane curvature
    k_0m : numpy.array_like
        laminate midplane curvature due to mechanical loading
    k_0t : numpy.array_like
        laminate midplane curvature due to thermal effects
    k_0h : numpy.array_like
        laminate midplane curvature due to hygral effects

    Methods
    -------

    """

    __name__ = 'Laminate'
    __protected__ = ['t', 'layup', 'N_t', 'N_h']
    __unprotected__ = ['_Laminate__plies']
    __updated__ = ['dT', 'dM', 'N_m', 'M_m']
    __slots__ = __protected__ + __unprotected__ + __updated__

    def __init__(self, *plies):
        """Initialize with a list of Lamina objects.

        Accepts a series of Lamina objects or a single string calling out the
        filepath to a csv file containg all the laminate data.
        """

        # create and zero out all properties needed for .__update()
        self.clear()
        self.__plies = list(plies)

        # import from csv if filepath is provided,otherwise, process as laminae
        try:
            self.from_csv(plies[0])

        except TypeError:
            if len(self.__plies) > 0:
                # make sure all plies are Lamina objects
                for ply in self.__plies:
                    self.check_lamina(ply)

                self.__update()

    def __setattr__(self, attr, val):
        """Check and set attributes."""

        # udpate laminate after updated properties are set
        if attr in self.__updated__:
            self.__sattr(attr, val)
            self.__update()

        # check if attribute is an unprotected value
        elif attr in self.__unprotected__:
            self.__sattr(attr, val)

        # don't set protected values
        elif attr in self.__protected__:
            raise AttributeError(self.__name__ + ".%s" % attr
                                 + "is a derived value and cannot be set")

    def __sattr(self, arg, val):
        """Directly set attribute."""
        super().__setattr__(arg, val)

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

    def __update(self):
        """Update the ply and laminate attributes based on laminate stackup."""

        # reset internal values
        self.__t = sum([ply.t for ply in self.__plies])
        self.__layup = [ply.theta for ply in self.__plies]
        self.__N_t = np.zeros((3, 1))   # thermal running loads
        self.__M_t = np.zeros((3, 1))   # thermal running moments
        self.__N_h = np.zeros((3, 1))   # hygroscopic running loads
        self.__M_h = np.zeros((3, 1))   # hygroscopic running moments

        # calculate ply bottom planes
        zk1 = self.__t / 2
        for i in range(1, int(len(self.__plies) + 1)):
            zk1 -= self.__plies[-i].t
            self.__plies[-i].zk1 = zk1

        # calculate the A, B, and D matricies; thermal and hygral loads
        self.__A = np.zeros((3, 3))
        self.__B = np.zeros((3, 3))
        self.__D = np.zeros((3, 3))

        for ply in self.__plies:
            # NASA-RP-1351, Eq (45) through (47)
            self.__A += ply.Qbar * (ply.t)
            self.__B += (1/2) * ply.Qbar * (ply.t * ply.z)
            self.__D += (1/3) * ply.Qbar * (ply.t**3 / 12 + ply.t * ply.z**2)

            # thermal running loads
            # NASA-RP-1351, Eq (92)
            ply.dT = self.__dT
            self.__N_t += mm(ply.Qbar, ply.e_tbar)*(ply.zk - ply.zk1)
            self.__M_t += (mm(ply.Qbar, ply.e_tbar)
                           * (ply.zk**2 - ply.zk1**2)/2)

            # hygroscopic running loads
            # NASA-RP-1351, Eq (93)
            ply.dM = self.__dM
            self.__N_h += mm(ply.Qbar, ply.e_hbar)*(ply.zk - ply.zk1)
            self.__M_h += (mm(ply.Qbar, ply.e_hbar)
                           * (ply.zk**2 - ply.zk1**2)/2)

        # intermediate star matrices
        # NASA-RP-1351 Eq (50)-(52a)
        A_star = np.linalg.inv(self.__A)
        B_star = - mm(A_star, self.__B)
        C_star = mm(self.__B, A_star)
        D_star = self.__D - mm(mm(self.__B, A_star), self.__B)

        D_prime = np.linalg.inv(D_star)
        C_prime = - mm(D_prime, C_star)
        B_prime = mm(B_star, D_prime)
        A_prime = A_star - mm(mm(B_star, D_prime), C_star)

        ABD_prime = np.vstack((np.hstack((A_prime, B_prime)),
                               np.hstack((C_prime, D_prime))))

        # the midplane strains
        # Jones, (4.108) and (4.109). Also see Section 4.5.4
        mech_load = np.vstack((self.__N_m, self.__M_m))
        thrm_load = np.vstack((self.__N_t, self.__M_t))
        hygr_load = np.vstack((self.__N_h, self.__M_h))
        self.__e_0m, self.__k_0m = np.vsplit(mm(ABD_prime, mech_load), 2)
        self.__e_0t, self.__k_0t = np.vsplit(mm(ABD_prime, thrm_load), 2)
        self.__e_0h, self.__k_0h = np.vsplit(mm(ABD_prime, hygr_load), 2)

        # add ply strains due to mechanical loading (in laminate orientation)
        # Jones, Eq (4.13)
        self.__e_0 = self.__e_0m + self.__e_0t + self.__e_0h
        self.__k_0 = self.__k_0m + self.__k_0t + self.__k_0h

        for ply in self.__plies:
            ply.e_mbar = mm(ply.T, (self.__e_0 + ply.z * self.__k_0))

        # calculate effective laminate properties
        # NASA, Section 5
        # NOTE: This is only valuable for symmetric laminates.
        if self.is_symmetric():
            # NASA, Eq. 64, 70, and 76
            A = self.__A
            self.__E1_eff = (
                A[0,0]
                + A[0,1] * ((A[1,2]*A[0,2] - A[0,1]*A[2,2])
                            / (A[1,1]*A[2,2] - A[1,2]**2))
                + A[0,2] * (-A[0,2]/A[2,2]
                            + ((A[1,2]*A[0,1]*A[2,2] - A[1,2]**2 * A[0,2])
                                / (A[1,1]*A[2,2]**2 - A[1,2]**2 * A[2,2])))
            ) / self.__t

            self.__E2_eff = (
                A[0,1] * ((A[0,2]*A[1,2] - A[0,1]*A[2,2])
                        / (A[0,0]*A[2,2] - A[0,2]**2))
                + A[1,1]
                + A[1,2] * (-A[1,2]/A[2,2]
                            + ((A[0,2]*A[0,1]*A[2,2] - A[0,2]**2 * A[1,2])
                                / (A[0,0]*A[2,2]**2 - A[0,2]**2 * A[2,2])))
            ) / self.__t

            self.__G12_eff = (
                A[2,2]
                - A[1,2] / A[1,1]
                + ((2*A[0,1]*A[0,2]*A[1,2]*A[1,1]
                    - (A[0,1]*A[1,2])**2
                    - (A[0,2]*A[1,1])**2)
                    / (A[0,0]*A[1,1]**2 - A[0,1]**2 * A[1,1]))
            ) / self.__t

        else:
            A = self.__A
            B = self.__B
            D = self.__D
            det_ABD = np.linalg.det(np.vstack((np.hstack((A, B)),
                                               np.hstack((B, D)))))

            # NASA, Eq. 84
            Astack = np.array([[A[1,1], A[1,2]], [A[1,2], A[2,2]]])
            Bstack = np.array([[B[0,1], B[1,1], B[1,2]],
                               [B[0,2], B[1,2], B[2,2]]])
            denom_mat = np.vstack((np.hstack((Astack, Bstack)),
                                   np.hstack((np.transpose(Bstack), D))))

            self.__E1_eff = (det_ABD / np.linalg.det(denom_mat)) / self.__t

            # NASA, Eq. 85
            Astack = np.array([[A[0,0], A[0,2]],
                               [A[0,0], A[2,2]]])
            Bstack = np.array([[B[0,0], B[0,1], B[0,2]],
                               [B[0,2], B[1,2], B[2,2]]])
            denom_mat = np.vstack((np.hstack((Astack, Bstack)),
                                   np.hstack((np.transpose(Bstack), D))))

            self.__E2_eff = (det_ABD / np.linalg.det(denom_mat)) / self.__t

            # NASA, Eq. 86
            Astack = np.array([[A[0,0], A[0,1]],
                               [A[0,1], A[1,1]]])
            Bstack = np.array([[B[0,0], B[0,1], B[0,2]],
                               [B[0,1], B[1,1], B[1,2]]])
            denom_mat = np.vstack((np.hstack((Astack, Bstack)),
                                   np.hstack((np.transpose(Bstack), D))))

            self.__G12_eff = (det_ABD / np.linalg.det(denom_mat)) / self.__t

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

        if type(ply) != Lamina:
            raise TypeError("Laminates may only contain Lamina objects")
        else:
            return True

    def clear(self):
        """Clear laminate of all plies."""

        self.__A = np.zeros((3, 3))
        self.__B = np.zeros((3, 3))
        self.__D = np.zeros((3, 3))
        self.__t = 0
        self.__layup = []
        self.__dT = 0
        self.__dM = 0
        self.__plies.clear()
        self.__N_m = np.zeros((3, 1))
        self.__M_m = np.zeros((3, 1))
        self.__N_t = np.zeros((3, 1))
        self.__M_t = np.zeros((3, 1))
        self.__N_h = np.zeros((3, 1))
        self.__M_h = np.zeros((3, 1))
        self.__E1_eff = 0
        self.__E2_eff = 0
        self.__G12_eff = 0

    def from_csv(self, inputFile, append=False):
        """Determine laminate properties from input CSV file.

        CSV must have headers:
            't' = thickness
            'theta' = angle of ply relative to laminate
            'E1' = modulus in lamina 1-direction
            'E2' = modulus in lamina 2-direction
            'nu12' = Poisson's ratio
            'G12' = shear modulus
            'a11' = coeff. of thermal expansion in 1-direction
            'a22' = coeff. of thermal expansion in 2-direction
            'b11' = coeff. of hygral expansion in 1-direction
            'b22' = coeff. of hygral expansion in 1-direction
        """

        if append == False:
            self.clear()

        for ply in pd.read_csv(inputFile).to_dict('records'):
            self.__plies.append(Lamina(t=ply['t'],
                                       theta=ply['theta'],
                                       E1=ply['E1'],
                                       E2=ply['E2'],
                                       nu12=ply['nu12'],
                                       G12=ply['G12'],
                                       a11=ply['a11'],
                                       a22=ply['a22'],
                                       b11=ply['b11'],
                                       b22=ply['b22']))

        self.__update()

    def pop(self, key):
        """Pop lamina out of Laminate."""

        self.__plies.pop(key)
        self.__update()

    def sort(self):
        """Sort the laminate plies by ID."""

        self.__plies.sort(key=lambda ply: ply.id)
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
                each.id = _
            del _
        self.__update()

    def is_balanced(self):
        """Return if True if laminate is balanced."""

        pass

    def is_symmetric(self):
        """Return if True if laminate is symmetric."""

        if len(self) % 2 == 0:  # even number of plies
            rng = int(len(self)/2)
        else:  # odd number of plies
            rng = int((len(self) - len(self)%2)/2)

        for i in range(0, rng):
            if self[i] != self[-(i+1)]:
                return False

        return True

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
        """numpy.array: mechanical running load applied to the laminate."""
        return self.__N_m

    @N_m.setter
    def N_m(self, newLoad):
        self.__N_m = newLoad
        self.__update()

    @property
    def M_m(self):
        """numpy.array: mechanical running moment applied to the laminate."""
        return self.__M_m

    @M_m.setter
    def M_m(self, newLoad):
        self.__M_m = newLoad
        self.__update()

    @property
    def N_t(self):
        """3x1 numpy.array: thermally induced running load."""
        return self.__N_t

    @property
    def M_t(self):
        """3x1 numpy.array: thermally induced running moment."""
        return self.__M_t

    @property
    def N_h(self):
        """3x1 numpy.array: hygroscopically induced running load."""
        return self.__N_h

    @property
    def M_h(self):
        """3x1 numpy.array: hygroscopically induced running moment."""
        return self.__M_h

    @property
    def dT(self):
        """float or int: change in temperature."""
        return self.__dT

    @dT.setter
    def dT(self, new_dT):
        self.__dT = new_dT
        self.__update()

    @property
    def dM(self):
        """float or int: change in moisture."""
        return self.__dM

    @dM.setter
    def dM(self, new_dM):
        self.__dM = new_dM
        self.__update()

    @property
    def E1_eff(self):
        """float: effective extensional modulus in the laminate 1-direction."""
        return self.__E1_eff

    @property
    def E2_eff(self):
        """float: ffective extensional modulus in the laminate 2-direction."""
        return self.__E2_eff

    @property
    def G12_eff(self):
        """float: effective in-plane shear modulus."""
        return self.__G12_eff

    @property
    def e_0(self):
        """3x1 numpy.array: total midplane strain on the laminate"""
        return self.__e_0

    @property
    def e_0m(self):
        """3x1 numpy.array: midplane strain due to mechanical loading"""
        return self.__e_0m

    @property
    def k_0m(self):
        """3x1 numpy.array: midplane curvature due to mechanical loading"""
        return self.__k_0m

    @property
    def e_0t(self):
        """3x1 numpy.array: midplane strain due to thermal loading"""
        return self.__e_0t

    @property
    def k_0t(self):
        """3x1 numpy.array: midplane curvature due to thermal loading"""
        return self.__k_0t

    @property
    def e_0h(self):
        """3x1 numpy.array: midplane strain due to hygral loading"""
        return self.__e_0h

    @property
    def k_0h(self):
        """3x1 numpy.array: midplane curvature due to thermal loading"""
        return self.__k_0h
