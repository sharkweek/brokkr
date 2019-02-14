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

from .lamina_subclass import Lamina, Ply
import numpy as np
from numpy import array, hstack, vsplit, vstack, zeros
from numpy.linalg import inv, det
import pandas as pd


class Laminate(dict):
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
    __protected__ = [
        'N_t', 'N_h', 'M_t', 'M_h' 'E1_eff', 'E2_eff', 'G12_eff', 'e_0', 'k_0'
    ]
    __unprotected__ = ['__locked']
    __updated__ = ['dT', 'dM', 'N_m', 'M_m']
    __slots__ = __protected__ + __unprotected__ + __updated__

    def __unlock(func):
        """Decorator function to allow for free attribute assignment."""

        def wrapper(self, *args, **kwargs):
            super().__setattr__('_Ply__locked', False)
            func(self, *args, **kwargs)
            super().__setattr__('_Ply__locked', True)

        return wrapper

    @__unlock
    def __init__(self, *plies):
        """Initialize with a list of Lamina objects.

        Accepts a series of Lamina objects or a single string calling out the
        filepath to a csv file containg all the laminate data.
        """

        # create and zero out all properties needed for .__update()
        self.clear()

        # init from CSV
        if len(plies) == 1 and type(plies[0]) == str:
            self.from_csv(plies[0])

        # standard init
        else:
            for each in plies:
                self.check_lamina(each)
            # convert plies to dict with ply ids as keys
            plies = {i+1:j for i, j in enumerate(plies)}
            super().__init__(plies)
            self.__update()

    def __setattr__(self, attr, val):
        """Set attribute."""

        if self.__locked:
            # udpate laminate after updated properties are set
            if attr in self.__updated__:
                super().__setattr__(attr, val)
                self.__update()

            # don't set protected values
            elif attr in self.__protected__:
                raise AttributeError(self.__name__ + ".%s" % attr
                                     + "is a derived value and cannot be set")

            # check if attribute is an unprotected value
            elif attr in self.__unprotected__:
                super().__setattr__(attr, val)

        else:
            super().__setattr__(attr, val)

    def __delitem__(self, key):
        """Delete a ply."""
        super().__delitem__(key)
        self.__update()

    def __setitem__(self, key, new_ply):
        """Set a specific ply to a new Lamina at the requested index."""
        if key in self:
            self.check_lamina(new_ply)
            super().__setitem__(key, new_ply)
            self.__update()
        else:
            raise KeyError(
                "Ply %s does not exist in laminate. To add a ply use 'append'"
                % key
            )

    def __repr__(self):
        """Set the representation of the laminate."""
        return "Laminate layup: " + self.layup.__repr__()

    @__unlock
    def __update(self):
        """Update the ply and laminate attributes based on laminate stackup."""

        for ply in self:
            self.check_lamina(ply)
            ply.laminate = self

        # reset internal values
        self.t = sum([ply.t for ply in self])
        self.N_t = zeros((3, 1))   # thermal running loads
        self.M_t = zeros((3, 1))   # thermal running moments
        self.N_h = zeros((3, 1))   # hygroscopic running loads
        self.M_h = zeros((3, 1))   # hygroscopic running moments

        # calculate ply bottom planes
        zk1 = self.t / 2
        for i in [x for x in sorted(self, reverse=True)]:  # count down the top
            zk1 -= self[i].t
            self[i].zk1 = zk1

        # calculate the A, B, and D matricies; thermal and hygral loads
        self.A = zeros((3, 3))
        self.B = zeros((3, 3))
        self.D = zeros((3, 3))

        for ply in self:
            # NASA-RP-1351, Eq (45) through (47)
            self.A += ply.Qbar * (ply.t)
            self.B += (1/2) * ply.Qbar * (ply.t * ply.z)
            self.D += (1/3) * ply.Qbar * (ply.t**3 / 12 + ply.t * ply.z**2)

            # thermal running loads
            # NASA-RP-1351, Eq (92)
            self.N_t += (ply.Qbar @ ply.e_tbar) * (ply.zk - ply.zk1)
            self.M_t += (ply.Qbar @ ply.e_tbar) * (ply.zk**2 - ply.zk1**2)/2

            # hygroscopic running loads
            # NASA-RP-1351, Eq (93)
            self.N_h += (ply.Qbar @ ply.e_hbar) * (ply.zk - ply.zk1)
            self.M_h += (ply.Qbar @ ply.e_hbar) * (ply.zk**2 - ply.zk1**2)/2

        # intermediate star matrices
        # NASA-RP-1351 Eq (50)-(52a)
        A_star = inv(self.A)
        B_star = - (A_star @ self.B)
        C_star = self.B @ A_star
        D_star = self.D - (self.B @ A_star @ self.B)

        D_prime = inv(D_star)
        C_prime = - (D_prime @ C_star)
        B_prime = B_star @ D_prime
        A_prime = A_star - (B_star @ D_prime @ C_star)

        ABD_prime = vstack((hstack((A_prime, B_prime)),
                            hstack((C_prime, D_prime))))

        # the midplane strains
        # Jones, (4.108) and (4.109). Also see Section 4.5.4
        mech_load = vstack((self.N_m, self.M_m))
        thrm_load = vstack((self.N_t, self.M_t))
        hygr_load = vstack((self.N_h, self.M_h))
        self.e_0m, self.k_0m = vsplit((ABD_prime @ mech_load), 2)
        self.e_0t, self.k_0t = vsplit((ABD_prime @ thrm_load), 2)
        self.e_0h, self.k_0h = vsplit((ABD_prime @ hygr_load), 2)

        # add ply strains due to mechanical loading (in laminate orientation)
        # Jones, Eq (4.13)
        self.e_0 = self.e_0m + self.e_0t + self.e_0h
        self.k_0 = self.k_0m + self.k_0t + self.k_0h

        for ply in self:
            ply.e_mbar = ply.T @ (self.e_0 + ply.z * self.k_0)

        # calculate effective laminate properties
        # NASA, Section 5
        # NOTE: This is only valuable for symmetric laminates.
        if self.is_symmetric():
            # NASA, Eq. 64, 70, and 76
            A = self.A
            self.E1_eff = (
                A[0,0]
                + A[0,1] * ((A[1,2]*A[0,2] - A[0,1]*A[2,2])
                            / (A[1,1]*A[2,2] - A[1,2]**2))
                + A[0,2] * (-A[0,2]/A[2,2]
                            + ((A[1,2]*A[0,1]*A[2,2] - A[1,2]**2 * A[0,2])
                                / (A[1,1]*A[2,2]**2 - A[1,2]**2 * A[2,2])))
            ) / self.t

            self.E2_eff = (
                A[0,1] * ((A[0,2]*A[1,2] - A[0,1]*A[2,2])
                        / (A[0,0]*A[2,2] - A[0,2]**2))
                + A[1,1]
                + A[1,2] * (-A[1,2]/A[2,2]
                            + ((A[0,2]*A[0,1]*A[2,2] - A[0,2]**2 * A[1,2])
                                / (A[0,0]*A[2,2]**2 - A[0,2]**2 * A[2,2])))
            ) / self.t

            self.G12_eff = (
                A[2,2]
                - A[1,2] / A[1,1]
                + ((2*A[0,1]*A[0,2]*A[1,2]*A[1,1]
                    - (A[0,1]*A[1,2])**2
                    - (A[0,2]*A[1,1])**2)
                    / (A[0,0]*A[1,1]**2 - A[0,1]**2 * A[1,1]))
            ) / self.t

        else:
            A = self.A
            B = self.B
            D = self.D
            det_ABD = det(vstack((hstack((A, B)),
                                  hstack((B, D)))))

            # NASA, Eq. 84
            Astack = array([[A[1,1], A[1,2]], [A[1,2], A[2,2]]])
            Bstack = array([[B[0,1], B[1,1], B[1,2]],
                            [B[0,2], B[1,2], B[2,2]]])
            denom_mat = vstack((hstack((Astack, Bstack)),
                                hstack((Bstack.T, D))))

            self.E1_eff = (det_ABD / det(denom_mat)) / self.t

            # NASA, Eq. 85
            Astack = array([[A[0,0], A[0,2]],
                            [A[0,0], A[2,2]]])
            Bstack = array([[B[0,0], B[0,1], B[0,2]],
                            [B[0,2], B[1,2], B[2,2]]])
            denom_mat = vstack((hstack((Astack, Bstack)),
                                hstack((Bstack.T, D))))

            self.E2_eff = (det_ABD / det(denom_mat)) / self.t

            # NASA, Eq. 86
            Astack = array([[A[0,0], A[0,1]],
                            [A[0,1], A[1,1]]])
            Bstack = array([[B[0,0], B[0,1], B[0,2]],
                            [B[0,1], B[1,1], B[1,2]]])
            denom_mat = vstack((hstack((Astack, Bstack)),
                                hstack((Bstack.T, D))))

            self.G12_eff = (det_ABD / det(denom_mat)) / self.t

    def append(self, newPly):
        """Add a new ply to the Laminate.

        Added plies are assumed to be placed on TOP SURFACE of existing
        laminate.
        """

        self.check_lamina(newPly)
        super().__setitem__(len(self)+1, newPly)
        self.__update()

    def check_lamina(self, ply):
        """Check if value is a ply and raises error if not."""

        if type(ply) != Lamina:
            raise TypeError("Laminates may only contain Lamina objects")
        else:
            return True

    def clear(self):
        """Clear laminate of all plies."""

        self.__protect(False)
        self.A = zeros((3, 3))
        self.B = zeros((3, 3))
        self.D = zeros((3, 3))
        self.t = 0
        self.dT = 0
        self.dM = 0
        self.N_m = zeros((3, 1))
        self.M_m = zeros((3, 1))
        self.N_t = zeros((3, 1))
        self.M_t = zeros((3, 1))
        self.N_h = zeros((3, 1))
        self.M_h = zeros((3, 1))
        self.E1_eff = 0
        self.E2_eff = 0
        self.G12_eff = 0
        super().__clear__()
        self.__protect(True)

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

        start = len(self)
        for i, rec in enumerate(pd.read_csv(inputFile).to_dict('records')):
            ply = Ply(t=rec['t'],
                      theta=rec['theta'],
                      E1=rec['E1'],
                      E2=rec['E2'],
                      nu12=rec['nu12'],
                      G12=rec['G12'],
                      a11=rec['a11'],
                      a22=rec['a22'],
                      b11=rec['b11'],
                      b22=rec['b22'])
            super().__setitem__(start + i + 1, ply)

        self.__update()

    def pop(self, key):
        """Pop lamina out of Laminate."""

        raise AttributeError(self.__name__ + " has no attribute 'pop'")

    def sort(self):
        """Sort the laminate plies by ID."""

        self.__plies.sort(key=lambda ply: ply.id)
        self.__update()

    def flip(self, renumber=False):
        """Flip the stacking sequence of the Lamina."""

        self.clear()
        vals = [self[x] for x in self]  # ply object references
        rev = sorted(self, reverse=True)  # reversed ply ids
        for i in rev:
            i-=1  # account for list indexing 0
            self[rev[i]] = vals[i]

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
        return [ply.theta for ply in self.__plies]

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
