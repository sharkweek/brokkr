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
from numpy import array, hstack, vsplit, vstack, zeros
from numpy.linalg import inv, det
import pandas as pd


class Laminate(dict):
    """
    Laminate(*plies)

    A composite laminate made up of multiple orthotropic lamina objects.

    Parameters
    ----------
    *plies : Lamina, Ply, or string
        lamina or plies that make up the laminate ordered from the bottom
        surface to the top surface. In creating the Laminate, either a
        sequence of objects or a file path to a CSV file containing all the
        appropriate data can be supplied.

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

    # lists of attributes for method filtering
    __name__ = 'Laminate'
    __baseattr__ = ['dT', 'dM', 'N_m', 'M_m']
    __calcattr__ = ['N_t', 'N_h', 'M_t', 'M_h', 'E1_eff', 'E2_eff', 'G12_eff',
                    'e_0m', 'e_0t', 'e_0h', 'k_0m', 'k_0t', 'k_0h', 'A', 'B',
                    'D', 't']
    __slots__ = __baseattr__ + __calcattr__ + ['__locked']

    def __unlock(func):
        """Decorate methods to unlock attributes.

        Parameters
        ----------
        update : bool
            Determines whether to run the __update() method after execution.

        Returns
        -------
        function
            An unprotected function.

        """

        def wrapper(self, *args, **kwargs):
            super().__setattr__('_Laminate__locked', False)
            func(self, *args, **kwargs)
            super().__setattr__('_Laminate__locked', True)
        return wrapper

    @__unlock
    def __init__(self, *plies):
        """Initialize with a list of Lamina or Ply objects.

        Accepts a series of Lamina or Ply objects or a single string calling
        out the filepath to a csv file containg all the laminate data.
        """

        # create and zero out all properties needed for .__update()
        self.clear()

        # init from CSV
        if len(plies) == 1 and type(plies[0]) == str:
            self.from_csv(plies[0])

        # standard init
        else:
            for ply in plies:
                # convert to Ply objects
                ply = self.check_ply(ply)
            # convert plies to dict with ply ids as keys
            plies = {i+1: j for i, j in enumerate(plies)}
            super().__init__(plies)
            self.__update()

    def __setattr__(self, attr, val):
        """Set attribute."""

        if self.__locked:
            # udpate laminate after updated properties are set
            if attr in self.__baseattr__:
                super().__setattr__(attr, val)
                self.__update()

            # don't set protected values
            elif attr in self.__calcattr__:
                raise AttributeError(self.__name__ + ".%s" % attr
                                     + "is a derived value and cannot be set")

        else:
            super().__setattr__(attr, val)

    def __delitem__(self, key):
        """Delete a ply."""
        super().__delitem__(key)
        self.__update()

    def __setitem__(self, key, new_ply):
        """Set a specific ply to a new Lamina at the requested index."""
        if key in self:
            new_ply = self.check_ply(new_ply)
            super().__setitem__(key, new_ply)
            self.__update()
        else:
            raise KeyError(
                "Ply %s does not exist in laminate. Use 'append' to add a ply"
                % key
            )

    def __repr__(self):
        """Set the representation of the laminate."""
        return "Laminate layup: " + self.layup.__repr__()

    @__unlock
    def __update(self):
        """Update the ply and laminate attributes based on laminate stackup."""

        # reset internal values
        self.t = sum([self[i].t for i in self])
        self.N_t = zeros((3, 1))   # thermal running loads
        self.M_t = zeros((3, 1))   # thermal running moments
        self.N_h = zeros((3, 1))   # hygroscopic running loads
        self.M_h = zeros((3, 1))   # hygroscopic running moments

        # calculate ply bottom planes
        zk1 = self.t / 2
        for i in [x for x in sorted(self, reverse=True)]:  # count down the top
            zk1 -= self[i].t
            super(Ply, self[i]).__setattr__('z', zk1 + self[i].t/2)

        # calculate the A, B, and D matricies; thermal and hygral loads
        self.A = zeros((3, 3))
        self.B = zeros((3, 3))
        self.D = zeros((3, 3))

        for i in self:
            ply = self[i]
            # NASA-RP-1351, Eq (45) through (47)
            self.A += ply.Qbar * (ply.t)
            self.B += (1/2) * ply.Qbar * (ply.zk**2 - ply.zk1**2)
            self.D += (1/3) * ply.Qbar * (ply.zk**3 - ply.zk1**3)

            # thermal running loads
            # NASA-RP-1351, Eq (92)
            ply._Ply__update()  # force ply to update
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
        for ply in self:
            e_m = self[ply].T @ (self.e_0m + self[ply].z * self.k_0m)
            super(Ply, self[ply]).__setattr__('e_m', e_m)

        # calculate effective laminate properties
        # NASA, Section 5
        def minor(arr, index):
            """Calculate the minor of a matrix."""
            return np.delete(
                np.delete(arr, index[0], axis=0), index[1], axis=1
            )

        # calculate effective moduli
        # NASA, Eq. 84
        numer = det(self.ABD)
        denom = det(minor(self.ABD, (0, 0)))
        self.E1_eff = (numer / denom) / self.t

        # NASA, Eq. 85
        denom = det(minor(self.ABD, (1, 1)))
        self.E2_eff = (numer / denom) / self.t

        # NASA, Eq. 86
        denom = det(minor(self.ABD, (2, 2)))
        self.G12_eff = (numer / denom) / self.t

    def append(self, newPly):
        """Add a new ply to the Laminate.

        Added plies are assumed to be placed on TOP SURFACE of existing
        laminate.
        """

        newPly = self.check_ply(newPly)
        super().__setitem__(len(self)+1, newPly)
        self.__update()

    def check_ply(self, ply):
        """Check if value is a ply and raises error if not."""

        if type(ply) == Ply:
            pass
        elif ply.laminate != self:
            super(Ply, ply).__setattr__('laminate', self)
        elif type(ply) == Lamina:
            ply = Ply.new_from_lamina(ply, laminate=self)
        else:
            raise TypeError("Laminates may only contain Ply or Lamina objects")

        return ply

    @__unlock
    def clear(self):
        """Clear laminate of all plies."""

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
        super().clear()

    def from_csv(self, inputFile, append=False):
        """Determine laminate properties from input CSV file.

        Parameters
        ----------
        inputFile : string
            file path to CSV file; the CSV must have headers:
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
        append : bool
            Flag indicating whether to clear the laminate before adding plies
            from the input file.

        Returns
        -------
        None

        Notes
        -----
        The CSV is assumed to have the plies listed in ascending order (i.e.
        the bottom ply is listed first).

        """

        if append is False:
            self.clear()

        start = len(self)
        for i, rec in enumerate(pd.read_csv(inputFile).to_dict('records')):
            ply = Ply(laminate=self,
                      t=rec['t'],
                      theta=rec['theta'],
                      E1=rec['E1'],
                      E2=rec['E2'],
                      nu12=rec['nu12'],
                      G12=rec['G12'],
                      a11=rec['a11'],
                      a22=rec['a22'],
                      b11=rec['b11'],
                      b22=rec['b22'])
            # user the super() method to avoid updating laminate for every ply
            super().__setitem__(start + i + 1, ply)

        self.__update()

    def pop(self, key, renumber=True):
        """Remove item from Laminate and return Ply."""
        p = super().pop(key)
        if renumber:
            ids = [i for i in sorted(self)]
            plies = [self[i] for i in sorted(self)]
            self.clear()

            for i in range(len(plies)):
                super().__setitem__(i+1, plies[i])

        self.__update()
        return p


    def flip(self):
        """Flip the stacking sequence of the Lamina."""
        ids = [i for i in sorted(self)]
        plies = [self[i] for i in sorted(self, reverse=True)]
        for i in range(len(plies)):
            self[ids[i]] = plies[i]

        self.__update()

    @property
    def is_balanced(self):
        """Return if True if laminate is balanced."""

        pass

    @property
    def is_symmetric(self):
        """Return if True if laminate is symmetric."""

        if len(self) % 2 == 0:  # even number of plies
            rng = int(len(self) / 2)
        else:  # odd number of plies
            rng = int((len(self) - len(self) % 2) / 2)

        for i in range(0, rng):
            if self[i] != self[-(i+1)]:
                return False

        return True

    @property
    def ABD(self):
        """ABD matrix."""
        return vstack((hstack((self.A, self.B)),
                       hstack((self.B, self.D))))

    @property
    def layup(self):
        """Laminate layup."""
        return [self[ply].theta for ply in self]

    @property
    def e_0(self):
        """Total midplane strain."""
        return self.e_0m + self.e_0t + self.e_0h

    @property
    def k_0(self):
        """Total laminate curvature."""
        return self.k_0m + self.k_0t + self.k_0h
