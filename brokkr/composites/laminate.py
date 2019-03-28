"""Defines the Laminate class for use with brokkr.

Notes
-----
See `brokkr.composites` documentation for relevant assumptions.

TODO:
-----
* [ ] add unyt units
* [ ] Create is_balanced() method
* [ ] create functions to import Lamina and Laminate properties from Nastran
  BDF
* [ ] fix calculations for effective laminate properties
* [ ] reorganize properties and setters to be adjacent
"""

from .lamina import Lamina, Ply
from brokkr._exceptions import (
    UnitDimensionError,
    CalculatedAttributeError,
)
from brokkr.config import USYS, UREG
from numpy import hstack, vsplit, vstack, zeros
from numpy.linalg import inv, det
import pandas as pd
from brokkr.mech_math import matrix_minor
from unyt import UnitSystem, unyt_array
from unyt.dimensions import (
    length,
    dimensionless,
    temperature,
    force,
    pressure
)

__all__ = ['Laminate']


class Laminate(dict):
    """A composite laminate made up of multiple orthotropic lamina objects.

    ``Laminate`` instances contain ``Ply`` instances and used for calculating
    the effects on individual lamina in a laminate stackup. All calculations
    adhere to classical laminate theory (CLT).

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
    dT : int or float
        change in temperature
    dM : int or float
        change in moisture
    N_m, N_t, N_h : 3x1 numpy.array_like
        mechanical, thermal, and hygroscopic running loads
    M_m, M_t, M_h : 3x1 numpy.array_like
        mechanical, thermal, and hygroscopic running moments
    Ex, Ey, Gxy : float
        effective laminate moduli
    e_0m, e_0t, e_0h : 3x1 numpy.array_like
        laminate midplane strain due to mechanical, thermal, and hygroscopic
        loading
    k_0m, k_0t, k_0h : numpy.array_like
        laminate midplane curvature due to mechanical, thermal, and hygroscopic
        loading

    """

    # lists of attributes for method filtering
    _param_dims = {
        'dT': (temperature,),
        'dM': (dimensionless,),
        'N_m': (force / length,),
        'M_m': (force * length / length,),
        'N_t': (force / length,),
        'N_h': (force / length,),
        'M_t': (force * length / length,),
        'M_h': (force * length / length,),
        'Ex': (pressure,),
        'Ey': (pressure,),
        'Gxy': (pressure,),
        'e_0m': (dimensionless,),
        'e_0t': (dimensionless,),
        'e_0h': (dimensionless,),
        'k_0m': (dimensionless,),
        'k_0t': (dimensionless,),
        'k_0h': (dimensionless,),
        'A': (pressure * length,),
        'B': (pressure * length**2,),
        'D': (pressure * length**3,),
        't': (length,)
    }

    _base_attr = tuple(_param_dims) + ('usys',)
    _calc_attr = ('N_t', 'N_h', 'M_t', 'M_h', 'Ex', 'Ey', 'Gxy', 'e_0m',
                  'e_0t', 'e_0h', 'k_0m', 'k_0t', 'k_0h', 'A', 'B', 'D', 't')
    __slots__ = _base_attr + _calc_attr + ('__locked',)

    def __unlock(locked_func):
        """Decorate methods to unlock attributes.

        Parameters
        ----------
        locked_func : bool
            The function to unlock.

        Returns
        -------
        function
            An unlocked function.

        """

        def unlocked_func(self, *args, **kwargs):
            super().__setattr__('_Laminate__locked', False)
            locked_func(self, *args, **kwargs)
            super().__setattr__('_Laminate__locked', True)
        return unlocked_func

    @__unlock
    def __init__(self, *plies, dT=0, dM=0, N_m=zeros((3, 1)),
                 M_m=zeros((3, 1)), usys=USYS):
        # create and zero out all properties needed for .__update()
        self.usys = usys
        self.dT = dT
        self.dM = dM
        self.N_m = N_m
        self.M_m = N_m

        # init from CSV
        if len(plies) == 1 and type(plies[0]) == str:
            self.from_csv(plies[0])

        # standard init
        else:
            # convert plies to dict with ply ids as keys
            ply_dict = {i+1: self.check_ply(j) for i, j in enumerate(plies)}
            super().__init__(ply_dict)
            self.__update()

    def __setattr__(self, name, attr):
        """Extend ``__setattr__`` to update instance when attributes change."""

        def check_dim(name, attr):
            """Validate dimensions."""
            if name in self._param_dims:
                correct_dim = self._param_dims.get(name)

                if hasattr(attr, 'units'):
                    if attr.units.dimensions not in correct_dim:
                        raise UnitDimensionError(
                            name,
                            ' OR '.join([str(i) for i in correct_dim])
                        )
                    else:  # force consistent units
                        attr.convert_to_base(self.usys)

                else:
                    attr *= self.usys[correct_dim[0]]  # first dim

            # validate unit system
            if name == 'usys':
                # validate unit system
                if type(attr) != UnitSystem:
                    raise TypeError(f"`{name}` must be a `UnitSystem`")
                elif attr['temperature'].base_offset != 0:
                    raise TypeError(
                        "Unit system must have an absolute temperature"
                        + " unit (e.g. R or K)"
                    )

            return attr

        if self.__locked:
            # don't set protected values
            if name in self._calc_attr:
                raise CalculatedAttributeError(name)

            # validate units
            if name in self._base_attr:
                attr = check_dim(name, attr)
                super().__setattr__(name, attr)
                self.__update()

            if name == 'usys':
                if attr != self.usys:
                    super().__setattr__(name, attr)
                    self.__update()

        else:
            attr = check_dim(name, attr)
            super().__setattr__(name, attr)

    def __delitem__(self, key):
        """Extend ``__delitem__`` to update instance when ply is removed."""
        super().__delitem__(key)
        self.__update()

    def __setitem__(self, key, new_ply):
        """Extend ``__setitem__`` to ensure added item is a ``Ply``."""
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
        """Return the string representation."""
        repr = 'Laminate:'
        for ply in self:
            repr += f'\n    Ply {ply}: {self[ply].theta.v}'
        return repr

    @__unlock
    def __update(self):
        """Update the ply and laminate attributes based on laminate stackup.

        Notes
        -----
        ``unyt_array``s do not yet support matrix multiplication. All matrix
        multiplication operations have to be performed on the base ``ndarray``s
        and then have the units applied afterward.

        """

        # reset internal values
        self.t = sum([self[i].t for i in self])
        self.N_t = zeros((3, 1))
        self.M_t = zeros((3, 1))
        self.N_h = zeros((3, 1))
        self.M_h = zeros((3, 1))

        # calculate ply bottom planes
        zk1 = self.t / 2
        for i in [x for x in sorted(self, reverse=True)]:  # count down the top
            zk1 -= self[i].t
            # bypass lock to avoid recursive updates to self and ply
            super(Ply, self[i]).__setattr__('z', zk1 + self[i].t/2)

        # calculate the A, B, and D matricies; thermal and hygroscopic loads
        self.A = zeros((3, 3))
        self.B = zeros((3, 3))
        self.D = zeros((3, 3))

        for i in self:
            ply = self[i]
            ply._Ply__update()  # force ply to update

            # NASA-RP-1351, Eq (45) through (47)
            self.A += ply.Qbar * ply.t
            self.B += (1/2) * ply.Qbar * (ply.zk**2 - ply.zk1**2)
            self.D += (1/3) * ply.Qbar * (ply.zk**3 - ply.zk1**3)

            # thermal running loads
            # NASA-RP-1351, Eq (92)
            e_tbar = ply.Tinv @ ply.e_t
            self.N_t += (ply.Qbar @ e_tbar) * (ply.zk - ply.zk1)
            self.M_t += (ply.Qbar @ e_tbar) * (ply.zk**2 - ply.zk1**2)/2

            # hygroscopic running loads
            # NASA-RP-1351, Eq (93)
            e_hbar = ply.Tinv @ ply.e_h
            self.N_h += (ply.Qbar @ e_hbar) * (ply.zk - ply.zk1)
            self.M_h += (ply.Qbar @ e_hbar) * (ply.zk**2 - ply.zk1**2)/2

        # intermediate star matrices
        # NASA-RP-1351 Eq (50)-(52a)
        A_star = inv(self.A.v)
        B_star = - (A_star @ self.B.v)
        C_star = self.B.v @ A_star
        D_star = self.D.v - (self.B.v @ A_star @ self.B.v)

        D_prime = inv(D_star)
        C_prime = - (D_prime @ C_star)
        B_prime = B_star @ D_prime
        A_prime = A_star - (B_star @ D_prime @ C_star)

        ABD_prime = vstack((hstack((A_prime, B_prime)),
                            hstack((C_prime, D_prime))))

        # the midplane strains
        # Jones, (4.108) and (4.109). Also see Section 4.5.4
        mech_load = vstack((self.N_m.v, self.M_m.v))
        thrm_load = vstack((self.N_t.v, self.M_t.v))
        hygr_load = vstack((self.N_h.v, self.M_h.v))
        self.e_0m, self.k_0m = vsplit((ABD_prime @ mech_load), 2)
        self.e_0t, self.k_0t = vsplit((ABD_prime @ thrm_load), 2)
        self.e_0h, self.k_0h = vsplit((ABD_prime @ hygr_load), 2)

        # assign mechanical ply strains (in ply orientation)
        # Jones, Eq (4.13)
        for ply in self:
            e_m = self[ply].T @ (self.e_0m + self[ply].z.v * self.k_0m)
            e_t = self[ply].T @ (self.e_0t + self[ply].z.v * self.k_0t)
            e_h = self[ply].T @ (self.e_0h + self[ply].z.v * self.k_0h)
            super(Ply, self[ply]).__setattr__('e_m', unyt_array(e_m))

        # calculate effective moduli
        # NASA, Eq. 84
        numer = det(self.ABD)
        denom = det(matrix_minor(self.ABD, (0, 0)))
        self.Ex = ((numer / denom) / self.t.v) * self.usys['pressure']

        # NASA, Eq. 85
        denom = det(matrix_minor(self.ABD, (1, 1)))
        self.Ey = ((numer / denom) / self.t.v) * self.usys['pressure']

        # NASA, Eq. 86
        denom = det(matrix_minor(self.ABD.v, (2, 2)))
        self.Gxy = ((numer / denom) / self.t.v) * self.usys['pressure']

    def append(self, new_ply, angle=0):
        """Add a new ply to the `Laminate`.

        Parameters
        ----------
        new_ply : `Lamina` or `Ply`
            the ply to be added to the laminate

        Notes
        -----
        Added plies are assumed to be placed on TOP SURFACE of existing
        laminate.

        """

        new_ply = self.check_ply(new_ply)
        new_ply.theta = angle
        super().__setitem__(len(self)+1, new_ply)
        self.__update()

    def check_ply(self, obj):
        """Check if object is a `Ply` or `Lamina` object.

        Parameters
        ----------
        obj : any
            the object to be type checked

        Returns
        -------
            Ply
                the original `Ply` or converted `Lamina` object

        Raises
        ------
            TypeError
                if `obj` is neither a `Ply` or `Lamina` object

        """

        if type(obj) == Ply:
            pass
        elif type(obj) == Lamina:
            obj = Ply.from_lamina(obj, laminate=self, theta=0)
        elif obj.laminate != self:
            super(Ply, obj).__setattr__('laminate', self)
        else:
            raise TypeError("Laminates may only contain Ply or Lamina objects")

        if obj.usys != self.usys:
            obj.usys = self.usys

        return obj

    @__unlock
    def clear(self):
        """Clear laminate of all plies."""

        super().clear()
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
        self.Ex = 0
        self.Ey = 0
        self.Gxy = 0

    def from_csv(self, inputFile, append=False):
        """Determine laminate properties from input CSV file.

        Parameters
        ----------
        inputFile : string
            file path to CSV file
        append : bool, optional
            flag indicating whether to clear the laminate before adding plies
            from the input file; default is ``False``

        Notes
        -----
        The CSV is assumed to have the plies listed in ascending order (i.e.
        the bottom ply is listed first) and must have the following headers:

        ========== ==============================================
        Header     Description
        ========== ==============================================
        ``t``      thickness
        ``theta``  angle of ply relative to laminate
        ``E1``     modulus in lamina 1-direction
        ``E2``     modulus in lamina 2-direction
        ``nu12``   Poisson's ratio
        ``G12``    shear modulus
        ``a11``    coeff. of thermal expansion in 1-direction
        ``a22``    coeff. of thermal expansion in 2-direction
        ``b11``    coeff. of hygroscopic expansion in 1-direction
        ``b22``    coeff. of hygroscopic expansion in 1-direction
        ========== ==============================================

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

    def pop(self, key, renumber=True, as_lamina=False):
        """Remove item from Laminate and return Ply.

        Parameters
        ----------
        key : int
            index of the ply to be removed and returned
        renumber : bool, optional
            If True, renumber remaining plies. (``default=True``)
        as_lamina : bool, optional
            If True, return the ply as a ``Lamina`` object, otherwise, return
            as a ``Ply`` object. (``default=False``)

        Returns
        -------
        Ply : Ply or Lamina
            Description of returned object.

        """
        p = super().pop(key)

        # renumber the remaining plies
        if renumber:
            plies = [self[i] for i in sorted(self)]
            self.clear()

            for i in range(len(plies)):
                super().__setitem__(i+1, plies[i])

        self.__update()

        if not as_lamina:
            return p
        else:
            return Lamina(p.t, p.E1, p.E2, p.nu12, p.G12, p.a11, p.a22, p.b11,
                          p.b22)

    def flip(self):
        """Flip the stacking sequence of the Lamina."""
        ids = [i for i in sorted(self)]
        plies = [self[i] for i in sorted(self, reverse=True)]
        for i in range(len(plies)):
            self[ids[i]] = plies[i]

        self.__update()

    def reorient(self, angle):
        """Reorient laminate

        The units of `angle` are assumed to be the same angle units as the
        laminate. If a sequence is provided, it is assumed to start with the
        first ply and iterate through to the last in order.

        Parameters
        ----------
        angle : float or iterable of floats
            the angle to rotate the laminate or the sequence of individual
            lamina orientations

        """

        if type(angle) == tuple or type(angle) == list:
            if len(angle) != len(self):
                raise ValueError("`angle` must have values for each lamina")
            else:
                for i, j in enumerate(angle):
                    self[i + 1].theta = angle
        else:
            for i in self:
                self[i].theta = self[i].theta.v + angle

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
        """The laminate ABD matrix.

        Notes
        -----
        Because A, B, and D matrices each have different units and
        ``unyt_array``s must have uniform units, the ABD matrix is returned as
        a dimensionless ``unyt_array``.
        """
        return unyt_array(vstack((hstack((self.A.v, self.B.v)),
                                  hstack((self.B.v, self.D.v)))))

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
