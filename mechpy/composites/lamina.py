"""Defines the Lamina class for use with mechpy.

Assumptions
-----------
* Lamina are generally orthotropic.
* Lamina properties and values are assumed to be in the lamina coordinate
  system
  (i.e. 12-plane).
* Theta is in degrees, measured from the laminate x-axis to the lamina 1-axis.
* The laminate z-direction is upward, away from the bottom surface.
* Reported strain values are in engineering strain.
* All loads and moments are running values supplied as force or moment per unit
  width.
* Unit systems are consistent (i.e. SI, US, or Imperial).
* Equations and symbol conventions are per NASA-RP-1351 'Basic Mechanics of
  Laminated Composite Plates' and Jones' 'Mechanics Of Composite Materials'.

TODO:
-----
* [ ] add functionality to register laminate property to forwad updates to
      Laminate classes
* [ ] add error checking for attribute types
* [ ] add a property for the compliance matrix Sk
* [ ] add failure indeces
* [ ] remove laminate orientation calculations and have them be returned as
      calculated properties
"""

import numpy as np
from numpy import array, zeros
from numpy.linalg import inv
from math import isnan, cos, sin, radians


class Lamina:
    """
    Lamina(t, E1, E2, nu12, G12, a11, a22, b11, b22)

    An individual lamina.

    Attributes
    ----------
    t : float
        lamina thickness
    E1, E2, G12 : float
        elastic moduli in the lamina 1-, 2-, and 12-directions
    nu12 : float
        Poisson's ratio in the 12-plane
    a11, a22 : float
        coefficients of thermal expansion (CTE) in the lamina 1- and 2-
        directions
    b11, b22 : float
        coefficients of hygral expansion (CTE) in the lamina 1- and 3-directions

    Assumptions
    -----------
    * Each lamina is assumed to be orthotropic.

    """

    __name__ = 'Lamina'
    __slots__ = ['t', 'E1', 'E2', 'nu12', 'G12', 'a11', 'a22', 'b11', 'b22']

    def __init__(self,
                 t,     # lamina thickness
                 E1,    # Young's modulus in 1-direction
                 E2,    # Young's modulus in 2-direction
                 nu12,  # Poisson's ratio in 12-plane
                 G12,   # shear modulus in 12-plane
                 a11,   # coeff. of thermal expansion in 1-direction
                 a22,   # coeff. of thermal expansion in 2-direction
                 b11,   # coeff. of moisture expansion in 1-direction
                 b22):  # coeff. of moisture expansion in 2-direction

        self.t = t
        self.E1 = E1
        self.E2 = E2
        self.nu12 = nu12
        self.G12 = G12
        self.a11 = a11
        self.a22 = a22
        self.b11 = b11
        self.b22 = b22

    @property
    def is_fully_defined(self):
        """Check if lamina is fully defined with proper attr data types."""

        if self.t == 0 or isnan(self.t):
            raise TypeError("lamina.tk must be a non-zero number")
        elif self.E1 == 0 or isnan(self.E1):
            raise TypeError("lamina.E1 must be a non-zero number")
        elif self.E2 == 0 or isnan(self.E2):
            raise TypeError("lamina.E2 must be a non-zero number")
        elif (self.nu12 >= 1) or isnan(self.nu12):
            raise TypeError("""lamina.nu12 must be less than or equal to 1""")
        elif self.G12 == 0 or isnan(self.G12):
            raise TypeError("lamina.G12 must be a non-zero number")
        else:
            return True


class Ply(Lamina):
    """A Ply for use in a Laminate."""

    __name__ = 'Ply'
    __baseattr__ = Lamina.__slots__ + ['theta']
    __calcattr__ = ['Q', 'Qbar', 'T', 'Tinv', 'e_t', 'e_h', 'z', 'e_m',
                    'laminate']
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
            super().__setattr__('_Ply__locked', False)
            func(self, *args, **kwargs)
            super().__setattr__('_Ply__locked', True)
        return wrapper

    @__unlock
    def __init__(self,
                 laminate,
                 t,
                 theta,
                 E1,
                 E2,
                 nu12,
                 G12,
                 a11,
                 a22,
                 b11,
                 b22):

        self.laminate = laminate
        self.z = 0
        self.theta = theta
        self.e_m = zeros((3, 1))
        self.e_t = zeros((3, 1))
        self.e_h = zeros((3, 1))

        super().__init__(t, E1, E2, nu12, G12, a11, a22, b22, b22)
        self.__update()

    def __setattr__(self, attr, val):
        if self.__locked:
            # udpate laminate after updated properties are set
            if attr in self.__baseattr__:
                super().__setattr__(attr, val)
                self.__update()

            # don't set protected values
            elif attr in self.__calcattr__:
                raise AttributeError(self.__name__ + ".%s" % attr
                                     + " is a derived value and cannot be set")

            # update the laminate
            if self.laminate:
                self.laminate._Laminate__update()

        else:
            super().__setattr__(attr, val)

    @__unlock
    def __update(self):
        """Update calculated attributes."""

        # on-axis reduced stiffness matrix, Q
        # NASA-RP-1351, Eq (15)
        nu21 = self.nu12 * self.E2 / self.E1  # Jones, Eq (2.67)
        q11 = self.E1 / (1 - self.nu12 * nu21)
        q12 = self.nu12 * self.E2 / (1 - self.nu12 * nu21)
        q22 = self.E2 / (1 - self.nu12 * nu21)
        q66 = self.G12

        self.Q = array([[q11, q12, 0], [q12, q22, 0], [0, 0, q66]])

        # the transformation matrix and its inverse
        # create intermediate trig terms
        m = cos(radians(self.theta))
        n = sin(radians(self.theta))

        # create transformation matrix and inverse
        self.T = array([[m**2, n**2, 2 * m * n],
                        [n**2, m**2, -2 * m * n],
                        [-m * n, m * n, m**2 - n**2]])
        self.Tinv = inv(self.T)

        # the transformed reduced stiffness matrix (laminate coordinate system)
        # Jones, Eq (2.84)
        self.Qbar = self.Tinv @ self.Q @ self.Tinv.T

        # thermal and hygral strains in laminate and lamina c-systems
        # NASA-RP-1351 Eq (90), (91), and (95)
        self.e_t = array([[self.a11], [self.a22], [0]]) * self.laminate.dT
        self.e_h = array([[self.b11], [self.b22], [0]]) * self.laminate.dM

    @staticmethod
    def new_from_lamina(lamina,
                        laminate,
                        t,
                        z,
                        theta):
        """Create a new Ply object from a Lamina object."""

        return Ply(laminate=laminate,
                   t=t,
                   theta=theta,
                   E1=lamina.E1,
                   E2=lamina.E2,
                   nu12=lamina.nu12,
                   G12=lamina.G12,
                   a11=lamina.a11,
                   a22=lamina.a22,
                   b11=lamina.b11,
                   b22=lamina.b22)

    @property
    def zk(self):
        """The vertical location of the lamina's top plane.

        zk is measured from the midplane of the laminate with the convention
        of positive pointing upward.
        """
        return self.z + self.t / 2

    @property
    def zk1(self):
        """The vertical location of the lamina's bottom plane.

        zk1 is measured from the midplane of the laminate with the convention
        of positive pointing upward. Note: zk1 is z_(k-1).
        """
        return self.z - self.t / 2

    @zk1.setter
    def zk1(self, new_zk1):
        self.z = new_zk1 + self.t / 2

    @property
    def e_tbar(self):
        """Transformed thermal strain."""
        return self.Tinv @ self.e_t  # NASA-RP-1351 Eq (90)

    @property
    def e_hbar(self):
        """Transformed hygral strain."""
        return self.Tinv @ self.e_h  # NASA-RP-1351 Eq (95)

    @property
    def e(self):
        """Total strain."""
        return self.e_m + self.e_t + self.e_h
