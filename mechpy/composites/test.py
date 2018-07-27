"""PlyVectors for stress and strain."""

import numpy as np
import math


class PlyStress:
    """A stress vector for a single lamina.

    Used to easily translate between ply and laminate coordinate systems.
    """

    def __init__(self, vector, angle, ply_csys=True, rad=False):
        """Initialize the instance.

        Default assumes that vector is in the ply c-sys and angle is in
        degrees.
        """

        # error check inputs for types
        if type(vector) == np.ndarray and vector.shape != (3, 1):
            raise TypeError("'vector' must be a 3x1 numpy array.")
        elif type(ply_csys) != bool:
            raise TypeError("'ply_csys' must be a boolean.")
        elif type(rad) != bool:
            raise TypeError("'rad' must be a boolean")
        try:
            float(angle)
        except TypeError:
            print("'angle' must be a number.")
            raise

        self.__vector = vector
        self.__angle = angle
        self.__ply_csys = ply_csys
        self.__rad = rad

        self._update(vector, angle, ply_csys, rad)

    def _update(self,
                vector=None,
                angle=None,
                ply_csys=True,  # default values are in ply c-sys
                rad=None):
        """Update PlyVector when values are changed."""

        # check if angle or units are changed
        if (vector != self.__vector or
            angle != self.__angle or
            rad != self.__rad):
            transform = True

        # check if angle has changed
        if angle is not None:
            self.__angle = angle
        if rad is not None:
            self.__rad = rad
        if vector is not None:
            self.__vector = vector

        # create trig terms for transformation matrix calcs
        if transform:
            if self.__rad:
                m = math.cos(self.__angle)
                n = math.sin(self.__angle)
            else:
                m = math.cos(math.radians(self.__angle))
                n = math.sin(math.radians(self.__angle))

            # create transformation matrix and inverse
            self.__T = np.array([[m**2, n**2, 2*m*n],
                                 [n**2, m**2, -2*m*n],
                                 [-m*n, m*n, m**2 - n**2]])
            self.__Tinv = np.linalg.inv(self.__T)

            # update ply and laminate vectors based on which is provided
            if ply_csys:
                self.__vector = np.matmul(self.__Tinv, self.__lam)
            else:
                self.__lam = np.matmul(self.__T, self.__ply)

    @property
    def angle(self):
        """Lamina orientation angle w.r.t. the laminate c-sys."""

        return self.__angle

    @angle.setter
    def angle(self, new_angle):
        """Lamina orientation angle w.r.t. the laminate c-sys."""

        # make sure angle is a number
        try:
            float(new_angle)
        except TypeError:
            print("'angle' must be a number.")
            raise

        self._update(angle=new_angle)

    @property
    def rad(self):
        """Units flag for orientation angle."""

        return self.__rad

    @rad.setter
    def rad(self, new_rad):
        """Units flag for orientation angle."""

        if type(new_rad) != bool:
            raise TypeError("'rad' must be a boolean")

        self._update(rad=new_rad)

    @property
    def ply(self):
        """Lamina vector w.r.t. the ply c-sys."""

        return self.__ply

    @ply.setter
    def ply(self, new_vector):
        """Lamina vector w.r.t. the ply c-sys."""

        if type(new_vector) == np.ndarray and new_vector.shape != (3, 1):
            raise AttributeError("'new_vector' must be a 3x1 numpy array.")

        self._update(vector=new_vector, ply_csys=True)

    @property
    def lam(self):
        """Lamina vector w.r.t. the laminate c-sys."""

        return self.__lam

    @lam.setter
    def lam(self, new_vector):
        """Lamina vector w.r.t. the laminate c-sys."""

        if type(new_vector) == np.ndarray and new_vector.shape != (3, 1):
            raise AttributeError("'new_vector' must be a 3x1 numpy array.")

        self._update(vector=new_vector, ply_csys=False)

    @property
    def Tmat(self):
        """The transformation matrix."""

        return self.__T

    @property
    def Tinv(self):
        """The inverse of the transformation matrix."""

        return self.__Tinv


class PlyStrain(PlyStress):
    """An strain vector for a single lamina.

    This class is similar to the PlyStress class, but allows for input either
    engineering or true strain and corrects appropriately for transformation
    between the ply and laminate coordinate systems using the relationship

                                     simple
                                     matrix
                                       |
                            |e1 | = |1 0 0|   |   e1  |
                            |e2 |   |0 1 0| x |   e2  |
                            |g12|   |0 0 2|   |g12 / 2|
                              |                   |
                          engineering            true
                            strain              strain

    Note: ALL transformations of strain using the transformation matrix and its
          inverse must be performed using true strain

    See 'Mechanics Of Composite Materials' by Robert M. Jones, Section 2.6.
     """

    # define the simple matrix
    r = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 2]])

    inv_r = np.linalg.inv(r)

    def __init__(self,
                 vector,
                 angle,
                 ply_csys=True,
                 rad=False,
                 eng_strn=True):
        """Initialize the instance."""

        self.__eng_strn = eng_strn
        self._update(vector, angle, ply_csys, rad, eng_strn)

    def _update(self,
                vector=None,
                angle=None,
                ply_csys=True,  # default values are in ply c-sys
                rad=None,
                eng_strn=True):  # default input is engineering strain
        """Update PlyVector when values are changed."""

        if type(eng_strn) is not bool:
            raise TypeError("'eng_strn' must be a boolean.")
        elif eng_strn != self.__eng_strn:
            self.__eng_strn = eng_strn

        # if vector is in engineering strain, convert to true strain
        if vector is not None:
            if self.__eng_strn:
                vector = np.matmul(self.inv_r, vector)
            else:
                vector = self.__vector
        else:
            vector = self.__vector

        super()._update(vector, angle, ply_csys, rad)

    @property
    def strain(self):
        """The engineering strain vector in the ply coordinate system."""

        return self.__ply
