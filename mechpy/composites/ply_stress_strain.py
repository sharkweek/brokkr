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

        self.__typecheck(vector, angle, ply_csys, rad)

        vector = np.array(vector).reshape((3, 1))

        self.__vector = 0
        self.__angle = 0
        self.__ply_csys = 0
        self.__rad = 0

        self.__update(vector, angle, ply_csys, rad)

    def __update(self,
                 vector=None,
                 angle=None,
                 ply_csys=True,  # default values are in ply c-sys
                 rad=None):
        """Update PlyStress when values are changed."""

        # check if angle has changed
        if angle is not None and angle != self.__angle:
            self.__angle = angle
            transform = True
        else:
            transform = False

        # check if angle units have changed
        if rad is not None and rad != self.__rad:
            self.__rad = rad
            transform = True
        else:
            transform = False

        # check if vector has changed
        if vector is not None and not np.array_equal(vector, self.__vector):
            self.__vector = vector
            transform = True
        else:
            transform = False

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

            # update the vector attribute
            # the vector attribute is ALWAYS stored in the ply c-sys
            if not ply_csys:
                self.__vector = np.matmul(self.__Tinv, vector)

        del transform, m, n

    def __typecheck(self,
                    vector=None,
                    angle=None,
                    ply_csys=None,
                    rad=None):
        """Check the data types for all values in the instance."""

        if vector is not None:
            if (type(vector) != list and type(vector) != np.ndarray):
                raise TypeError("""'vector' must be a list of format [x, y, z]
                                or a 3x1 numpy array.""")
            elif len(vector) != 3:
                raise TypeError("""'vector' must be a list of format [x, y, z]
                                or a 3x1 numpy array.""")

        if ply_csys is not None and type(ply_csys) != bool:
            raise TypeError("'ply_csys' must be a boolean.")

        if rad is not None and type(rad) != bool:
            raise TypeError("'rad' must be a boolean")

        if angle is not None:
            try:
                float(angle)
            except TypeError:
                print("'angle' must be a number.")
                raise

    @property
    def angle(self):
        """Lamina orientation angle w.r.t. the laminate c-sys."""

        return self.__angle

    @angle.setter
    def angle(self, new_angle):
        """Lamina orientation angle w.r.t. the laminate c-sys."""

        self.__typecheck(angle=new_angle)
        self.__update(angle=new_angle)

    @property
    def rad(self):
        """Units flag for orientation angle."""

        return self.__rad

    @rad.setter
    def rad(self, new_rad):
        """Units flag for orientation angle."""

        self.__typecheck(rad=new_rad)
        self.__update(rad=new_rad)

    @property
    def ply_s(self):
        """Lamina stress w.r.t. the ply c-sys."""

        return self.__vector

    @ply_s.setter
    def ply_s(self, new_vector):
        """Lamina stress w.r.t. the ply c-sys."""

        self.__typecheck(vector=new_vector)

        # convert to numpy array if a list is provided
        if type(new_vector) == list:
            new_vector = np.array(new_vector).reshape((3, 1))

        self.__update(vector=new_vector)

    @property
    def lam_s(self):
        """Lamina stress w.r.t. the laminate c-sys."""

        return np.matmul(self.__T, self.__vector)

    @lam_s.setter
    def lam_s(self, new_vector):
        """Lamina stress w.r.t. the laminate c-sys."""

        self.__typecheck(vector=new_vector)

        if type(new_vector) == list:
            new_vector = np.array(new_vector).reshape((3, 1))

        self.__update(vector=new_vector, ply_csys=False)

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

    This class is similar to the PlyStress class, but corrects for using
    engineering strain instead of tensor strain. The relationship between them
    is shown below.

                                     simple
                                   matrix (r)
                                       |
                            |e1 | = |1 0 0|   |   e1  |
                            |e2 |   |0 1 0| x |   e2  |
                            |g12|   |0 0 2|   |g12 / 2|
                              |                   |
                          engineering           tensor
                            strain              strain

    Note: ALL transformations of strain using the transformation matrix must be
    performed using tensor strain.

    See 'Mechanics Of Composite Materials' by Robert M. Jones, Section 2.6.
     """

    # define the simple matrix
    r = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 2]])

    r_inv = np.linalg.inv(r)

    def __init__(self,
                 vector,
                 angle,
                 ply_csys=True,
                 rad=False):
        """Initialize the instance."""

        # make sure vector is the correct type
        self._PlyStress__typecheck(vector=vector)

        vector = np.array(vector).reshape((3, 1))

        vector = np.matmul(self.r_inv, vector)  # Jones Eq 2.78 and 2.79

        super().__init__(vector, angle, ply_csys, rad)

    @property
    def ply_s(self):
        """Lamina strain w.r.t. the ply c-sys."""

        # return the strain as engineering strain
        return np.matmul(self.r, self._PlyStress__vector)

    @ply_s.setter
    def ply_s(self, new_strain):
        """Lamina strain w.r.t. the ply c-sys."""

        self._PlyStress__typecheck(vector=new_strain)

        # convert to numpy array if a list is provided
        if type(new_strain) == list:
            new_strain = np.array(new_strain).reshape((3, 1))

        self._PlyStress__vector = np.matmul(self.r_inv, new_strain)

    @property
    def lam_s(self):
        """Lamina strain w.r.t. the ply c-sys."""

        # transform to laminate c-sys
        eng_strn = np.matmul(self._PlyStress__T,
                             self._PlyStress__vector)

        # return as engineering strain
        return np.matmul(self.r, eng_strn)

    @lam_s.setter
    def lam_s(self, new_strain):
        """Lamina strain w.r.t. the laminate c-sys."""

        self._PlyStress__typecheck(vector=new_strain)

        # convert to numpy array if a list is provided
        if type(new_strain) == list:
            new_strain = np.array(new_strain).reshape((3, 1))

        # convert to tensor strain
        new_strain = np.matmul(self.r_inv, new_strain)

        self._PlyStress__update(vector=new_strain)
