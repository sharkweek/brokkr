"""Stress Vector"""

from numpy import ndarray, array, zeros
from pint import UnitRegistry

ureg = UnitRegistry()


class Stress(ndarray):
    """A stress vector object. The vector follows standard tensor conventions.

    Parameters
    ----------
    *components : numpy.array_like, float, or int
        x, y, z, yz, zx, and xy components of stress tensor

    Attributes
    ----------
    sigma1 : int or float
        normal stress in the 1-direction
    sigma2 : int or float
        normal stress in the 2-direction
    sigma3 : int or float
        normal stress in the 3-direction
    tau23 : int or float
        shear stress in the 23-plane (1-normal)
    tau31 : int or float
        shear stress in the 31-plane (2-normal)
    tau12 : int or float
        shear stress in the 12-plane (3-normal)
    principal : numpy.ndarray
        principal stresses

    Notes
    -----
    * Stresses may be input as either a series of numerical values or else a
    6x1 numpy array_like object.
    * See documentation for base class numpy.ndarray for more information.

    Examples
    --------
    >>> x = Stress(1, 2, 3, 4, 5, 6)
    >>> print(x)
    Stress([[1., 6., 5.],
            [6., 2., 4.],
            [5., 4., 3.]])

    """

class Stress(ndarray):

    def __new__(cls, str):
        arr = zeros((3, 3))
        return super(Stress, cls).__new__(cls, shape=arr.shape,
                                          dtype=arr.dtype)

    def __init__(self, str):
        if len(str) == 6:
            self[0][0] = str[0]
            self[1][1] = str[1]
            self[2][2] = str[2]
            self[1][2] = str[3]
            self[2][1] = str[3]
            self[0][2] = str[4]
            self[2][0] = str[4]
            self[0][1] = str[5]
            self[1][0] = str[5]
        else:
            raise TypeError("Ah-ah-ah... you didn't say the magic word...")

    @property
    def Voigt(self):
        return array([[self[0][0]],
                      [self[1][1]],
                      [self[2][2]]
                      [self[1][2]]
                      [self[0][2]]
                      [self[0][1]]])

    @property
    def principal(self):
        """Principal stresses."""
        pass
