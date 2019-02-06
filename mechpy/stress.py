"""Stress Vector"""

from numpy import ndarray, array


class Stress(ndarray):
    """
    Stress(*components)

    A stress vector object. The vector follows standard tensor conventions.

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
    Stress([[1],
            [2],
            [3],
            [4],
            [5],
            [6]])

    >>> y = Stress([1, 2, 3, 4, 5, 6])
    >>> print(y)
    Stress([[1],
            [2],
            [3],
            [4],
            [5],
            [6]])
    """

    def __new__(cls, *stresses):
        try:
            new = array(stresses).flatten().reshape(6, 1)
        except ValueError:
            print("`stresses` must be six values or else a 6x1 array-like obj")
        
        return ndarray.__new__(cls, shape=new.shape, dtype=new.dtype)

    @property
    def principal(self):
        """Principal stresses."""

    @property
    def sigma1(self):
        """Normal Stress in the x-direction."""
        return self[0]

    @property
    def sigma2(self):
        """Normal Stress in the y-direction."""
        return self[1]
    
    @property
    def sigma3(self):
        """Normal Stress in the z-direction."""
        return self[2]

    @property
    def tau23(self):
        """Shear stress in the yz-plane (x-normal)."""
        return self[3]

    @property
    def tau31(self):
        """Shear stress in the zx-plane (y-normal)."""
        return self[4]

    @property
    def tau12(self):
        """Shear stress in the xy-plane (z-normal)."""
        return self[5]
