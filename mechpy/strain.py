from numpy import array, ndarray

class Strain(ndarray):
    """
    Strain(*components, type='eng')

    A strain vector object. The vector follows standard tensor conventions.

    Parameters
    ----------
    *components : numpy.array_like, float, or int
        x, y, z, yz, zx, and xy components of strain tensor
    type : string
        'engr' for engineering strain
        'tens' for tensor strain

    Attributes
    ----------
    vector : numpy.ndarray
        6x1 array containing all strain values
    sigma1 : int or float
        normal strain in the 1-direction
    sigma2 : int or float
        normal strain in the 2-direction
    sigma3 : int or float
        normal strain in the 3-direction
    tau23 : int or float
        shear strain in the 23-plane (1-normal)
    tau31 : int or float
        shear strain in the 31-plane (2-normal)
    tau12 : int or float
        shear strain in the 12-plane (3-normal)
    principal : numpy.ndarray
        principal strains

    Notes
    -----
    * Strain may be input as either a series of numerical values or else a
      6x1 Numpy array_like object.
    * Because Strain is a subclass of the numpy.ndarray, see Numpy docs for
      numpy.ndarray for documentation.

    Examples
    --------
    >>> x = Strain(1, 2, 3, 4, 5, 6)
    >>> print(x)
    Strain([[1],
            [2],
            [3],
            [4],
            [5],
            [6]])

    >>> y = Strain([1, 2, 3, 4, 5, 6])
    >>> print(y)
    Strain([[1],
            [2],
            [3],
            [4],
            [5],
            [6]])
    """

    def __new__(cls, *strains, stype='engr'):

        # make sure six values are given
        try:
            new = array(strains).flatten().reshape(6, 1).view(cls)
        except ValueError:
            print("`strains` must be six values or else a 6x1 array-like obj")

        # assign strain type
        if stype != 'engr' and stype != 'tens':
            raise ValueError("`stype` must be 'engr' or 'tens'.")
        else:
            new.__stype = stype

        return ndarray.__new__(cls, shape=new.shape, dtype=new.dtype)

    @property
    def stype(self):
        """Type of strain (i.e. engineering or tensor)."""
        return self.__stype

    @stype.setter
    def stype(self, new_stype):
        if self.__stype == new_stype:
            pass

        # convert tensor strain to engineering strain
        elif new_stype == 'engr':
            self *= array((1, 1, 1, 2, 2, 2)).reshape(6, 1)

        # convert engineering strain to tensor
        elif new_stype == 'tens':
            self *= array((1, 1, 1, 0.5, 0.5, 0.5)).reshape(6, 1)

        else:
            raise ValueError("`stype` must be 'engr' or 'tens'.")
