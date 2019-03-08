"""Define vector classes."""

from numpy import array
from unyt.dimensions import force, length
from unyt.array import unyt_array

__all__ = ['Vector', 'Force', 'Moment', 'Displacement']


class Vector(unyt_array):
    """Base 3x1 vector array.

    Parameters
    ----------
    matrix : numpy.array_like
        an array-like object with three items
    units : string or unyt.Unit object
        the physical units for the vector

    Notes
    -----
    If ``units`` is supplied as a string, it must conform to the same
    formatting rules as the ``unyt.array.unyt_array.to()`` method.

    """

    def __new__(cls, matrix, units):
        """Create `Vector` instance."""
        new = array(matrix)
        if new.size != 3:
            raise TypeError("`matrix` must have six values")
        else:
            vector = unyt_array(new, units).view(cls)
            vector.resize((3, 1))
            return vector

    @property
    def resultant(self):
        """The resultant magnitude."""

        return (self[0]**2 + self[1]**2 + self[2]**2)**(0.5)

    def has_dimension(self, dim):
        """Check if dimensions of vector match `dim`.

        Parameters
        ----------
        dim : ``unyt.dimensions`` dimension
            the dimension to compare the instance against

        Returns
        -------
        bool
            True if instance's dimension matches ``dim``; False otherwise.

        """

        return self.units.dimensions == dim


class Force(Vector):
    """A 3D force vector."""

    def __new__(cls, forces, units='lbf'):
        """Create `Force` instance."""
        new = super().__new__(cls, forces, units)
        # check for force dim
        if new.has_dimension(force):
            return new
        else:
            raise TypeError("`units` must be force units")


class Moment(Vector):
    """A 3D moment vector."""

    def __new__(cls, moments, units='inch*lbf'):
        """Create `Moment` instance."""
        new = super().__new__(cls, moments, units)
        # check for moment dim
        if new.has_dimension(length * force):
            return new
        else:
            raise TypeError("`units` must be moment units")


class Displacement(Vector):
    """A 3D distance vector."""

    def __new__(cls, displacements, units='inch'):
        """Create `Displacement` instance."""
        new = super().__new__(cls, displacements, units)
        # check for length dim
        if new.has_dimension(length):
            return new
        else:
            raise TypeError("Units must be length units.")
