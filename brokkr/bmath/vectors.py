"""Define vector classes."""

from numpy import array
from unyt.dimensions import force, length
from unyt.array import unyt_array

from brokkr.config import DEFAULT_USYS
from brokkr.core.bases import BaseVector

__all__ = ['ForceVector', 'MomentVector', 'DisplacementVector']


class ForceVector(BaseVector):
    """A 3D force vector."""

    def __new__(cls, forces, units='lbf'):
        """Create `ForceVector` instance."""
        new = super().__new__(cls, forces, units, usys)
        # check for force dim
        if new.has_dimension(force):
            return new
        else:
            raise AttributeError("`units` must be force units")


class MomentVector(BaseVector):
    """A 3D moment vector."""

    def __new__(cls, moments, units='inch*lbf'):
        """Create `MomentVector` instance."""
        new = super().__new__(cls, moments, units)
        # check for moment dim
        if new.has_dimension(length * force):
            return new
        else:
            raise AttributeError("`units` must be moment units")


class DisplacementVector(BaseVector):
    """A 3D distance vector."""

    def __new__(cls, displacements, units='inch'):
        """Create `DisplacementVector` instance."""
        new = super().__new__(cls, displacements, units)
        # check for length dim
        if new.has_dimension(length):
            return new
        else:
            raise AttributeError("Units must be length units.")
