"""Define stress/strain tensors."""

from numpy import array, empty
from unyt.dimensions import pressure, length
from unyt.array import unyt_array

__all__ = ['BaseTensor', 'StrainTensor', 'StressTensor']


class BaseTensor(unyt_array):
    """Base 3x3 symmetric tensor.

    Parameters
    ----------
    matrix : numpy.array_like
        an array-like object with six values
    units : string or unyt.Units
        the physical units for the tensor

    """

    def __new__(cls, matrix, units):
        """Create `BaseTensor` instance."""
        if array(matrix).size != 6:
            raise TypeError("`matrix` must have six values")
        else:
            new = unyt_array(empty((3, 3)), units).view(cls)
            return new

    def __init__(self, matrix):
        """Initialize `BaseTensor` instance."""
        for i, j in enumerate(array(matrix).ravel()):
            self.set_v_item(i, j)

    @property
    def voigt(self):
        """The Voigt representation of the array."""
        # create the indexing tuple for returning specific items in self
        voigt = [0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]

        return self[voigt].reshape((6, 1))

    def set_v_item(self, vindex, new_val):
        """Set an item in-place using Voigt index.

        Parameters
        ----------
        vindex : int
            the Voigt index of the item to set
        new_val : float
            the new value for the specified ``vindex``

        """

        i = [((0, 0),),
             ((1, 1),),
             ((2, 2),),
             ([1, 2], [2, 1]),
             ([2, 0], [0, 2]),
             ([0, 1], [1, 0])]

        for j, k in i[vindex]:
            self[j, k] = new_val

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


class StressTensor(BaseTensor):
    """Stress tensor.

    A subset of the tensor class, the stress tensor requires that ``units`` be
    in pressure and provides additional, stress-specific functionality.

    """

    def __new__(cls, stresses, units='psi'):
        """Create `StressTensor` instance."""
        new = super().__new__(cls, stresses, units)
        # check for pressure dim
        if new.has_dimension(pressure):
            return new
        else:
            raise TypeError("`units` must be pressure units.")

    def vonMises(self):
        r"""The equivalent von Mises stress of the stress tensor.

        The von Mises stress :math:`\sigma_{v}` is defined:

        .. math::
            \sigma_{v} = \sqrt{
            \frac{1}{2} \left[ \left( \sigma_{11} - \sigma_{22} \right) ^ 2
            + \left( \sigma_{22} - \sigma_{33} \right) ^ 2
            + \left( \sigma_{33} - \sigma_{11} \right) ^ 2
            + 6 \left( \tau_{12}^2 + \tau_{23}^2 + \tau_{31}^2 \right) \right]
            }

        """
        s = self.voigt
        return ((1 / 2) * ((s[0] - s[1])**2
                           + (s[1] - s[2])**2
                           + (s[2] - s[0])**2
                           + 6 * (s[3]**2 + s[4]**2 + s[5]**2))) ** (1 / 2)

    @property
    def principal(self):
        r"""Array of principal stresses.

        The array format is:

        .. math::
           \left[ \begin{matrix}
           (\sigma_{xy})_1 & (\sigma_{xy})_2 & (\tau_{xy})_\text{max} \\
           (\sigma_{yz})_1 & (\sigma_{yz})_2 & (\tau_{yz})_\text{max} \\
           (\sigma_{zx})_1 & (\sigma_{zx})_2 & (\tau_{zx})_\text{max}
           \end{matrix} \right]

        """
        s = [self.voigt[i] for i in range(6)]

        # calculate centers
        cxy = (s[0] + s[1]) / 2
        cyz = (s[1] + s[2]) / 2
        czx = (s[2] + s[0]) / 2

        # calculate radii
        rxy = ((s[0] - s[1])**2 + (2*s[3])**2)**(1/2) / 2
        ryz = ((s[1] - s[2])**2 + (2*s[4])**2)**(1/2) / 2
        rzx = ((s[2] - s[0])**2 + (2*s[5])**2)**(1/2) / 2

        return array([[cxy + rxy, cxy - rxy, rxy],
                      [cyz + ryz, cyz - ryz, ryz],
                      [czx + rzx, czx - rzx, rzx]])


class StrainTensor(BaseTensor):
    """A strain tensor.

    A subset of the tensor class, the strain tensor requires that ``units`` be
    in units of length/length and provides additional, strain-specific
    functionality.

    Parameters
    ----------
    stype : {'engineering', 'tensor'}
        the type of shear strain

    """

    def __new__(cls, strains, units='inch/inch', stype='engineering'):
        """Create `StrainTensor` instance object."""

        new = super().__new__(cls, strains, units)
        new.__stype = stype

        # check dimensions
        if new.has_dimension(length / length):
            return new
        else:
            raise TypeError(
                "`units` must be in strain (length/length or dimensionsless)"
            )

    @property
    def stype(self):
        r"""StrainTensor type.

        Must be 'engineering' or 'tensor' strain.

        Notes
        -----
        Only affects shear strains (i.e. :math:`\varepsilon_4`,
        :math:`\varepsilon_5`, and :math:`\varepsilon_6`). Using Voigt
        notation, the relationship is defined:

        .. math::
           \underbrace{ \left[ \begin{matrix}
           \varepsilon_1 \\ \varepsilon_2 \\ \varepsilon_3 \\ \varepsilon_4 \\
           \varepsilon_5 \\ \varepsilon_6
           \end{matrix} \right] }_\text{contracted notation} =
           \underbrace{ \left[ \begin{matrix}
           \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\
           \gamma_{23} \\ \gamma_{31} \\ \gamma_{12}
           \end{matrix} \right] }_\text{engineering strain} =
           \underbrace{ \left[ \begin{matrix}
           \varepsilon_{11} \\ \varepsilon_{22} \\ \varepsilon_{33} \\
           2 \varepsilon_{23} \\ 2 \varepsilon_{31} \\ 2 \varepsilon_{12}
           \end{matrix} \right] }_\text{engineering strain}

        """

        return self.__stype

    @stype.setter
    def stype(self, new_stype):
        if new_stype not in ['engineering', 'tensor']:
            raise TypeError("`stype` must be 'engineering' or 'tensor'")

        elif new_stype != self.__stype:
            self.__stype = new_stype

            # to modify in-place, individual values have to be modified
            if new_stype == 'engineering':
                for i in (3, 4, 5):
                    self.set_v_item(i, self.voigt[i] * 2)

            elif new_stype == 'tensor':
                for i in (3, 4, 5):
                    self.set_v_item(i, self.voigt[i] / 2)
