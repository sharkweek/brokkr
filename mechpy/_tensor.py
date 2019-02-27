"""Define stress/strain tensors."""

from numpy import array, empty
import unyt as u
from unyt.array import unyt_array

# set unit regsitry
us_unit_system = u.UnitSystem(name='US',
                              length_unit='inch',
                              mass_unit='lb',
                              time_unit='s',
                              temperature_unit='degF',
                              angle_unit='deg')

# add strain unit
ureg = u.UnitRegistry()
ureg.add('strain',
         base_value=1.0,
         dimensions=u.dimensions.dimensionless,
         tex_repr=r"\rm{\varepsilon}",
         prefixable=True)
u.unit_registry.default_unit_registry = ureg


class Tensor(unyt_array):
    """Base 3x3 symmetric tensor.

    Parameters
    ----------
    matrix : numpy.array_like
        an array-like object with six values
    units : string or unyt.Units
        the physical units for the tensor

    """

    def __new__(cls, matrix, units):
        """Create `Tensor` instance."""
        if array(matrix).size != 6:
            raise TypeError("`matrix` must have six values")
        else:
            new = unyt_array(empty((3, 3)), units, ureg).view(cls)
            return new

    def __init__(self, matrix):
        """Initialize `Tensor` instance."""
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


class Stress(Tensor):
    """Stress tensor.

    A subset of the tensor class, the stress tensor requires that units be
    a pressure and provides additional, stress-specific functionality.
    Acceptable pressure units are ``'Ba'``, ``'Pa'``, ``'atm'``, and ``'psi'``.

    """

    def __new__(cls, stresses, units='psi'):
        """Create `Stress` instance."""
        # check for pressure dim
        if units not in ureg.list_same_dimensions(u.psi):
            raise TypeError(
                "Units must be pressure untis: 'Ba', 'Pa', 'atm', or 'psi'"
            )

        new = super().__new__(cls, stresses, units)
        return new

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


class Strain(Tensor):
    """Strain tensor."""

    def __new__(cls, strains, stype='engr'):
        """Create `Strain` instance object."""

        new = super().__new__(cls, strains, 3, 'strain')
        new.stype = stype
        return new

    @property
    def stype(self):
        r"""Strain type.

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
