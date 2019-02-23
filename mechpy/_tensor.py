from numpy import array, empty
import unyt as u
from unyt.array import unyt_array

# set unit regsitry
us_unit_system = u.UnitSystem(name='US',
                              length_unit='inch',
                              )

# add strain unit
ureg = u.UnitRegistry()
ureg.add('strain',
         base_value=1.0,
         dimensions=u.dimension.dimensionless,
         text_repr=r"\rm{\varepsilon}",
         prefixable=True)


class Tensor(unyt_array):
    """Base tensor."""

    def __new__(cls, matrix, dim, units):
        """Create `Tensor` instance."""
        return unyt_array(empty((dim, dim)), units, ureg).view(cls)


class Stress(Tensor):
    """Stress tensor."""

    def __new__(cls, stresses):
        """Create `Stress` instance."""
        return super().__new__(cls, stresses, 3, 'psi')

    def __init__(self, stresses):
        """Initialize `Stress` instance."""
        self[0][0] = stresses[0]
        self[1][1] = stresses[1]
        self[2][2] = stresses[2]
        self[1][2] = stresses[3]
        self[2][1] = stresses[3]
        self[0][2] = stresses[4]
        self[2][0] = stresses[4]
        self[0][1] = stresses[5]
        self[1][0] = stresses[5]

    def principal(self):
        #http://www.engapplets.vt.edu/Mohr/java/nsfapplets/MohrCircles2-3D/Theory/theory.htm
        pass


class Strain(Tensor):
    """Strain tensor."""

    def __new__(cls, strains, stype='engr'):
        """Create `Strain` instance object."""
        new = super().__new__(cls, strains, 3, 'strain')
        new.stype = stype
        return new

    def __init__(self, strains):
        """Initialize `Stress` instance."""
        self[0][0] = strains[0]
        self[1][1] = strains[1]
        self[2][2] = strains[2]
        self[1][2] = strains[3]
        self[2][1] = strains[3]
        self[0][2] = strains[4]
        self[2][0] = strains[4]
        self[0][1] = strains[5]
        self[1][0] = strains[5]
