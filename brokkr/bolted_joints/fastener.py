"""A fastener class for bolted joint calculations."""

from unyt import (
    unyt_array
)
from unyt.dimensions import (
    length,
    pressure
)

from brokkr.config import DEFAULT_USYS

__all__ = ['Fastener']

class Fastener:
    """A fastener.

    Attributes
    ----------
    xyz : 3x1 unyt.unyt_array
        location of the fastener (``dim: length``)
    dia : float
        shakn diameter (``dim: length``)
    E : float
        Young's modulus of the fastener material (``dim: (mass)
        /((length)*(time)**2)``)
    G : float
        shear modulus of the fastener material (``dim: (mass)
        /((length)*(time)**2)``)
    length : float
        shank length (``dim: length``)
    flex_method : {'Huth-bolt-metal', 'Huth-rivet-metal', 'Huth-bolt-graphite',
                   'Boeing', 'Swift', 'Tate', 'Grumman'}
        method for calculating fastener flexibility
        (``default='Huth-bolt-metal'``)
    usys : unyt.UnitSystem
        fastener unit system (``default=DEFAULT_USYS``)
    """

    _flex_types = ['Huth-bolt-metal', 'Huth-rivet-metal', 'Huth-bolt-graphite',
                   'Boeing', 'Swift', 'Tate', 'Grumman']

    _param_dims = {
        'xyz': (length, ),
        'dia': (length, ),
        'E': (pressure, ),
        'G': (pressure, ),
        'length': (length, )
    }

    _param_limits = {
        'dia': {'mn': 0, 'mx': None, 'condition': 'g'},
        'E': {'mn': 0, 'mx': None, 'condition': 'g'},
        'G': {'mn': 0, 'mx': None, 'condition': 'g'},
        'length': {'mn': 0, 'mx': None, 'condition': 'g'}
    }

    __slots__ = ['xyz', 'dia', 'E', 'G', 'length', 'c', 'flex_method', 'usys']

    def __init__(self,
                 xyz,
                 dia,
                 E,
                 G,
                 length,
                 flex_method='Huth-bolt-metal',
                 usys=DEFAULT_USYS):
        """Create a `Fastener`."""

        self.xyz = xyz
        self.dia = dia
        self.E = E
        self.G = G
        self.length = length

    @property
    def x(self):
        """The x-coordinate of the fastener."""

        return
