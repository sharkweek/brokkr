"""A fastener class for bolted joint calculations."""

import math
import numpy as np


class Fastener:
    """Fastener class containing an ID, location in space, and load vector."""

    def __init__(self,
                 ID,
                 xyz,
                 diameter,
                 E,
                 G,
                 length,
                 wt=1,
                 axis=[0, 0, 0]):
        """Initialize the instance."""

        self.ID = ID
        self.xyz = xyz
        self.diameter = diameter
        self.E = E
        self.G = G
        self.length = length
        self.force = [0, 0, 0]
        self.wt = wt
        self.axis = np.array(axis) / np.linalg.norm(axis)  # unit vector

    def __repr__(self):
        """Return the "official" Fastener string representation."""

        return "Fastener ID: %s\nLocation: %s\nForce: %s\nDiameter: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.diameter)

    @property
    def area(self):
        """Cross-sectional area of the fastener."""

        return (self.diameter ** 2) * math.pi / 4

    @area.setter
    def area(self, a):
        """Determine the diameter if the area is set manually."""

        self.diameter = math.sqrt(4 * a / math.pi)

    @property
    def stiffness(self):
        """Return stiffness of the fastener in each direction.

        The stiffness is calculated based on the diameter and material
        moduli of the fastener. The x-direction is the shaft axis.

        Note: this does not account for total behavior of the joints.
        """

        return {'x': self.E * self.area / self.length,
                'y': self.G * self.area / self.length,
                'z': self.G * self.area / self.length}

    @property
    def compliance(self):
        """Return the compliance of the fastener in each direction.

        The compliance in each direction is the inverse of the stiffness.
        """

        return {'x': self.length / (self.E * self.area),
                'y': self.length / (self.G * self.area),
                'z': self.length / (self.G * self.area)}
