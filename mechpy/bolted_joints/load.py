"""A load class."""

class Load:
    """Load class containing location, force, and moment vectors."""

    def __init__(self,
                 xyz=[0, 0, 0],
                 force=[0, 0, 0],
                 moment=[0, 0, 0]):
        """Initialize a Load instance."""

        self.xyz = xyz
        self.force = force
        self.moment = moment

    def __repr__(self):
        """Return the "official" Load string representation."""

        return "Load ID: %s\nLocation: %s\nForce: %s\nMoment: %s" \
               % (self.ID,
                  self.xyz,
                  self.force,
                  self.moment)