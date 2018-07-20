"""Calculate the flexibility of a fastener in a lap joint."""

import math


class Bolt(object):
    """A class representing a single fastner.

    Attributes:
        dia (float): shank diameter
        bolt_type (str): fastener type. Maybe either 'rivet' or 'bolt'
        e (float): Young's modulus
        c (float): flexibility
        k (float): stiffness
    """

    def __init__(self,
                 dia,
                 bolt_type,
                 e,
                 c=1):
        """Initialize the instance of Bolt."""

        self.dia = dia
        self.bolt_type = bolt_type
        self.e = e
        # optionals
        self.c = c  # sets c using c.setter

    @property
    def c(self):
        """Return the fastener flexibility based on the chosen method."""

        return self._c

    @c.setter
    def c(self, new_c):
        """Manually set the flexibility."""

        self._c = new_c
        self._flex_method = 'Manual'

    @property
    def k(self):
        """Return the fastener stiffness based on the flexibility method."""

        return 1/self._c

    @k.setter
    def k(self, new_k):
        """Manually set the stiffness."""

        self._c = 1 / new_k
        self._flex_method = 'Manual'

    @property
    def flex_method(self):
        """Return the fastener flexibility method."""

        return self._flex_method

    @flex_method.setter
    def flex_method(self, _):
        """Keep user from setting the flexibility method manually."""

        raise AttributeError("The flexibility method cannot be set mannually.")

    def swift(self, t1, t2, e1, e2):
        """Calculate the bolt stiffness using the Swift/Douglas method."""

        self._flex_method = 'Swift'

        # calculate flexibility
        self._c = (5 / (self.dia * self.e)
                   + 0.8 * (1 / (t1 * e1) + (1 / (t2 * e2))))

    def tate(self, t1, t2, e1, e2, nu_b):
        """Tate and Rosenfeld method"""
        self._flex_method = 'Tate/Rosenfeld'

        # calculate flexibility
        self._c = (1 / (self.e * t1) + 1 / (self.e * t2) + 1 / (e1 * t1)
                   + 1 / (e2 * t2)
                   + 32 / (9 * self.e * math.pi * self.dia**2)
                   * (1 + nu_b) * (t1 + t2)
                   + 8 / (5 * self.e * math.pi * self.dia**4)
                   * (t1**3 + 5 * t1**2 * t2 +
                      5 * t1 * t2**2 + t2**3))

    def boeing(self, t1, t2, e1, e2):
        """Boeing method"""
        self._flex_method = 'Boeing'

        self._c = (2 ** (t1 / self.dia) ** 0.85 / t1
                   * (1 / e1 + 3 / (8 * self.e))
                   + 2 ** (t2 / self.dia) ** 0.85 / t2
                   * (1 / e2 + 3 / (8 * self.e)))

    def huth(self, t1, t2, e1, e2, mat1, mat2, lap_type='single'):
        """Huth/Airbus method"""
        self._flex_method = 'Huth/Airbus'

        # single or double shear
        if lap_type == 'single':
            n = 1
        else:
            n = 2

        # determine a and b values based on joint configuration
        if (mat1 == 'metal' and mat2 == 'metal' and self.bolt_type == 'bolt'):
            a = 2/3
            b = 3.0

        elif (mat1 == 'metal' and mat2 == 'metal' and
              self.bolt_type == 'rivet'):
            a = 2/5
            b = 2.2

        elif (mat1 == 'composite' and mat2 == 'composite' and
              self.bolt_type == 'bolt'):
            a = 2/3
            b = 4.2

        # calculate flexibility
        self._c = (((t1 + t2) / (2 * self.dia))**a * (b / n)
                   * (1 / (t1 * e1)
                      + 1 / (n * t2 * e2)
                      + 1 / (2 * t1 * self.e)
                      + 1 / (2 * n * t2 * self.e)))

    def vought(self, t1, t2):
        """Vought method

        Note: only applies for aluminum sheets joined by steel fasteners"""
        self._flex_method = 'Vought'

        # determine overall grip length
        grip = t1 + t2

        # calculate epsilon based on grip-length-to-diameter ratio
        if (grip / self.dia) < 0.65:
            epsilon = 1

        elif (grip / self.dia) > 0.9:
            epsilon = 1.29 * (grip / self.dia)

        # calculate flexibility
        self._c = 56 * (epsilon * (grip / self.dia) * (1 / t1)
                        + epsilon * (grip / self.dia) * (1 / t2))

    def grumman(self, t1, t2, e1, e2):
        """Grumman method"""
        self._flex_method = 'Grumman'

        # calculate flexibility
        self._c = ((t1 + t2)**2 / (self.e * self.dia**3)
                   + 3.7 * (1 / (e1 * t1) + 1 / (e2 + t2)))
