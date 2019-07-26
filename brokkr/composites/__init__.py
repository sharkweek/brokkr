"""Tools for composite laminates.

When no units are supplied, they are assumed to be in the master unit system
``brokkr.config.DEFAULT_USYS``.

All calculations are performed by using classical laminate theory (CLT) under
the following assumptions:

* Lamina are generally orthotropic.
* Lamina properties and values are assumed to be in the lamina coordinate
  system (i.e. 12-plane).
* Lamina angles (``Lamina.theta``) are in degrees, measured from the laminate
  x-axis to the lamina 1-axis.
* The laminate z-direction is upward, away from the bottom surface.
* Reported strain values are in engineering strain.
* All loads and moments are running values supplied as force or moment per unit
  width.
* Unit systems are consistent (i.e. SI, US, or Imperial).
* Equations and symbol conventions are per NASA-RP-1351 [#1]_ and Jones [#2]_.

References
----------
.. [#1] Nettles, A.T., "Basic Mechanics of Laminated Composite Plates", NASA
   Marshall Space Flight Center, Alabama, NASA-RP-1351, October 1994
.. [#2] Jones, Robert M. "Mechanics Of Composite Materials", Taylor & Francis,
   Inc., 1999

"""

from brokkr.composites.lamina import Lamina, Ply
from brokkr.composites.laminate import Laminate

__all__ = ['Lamina', 'Ply', 'Laminate']
