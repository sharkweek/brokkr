"""Tools for composite laminates.

Notes
-----
All calculations are performed by using classical laminate theory (CLT) under
the following assumptions:

* Lamina are generally orthotropic.
* Lamina properties and values are assumed to be in the lamina coordinate
  system (i.e. 12-plane).
* Lamina angles (`Lamina.theta`) are in degrees, measured from the laminate
  x-axis to the lamina 1-axis.
* The laminate z-direction is upward, away from the bottom surface.
* Reported strain values are in engineering strain.
* All loads and moments are running values supplied as force or moment per unit
  width.
* Unit systems are consistent (i.e. SI, US, or Imperial).
* Equations and symbol conventions are per NASA-RP-1351 'Basic Mechanics of
  Laminated Composite Plates' and Jones' 'Mechanics Of Composite Materials'.

References
----------
.. [1] NASA-RP-1351, "Basic Mechanics of Laminated Composite Plates"
.. [2] Jones, Robert M. "Mechanics Of Composite Materials"

"""

from .lamina import Lamina
from .laminate import Laminate
