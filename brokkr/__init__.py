import brokkr.bolted_joints
import brokkr.composites
from brokkr.tensor import *
from brokkr.vector import *
import brokkr.mech_math
from unyt import UnitSystem, lbf, psi, inch, ft
from unyt.dimensions import length, force, dimensionless
from unyt.unit_registry import default_unit_registry

UREG = default_unit_registry
UREG.add(
    symbol='strain',
    base_value=1.0,
    dimensions=length / length,
    tex_repr=r"\rm{\varepsilon}",
    prefixable=True
)
# add percent for coefficient of hygroscopic expansion
UREG.add(
    symbol='percent',
    base_value=1.0,
    dimensions=dimensionless,
    tex_repr=r'\rm{%}',
    prefixable=False
)
# add torque units
UREG.add(
    symbol='N_m',
    base_value=1.0,
    dimensions=force * length,
    tex_repr=r'\rm{N \cdot m}',
    prefixable=True
)
UREG.add(
    symbol='in_lbf',
    base_value=(1 * inch * lbf).to('N*m').v.item(0),
    dimensions=force * length,
    tex_repr=r'\rm{in \cdot lb_{F}}',
    prefixable=True
)
UREG.add(
    symbol='ft_lbf',
    base_value=(1 * ft * lbf).to('N*m').v.item(0),
    dimensions=force * length,
    tex_repr=r'\rm{ft \cdot lb_{f}}'
)

us_unit_system = UnitSystem(
    name='us',
    length_unit='inch',
    mass_unit='lb',
    time_unit='s',
    temperature_unit='degF',
    angle_unit='degree',
    registry=UREG
)
us_unit_system['force'] = lbf
us_unit_system['pressure'] = psi

si_unit_system = UnitSystem(
    name='si',
    length_unit='m',
    mass_unit='kg',
    time_unit='s',
    temperature_unit='degC',
    angle_unit='rad',
    registry=UREG
)

USYS = us_unit_system
