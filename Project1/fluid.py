import CoolProp.CoolProp as CP
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity

class Fluid():
    def __init__(self, fluid) -> None:
        #TODO: check that the fluid is actually valid
        self.name = fluid
        self.properties_dict = {
            "pressure": ("P", ureg.pascal),
            "temperature": ("T", ureg.kelvin),
            "density": ("D", ureg.kilogram / ureg.meter ** 3),
            "enthalpy": ("H", ureg.joule / ureg.kilogram),
            "entropy": ("S", ureg.joule / (ureg.kilogram * ureg.kelvin)),
            "cpmass": ("Cpmass", ureg.joule / (ureg.kilogram * ureg.kelvin)),
            "cvmass": ("Cvmass", ureg.joule / (ureg.kilogram * ureg.kelvin)),
            "sound_speed": ("A", ureg.meter / ureg.second),
            "viscosity": ("V", ureg.pascal * ureg.second),
            "prandtl": ("Prandtl", ureg.dimensionless)
        }
    
    def __getattr__(self, property_name):

        if property_name in self.properties_dict.keys():
            def property_method(**kwargs):
                if len(kwargs) != 2:
                    raise ValueError("Two state properties are required")

                states = list(kwargs.keys())
                values = list(kwargs.values())
                try:
                    return Q_(CP.PropsSI(self.properties_dict[property_name][0], states[0], values[0], states[1], values[1], self.name),
                              self.properties_dict[property_name][1])
                except Exception as e:
                    raise ValueError(f"Error in calculating property: {e}")

            return property_method

        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{property_name}'")