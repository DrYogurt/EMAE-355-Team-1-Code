import CoolProp.CoolProp as CP
from pint import UnitRegistry



class Fluid():
    def __init__(self, fluid, ureg=None) -> None:
        if ureg is None:
            self.ureg = UnitRegistry()
        else:
            self.ureg = ureg
        self.Q_ = self.ureg.Quantity
        #TODO: check that the fluid is actually valid
        self.name = fluid
        self.properties_dict = {
            "pressure": ("P", self.ureg.pascal),
            "temperature": ("T", self.ureg.kelvin),
            "density": ("D", self.ureg.kilogram / self.ureg.meter ** 3),
            "enthalpy": ("H", self.ureg.joule / self.ureg.kilogram),
            "entropy": ("S", self.ureg.joule / (self.ureg.kilogram * self.ureg.kelvin)),
            "cpmass": ("Cpmass", self.ureg.joule / (self.ureg.kilogram * self.ureg.kelvin)),
            "cvmass": ("Cvmass", self.ureg.joule / (self.ureg.kilogram * self.ureg.kelvin)),
            "sound_speed": ("A", self.ureg.meter / self.ureg.second),
            "viscosity": ("V", self.ureg.pascal * self.ureg.second),
            "prandtl": ("Prandtl", self.ureg.dimensionless)
        }
    
    def __getattr__(self, property_name):

        if property_name in self.properties_dict.keys():
            def property_method(**kwargs):
                if len(kwargs) != 2:
                    raise ValueError("Two state properties are required")

                states = list(kwargs.keys())
                values = list(kwargs.values())
                values = [value if type(value) is float else value.magnitude for value in values]
                try:
                    property = CP.PropsSI(self.properties_dict[property_name][0], states[0], values[0], states[1], values[1], self.name)
                    return self.Q_(property, self.properties_dict[property_name][1])
                except Exception as e:
                    raise ValueError(f"Error in calculating property: {e}")

            return property_method

        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{property_name}'")