import CoolProp.CoolProp as CP



class Fluid():
    def __init__(self, fluid) -> None:
        #TODO: check that the fluid is actually valid
        self.name = fluid
        self.properties_dict = {
            "pressure": "P",
            "temperature": "T",
            "density": "D",
            "enthalpy": "H",
            "entropy": "S",
            "cpmass": "Cpmass",
            "cvmass": "Cvmass",
            "sound_speed": "A",
            "viscosity": "V",
            "prandtl": "Prandtl"
        }
    
    def __getattr__(self, property_name):

        if property_name in self.properties_dict.keys():
            def property_method(**kwargs):
                if len(kwargs) != 2:
                    raise ValueError("Two state properties are required")

                states = list(kwargs.keys())
                values = list(kwargs.values())
                try:
                    property = CP.PropsSI(self.properties_dict[property_name], states[0], values[0], states[1], values[1], self.name)
                    return property
                except Exception as e:
                    raise ValueError(f"Error in calculating property: {e}")

            return property_method

        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{property_name}'")