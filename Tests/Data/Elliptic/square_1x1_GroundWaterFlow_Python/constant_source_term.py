
import OpenGeoSys

# Constant source term
class ConstantSourceTerm(OpenGeoSys.SourceTerm):
    def getFlux(self, t, coords, primary_vars):
        x, y, z = coords
        value = 1.0
        return (value, 0.0)

# instantiate source term object referenced in OpenGeoSys' prj file
constant_source_term = ConstantSourceTerm()
