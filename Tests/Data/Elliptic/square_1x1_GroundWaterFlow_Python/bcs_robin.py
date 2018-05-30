
import OpenGeoSys

class BC(OpenGeoSys.BoundaryCondition):
    # def __init__(self):
    #     super().__init__()
    #     self._state = 0

    # def getDirichletBCValue(self, t, coords, node_id, primary_vars):
    #     print(">D", t, coords, node_id, primary_vars)
    #     return (True, 0.0)

    def getFlux(self, t, coords, primary_vars):
        p = primary_vars[0]
        p0 = 0.25
        a = 1.0
        # print(">N", t, coords, primary_vars, type(primary_vars))
        return (True, a * (p0 - p), [ -a ])

bc = BC()
