import OpenGeoSys

def get_y_top(t):
    return 1.0 - 0.09 * t

class BC(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        y = coords[1]

        # return (True, -0.1 * y)

        y_deformed = y + primary_vars[1]
        y_top = get_y_top(t)
        print("y_top", y_top)
        if y_deformed > y_top - 1e-6:
            print("XXXXX def!", y, y_top, y_top - y)
            # return (True, y_top - y)
            return (True, y_top - y)
            # return (True, -0.01 * y)

        return (False, 0.0)

    # def getFlux(self, t, coords, primary_vars):
    #     print(">N", t, coords, primary_vars, type(primary_vars))
    #     return (True, 0.0, [ 1.0, 1.0 ])

bc = BC()
