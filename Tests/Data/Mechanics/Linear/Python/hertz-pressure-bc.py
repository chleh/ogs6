import OpenGeoSys

NODE_RELEASE_FRACTION = 0.5

class BC(OpenGeoSys.BoundaryCondition):
    def __init__(self):
        super(BC, self).__init__()
        self._restricted_node_xs = set()
        self._t_old = -1.0

    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        if self._t_old < t:
            # first iteration in current timestep
            print("restricted nodes in prev. ts", self._restricted_node_xs)
            self._restricted_node_xs.clear()
            self._t_old = t

        x, y, z = coords

        y_deformed = y + primary_vars[1]
        y_top = BC.get_y_top(t)

        print("y_top", y_top)
        if y_deformed > y_top:
            res = (True, y_top - y)
            self._restricted_node_xs.add(x)
        elif y_deformed > y_top - 1e-8:
            print("XXXXX def!", y, y_top, y_top - y)
            if self.is_new_restricted_node(x):
                res = (True, y_top - y)
                self._restricted_node_xs.add(x)
            elif self.is_relaxed_restricted_node(x):
                res = (False, 0.0)
            else:
                res = (True, y_top - y)
        else:
            res =(False, 0.0)

        return res

    @staticmethod
    def get_y_top(t):
        return 1.0 - 0.005 * t

    def is_new_restricted_node(self, x):
        return len(self._restricted_node_xs) == 0 or x > max(self._restricted_node_xs)

    def is_relaxed_restricted_node(self, x):
        return len(self._restricted_node_xs) != 0 and x > (1.0 - NODE_RELEASE_FRACTION) * max(self._restricted_node_xs)


bc = BC()
