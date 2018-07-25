import OpenGeoSys

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
        if y_deformed >= y_top:
            print("XXXXX def!", y, y_top, y_top - y)
            if self.is_new_restricted_node(x):
                res = (True, y_top - y)
                self._restricted_node_xs.add(x)
            elif self.is_rightmost_restricted_node(x):
                res = (False, 0.0)
            else:
                res = (True, y_top - y)
        else:
            res =(False, 0.0)

        self._t_old = t
        return res

    @staticmethod
    def get_y_top(t):
        return 1.0 - 0.005 * t

    def is_first_iteration(self, x):
        # TODO better detection
        # return x in self._restricted_node_xs
        return len(self._restricted_node_xs) == 0 or x > max(self._restricted_node_xs)

    def is_new_restricted_node(self, x):
        return len(self._restricted_node_xs) == 0 or x > max(self._restricted_node_xs)

    def is_rightmost_restricted_node(self, x):
        return len(self._restricted_node_xs) != 0 and x == max(self._restricted_node_xs)


bc = BC()
