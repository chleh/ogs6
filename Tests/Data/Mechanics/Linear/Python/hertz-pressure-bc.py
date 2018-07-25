import OpenGeoSys

NODE_RELEASE_FRACTION = 0.5

class BC(OpenGeoSys.BoundaryCondition):
    def __init__(self):
        super(BC, self).__init__()
        self._restricted_nodes = {}
        self._first_node = None
        self._iteration = 0
        self._t_old = -1.0

    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        if self._t_old < t:
            # first iteration in current timestep
            print("restricted nodes in prev. ts", self._restricted_nodes)
            self._restricted_nodes.clear()
            self._t_old = t
            self._iteration = 0
        if self._first_node is None:
            self._first_node = node_id
        if self._first_node == node_id:
            self._iteration += 1

        x, y, z = coords

        y_deformed = y + primary_vars[1]
        y_top = BC.get_y_top(t)

        # print("y_top", y_top)
        if y_deformed > y_top:
            res = (True, y_top - y)
            self._restricted_nodes[node_id] = x
            print(f"[BC] {self._iteration:2} restr: y_deformed > y_top")
        elif y_deformed > y_top - 1e-8:
            # print("XXXXX def!", y, y_top, y_top - y)
            if self.is_new_restricted_node(node_id):
                res = (True, y_top - y)
                self._restricted_nodes[node_id] = x
                print(f"[BC] {self._iteration:2} new restr node")
            elif self.is_relaxed_restricted_node(x):
                res = (False, 0.0)
                print(f"[BC] {self._iteration:2} relax restr node")
            else:
                print(f"[BC] {self._iteration:2} generally y_deformed > y_top - 1e-8")
                res = (True, y_top - y)
        else:
            res =(False, 0.0)

        return res

    @staticmethod
    def get_y_top(t):
        return 1.0 - 0.005 * t

    def is_new_restricted_node(self, node_id):
        return node_id not in self._restricted_nodes
        # len(self._restricted_nodes) == 0 or \
        # x > max(self._restricted_nodes.items(), key=lambda p: p[1])[1]

    def is_relaxed_restricted_node(self, x):
        return len(self._restricted_nodes) != 0 and \
                x > (1.0 - NODE_RELEASE_FRACTION) * max(self._restricted_nodes.items(), key=lambda p:p[1])[1]


bc = BC()
