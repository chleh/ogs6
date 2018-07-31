from __future__ import print_function

import OpenGeoSys

NODE_RELEASE_FRACTION = 0.5

class BC(OpenGeoSys.BoundaryCondition):
    def __init__(self):
        super(BC, self).__init__()
        self._restricted_nodes = {}
        self._first_node = None
        self._iteration = -1
        self._t_old = -1.0
        self._flicker_nodes = {}

    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        if self._t_old < t:
            # first iteration in current timestep
            print("restricted nodes in prev. ts", self._restricted_nodes)
            self._restricted_nodes.clear()
            self._t_old = t
            self._iteration = -1

        # detect iteration number
        # iteration 0: apply Dirichlet BCs to the solution vector prior to
        # assembly (due to Newton-Raphson solver implementation)
        # iteration 1 onwards: apply Dirichlet BCs to the linearized equation
        # system after assembly
        if self._first_node is None:
            self._first_node = node_id
        if self._first_node == node_id:
            self._iteration += 1

        x, y, z = coords

        y_deformed = y + primary_vars[1]
        y_top = BC.get_y_top(t)

        if y_deformed > y_top:
            res = (True, y_top - y)
            self._restricted_nodes[node_id] = x
            print("[BC] {it:2} {n:5} restr: y_deformed > y_top".format(
                it=self._iteration, n=node_id))
            self.restrict(node_id)

        elif y_deformed > y_top - 1e-8:
            if self.is_new_restricted_node(node_id):
                res = (True, y_top - y)
                self._restricted_nodes[node_id] = x
                print("[BC] {self._iteration:2} {node_id:5} new restr node".format(
                    it=self._iteration, n=node_id))

            elif self._iteration <= 1:
                res = (True, y_top - y)
                print("[BC] {self._iteration:2} {node_id:5} keep restr node".format(
                    it=self._iteration, n=node_id))

            elif self.is_relaxed_restricted_node(x) and not self.flickers(node_id):
                res = (False, 0.0)
                print("[BC] {self._iteration:2} {node_id:5} relax restr node".format(
                    it=self._iteration, n=node_id))
                self.relax(node_id)

            else:
                print("[BC] {self._iteration:2} {node_id:5} generally y_deformed > y_top - 1e-8".format(
                    it=self._iteration, n=node_id))
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

    def restrict(self, node_id):
        if self._iteration <= 1: return
        try:
            self._flicker_nodes[node_id][0] += 1
        except:
            self._flicker_nodes[node_id] = [1, 0]

    def relax(self, node_id):
        if self._iteration <= 1: return
        try:
            self._flicker_nodes[node_id][1] += 1
        except:
            self._flicker_nodes[node_id] = [0, 1]

    def flickers(self, node_id):
        if self._iteration <= 1: return False

        try:
            data = self._flicker_nodes[node_id]
            return data[0] >= 2 and data[1] >= 2
        except:
            return False


bc = BC()
