#!/usr/bin/python

import numpy as np
np.set_printoptions(precision=3)

import ogs.solids as solids
import ogs.ProcessLib as PL

for k, v in solids.all_types.items():
    print(k, v)

ConstParam = PL.all_types["ConstantParameter<double>"]
MatProps   = solids.all_types["LinearElasticIsotropic<3>::MaterialProperties"]
MatModel   = solids.all_types["LinearElasticIsotropic<3>"]

E = 100e9
nu = 0.25
E_param  = ConstParam("youngs_modulus", E)
nu_param = ConstParam("poissons_ratio", nu)

props = MatProps(E_param, nu_param)
lin_el_iso3 = MatModel(props)

mat_state = lin_el_iso3.createMaterialStateVariables()


# constant params ==> actual time and position don't matter
t = 0.0
pos = PL.SpatialPosition()

dt = 1.0
eps_prev   = np.zeros(6)
eps        = np.array([0.01, 0, 0, 0, 0, 0])
sigma_prev = np.zeros(6)

sigma      = np.empty(6)
C          = np.empty((6, 6))

res = lin_el_iso3.integrateStressPythonUniqueName(
        t, pos, dt, eps_prev, eps, sigma_prev, mat_state,
        # out-variables (mat_state will be overwritten)
        sigma, mat_state, C)


assert res

print("# stress computation succeeded")
print("sigma\n", sigma)
print("tangent modulus\n", C)


# Check results

a = 1 - 2*nu
C_ref = np.matrix([
    1-nu, nu,   nu,   0, 0, 0,
    nu,   1-nu, nu,   0, 0, 0,
    nu,   nu,   1-nu, 0, 0, 0,
    0,    0,    0,    a, 0, 0,
    0,    0,    0,    0, a, 0,
    0,    0,    0,    0, 0, a
    ]).reshape(6, 6) * E / (1+nu) / a

assert (C_ref == C).all()

# make column vectors
sigma = np.atleast_2d(sigma).T
sigma_prev = np.atleast_2d(sigma_prev).T
eps = np.atleast_2d(eps).T
eps_prev = np.atleast_2d(eps_prev).T
assert (sigma - sigma_prev == C_ref * (eps - eps_prev)).all()

print("# test succeeded")  # because assertions didn't fail
