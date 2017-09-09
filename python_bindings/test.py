#!/usr/bin/python

import numpy as np

import ogs.solids as solids
import ogs.ProcessLib as PL

for k, v in solids.all_types.items():
    print(k, v)

ConstParam = PL.all_types["ConstantParameter<double>"]
MatProps   = solids.all_types["LinearElasticIsotropic<3>::MaterialProperties"]
MatModel   = solids.all_types["LinearElasticIsotropic<3>"]

E  = ConstParam("youngs_modulus", 100e9)
nu = ConstParam("poissons_ratio", 0.25)

props = MatProps(E, nu)
lin_el_iso3 = MatModel(props)

mat_state = lin_el_iso3.createMaterialStateVariables()

help(lin_el_iso3)

# constant params ==> actual time and position don't matter
t = 0.0
pos = PL.SpatialPosition()
dt = 1.0
eps_prev = np.array([ 0 ] * 6)
eps = np.array([ 0 ] * 6)
sigma_prev = np.array([ 0 ] * 6)

res = lin_el_iso3.integrateStressPythonUniqueName(t, pos, dt, eps_prev, eps, sigma_prev, mat_state)

print(res)
