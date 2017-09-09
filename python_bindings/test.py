#!/usr/bin/python

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
