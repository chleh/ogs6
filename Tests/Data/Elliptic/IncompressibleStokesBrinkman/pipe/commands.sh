#!/bin/sh

OGSPATH=/home/lehmannc/prog/ogs6/github-chleh-PRs/build-release-eigenlis/bin

"$OGSPATH/generateStructuredMesh" \
    -e quad --lx 1 --ly 1 --nx 20 --ny 20 -o pipe.vtu

"$OGSPATH/createQuadraticMesh" -i pipe.vtu -o pipe_quadratic.vtu

false && /home/lehmannc/prog/ogs-utils/pre/BCs/inhomogeneous-mass-flux/prepare-inhomogeneous-mass-flux-BC.py \
    --input pipe_quadratic.vtu \
    --output pipe_inlet.vtu \
