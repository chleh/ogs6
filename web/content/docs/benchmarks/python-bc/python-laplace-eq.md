+++
project = "https://github.com/ufz/ogs-data/blob/master/Elliptic/square_1x1_GroundWaterFlow_Python/square_1e3_laplace_eq.prj"
author = "Christoph Lehmann"
date = "2018-06-01T14:16:55+02:00"
title = "Manufactured Solution for Laplace's Equation with Python"
weight = 156

[menu]
  [menu.benchmarks]
    parent = "python-bc"

+++

{{< data-link >}}

## Problem description

We solve Laplace's Equation in 2D on a $1 \times 1$ square domain.
The objective of this test is to check whether both essential and natural
boundary conditions are implemented correctly in OpenGeoSys's Python boundary
conditions. The boundary conditions will be specified below.

## Weak form

Laplace's equation is
$$
\[
- \mathop{div} (a \mathop{grad} u) = 0
\]
$$
The weak form is derived as usual by multiplying with a test function $v$ and
integrating over the domain $\Omega$:
$$
\[
- \int_\Omega v \mathop{div} (a \mathop{grad} u) \, \mathrm d\Omega
\,,
\]
$$
which can be transformed further to
$$
\[
\int_\Omega a \mathop{grad} v \cdot \mathop{grad} u \, \mathrm d\Omega
= \int_\Omega \mathop{div} (v a \mathop{grad} u) \, \mathrm d\Omega
= \int_(\Gamma_\mathrm{N}) v a \mathop{grad} u \cdot n \, \mathrm d\Gamma
\,,
\]
$$
where in the second equality Gauss's theorem has been applied.
As usual, the domain boundary $\partial\Omega = \Gamma_\mathrm{D} \cup \Gamma_\mathrm{N}$ is subdivided
into the dirichlet and the Neumann boundary and $v$ vanishes on
$\Gamma_\mathrm{D}$.
The r.h.s. of the above equation is the total flux associated with $u$ flowing
**into** the domain $\Omega$ through $\Gamma_\mathrm{N}$:
$-a \grad u$ is the flux density and $-n$ is the inwards directed surface
normal.

The weak form just derived is implemented (after FEM discretization) in  the
groundwater flow process in OpenGeoSys.
Note that for the application of Neumann boundary conditions, it is necessary to
know whether the flux has to be applied with a positive or a negative sign!

## Analytial solution

The coefficient $a$ of Laplace's equation is taken to be unity.
By differentiation it can be easily checked that
$$
\[
u(x, y) = \sin(bx) \sinh(by)
\]
$$
solves Laplace's equation inside $\Omega$ for any $b$.
In this example we set $b = \tfrac 23 \pi$.

As boundary conditions we apply Dirichlet BCs at the top, left and bottom of the
domain with values from $u(x,y)|_{\Gamma_\mathrm{D}}$.
On the right boundary of the domain a Neumann BC is applied.
There $n = (1, 0)$, which implies that $a \mathop{grad} u \cdot n
= a \partial u / \partial x$.


## Results

The numerical result obtained from OpenGeoSys is:

{{< img src="../python_laplace_eq_solution.png" >}}

The absolute difference between the analytical and numerical solutions is
smaller than $4 \cdot 10^{-4}$:

{{< img src="../python_laplace_eq_diff.png" >}}
