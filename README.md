# DBD-Plasma-CFD-Code

This repository contains the fortran code to simulate a flow through a channel (flow over a square cylinder is also present inside the code).

# CFD Equations
We solve the u, v velocities, pressure and temperature equations.

# Discretization
All the equations in time are pseudo implicit.
Convection is discretized by using the Hybrid Scheme (upwind and central differencing).
Diffusion is discretized by using the central differencing scheme.