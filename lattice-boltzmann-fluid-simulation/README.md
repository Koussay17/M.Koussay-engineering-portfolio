# Lattice Boltzmann Fluid Simulation  
### Numerical Modeling of Fluid Flow Using the Lattice Boltzmann Method (LBM)

This repository contains a numerical implementation of the **Lattice Boltzmann Method (LBM)** used to simulate fluid flow dynamics.

The project focuses on modeling fluid behavior using a mesoscopic approach based on particle distribution functions rather than directly solving the Navier–Stokes equations.

---

## Scientific Context

The Lattice Boltzmann Method is a modern computational technique widely used in:

- Computational Fluid Dynamics (CFD)
- Aerodynamics
- Multiphase flow
- Porous media flow
- Heat transfer
- Microfluidics
- Turbulence modeling

Unlike classical CFD methods, LBM models fluid motion through particle distribution functions evolving on a discrete lattice.

---

## Physical Model

The simulation is based on the discrete Boltzmann equation:

f_i(x + c_i Δt, t + Δt) = f_i(x, t) + Ω_i

Where:

- f_i — particle distribution function  
- c_i — discrete velocity direction  
- Ω_i — collision operator  
- Δt — time step  

Macroscopic fluid quantities are obtained from:

Density:

ρ = Σ f_i  

Velocity:

u = (1/ρ) Σ f_i c_i  

---

## Numerical Method

The simulation follows the standard LBM workflow:

1. Initialization  
2. Collision step  
3. Streaming step  
4. Boundary condition handling  
5. Macroscopic variable computation  

The implementation typically uses:

- D2Q9 lattice model  
- explicit time stepping  
- local collision operators  
- discrete velocity sets  

---

## Repository Structure

```text
lattice-boltzmann-fluid-simulation/

README.md

lbm.py

requirements.txt

results/
    velocity_field.png
    density_field.png
