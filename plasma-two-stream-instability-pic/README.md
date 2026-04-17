# Plasma Two-Stream Instability Simulation  
### Particle-In-Cell (PIC) Numerical Model

This repository contains a numerical implementation of a **Particle-In-Cell (PIC)** simulation used to model the **two-stream instability** in plasma physics.

The simulation computes the motion of charged particles under self-consistent electric fields obtained by solving the **Poisson equation**, forming a simplified **Poisson–Maxwell plasma model**.

---

## Scientific Context

The two-stream instability is a fundamental phenomenon in plasma physics and kinetic theory.

It occurs when two particle populations move in opposite directions, leading to:

- wave growth
- energy transfer
- plasma instability
- nonlinear structure formation

This instability is widely studied in:

- plasma physics
- space physics
- fusion research
- astrophysics
- beam dynamics
- computational physics

---

## Physical Model

The simulation solves the motion of particles using:

- Newton's equation of motion
- Poisson equation for electric potential
- periodic boundary conditions
- electrostatic approximation

The electric field is computed from particle density using:

Laplacian(phi) = charge_density

The field is then interpolated back to particle positions to compute acceleration.

This numerical approach is known as the:

Particle-In-Cell method

---

## Numerical Method

The simulation uses:

- Particle-In-Cell (PIC) algorithm
- finite difference operators
- sparse matrix linear solvers
- time integration scheme (kick-drift-kick)

The implementation includes:

- particle motion update
- density interpolation
- Poisson solver
- electric field computation
- real-time visualization

---

## Repository Structure

```text
plasma-two-stream-instability-pic/

README.md

pic.py

requirements.txt

pic.png
