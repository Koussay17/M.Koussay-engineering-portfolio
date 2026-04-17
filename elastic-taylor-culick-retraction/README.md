# Elastic Taylor–Culick Retraction  
### Dynamics of Elastic Liquid Film Retraction

This repository contains a numerical and physical study of the **Taylor–Culick retraction phenomenon** in elastic or viscoelastic liquid films.

The project focuses on understanding the physics governing the rapid retraction of thin liquid films after rupture, with particular attention to the role of elasticity and surface tension.

---

## Scientific Context

The Taylor–Culick velocity describes the characteristic speed at which a thin liquid film retracts after being punctured.

This phenomenon is fundamental in:

- fluid mechanics
- multiphase flow
- thin film dynamics
- surface tension driven flows
- viscoelastic fluid physics
- instability and rupture mechanics

Understanding film retraction dynamics is essential in:

- coating processes
- spray formation
- droplet generation
- industrial fluid systems
- material science
- microfluidics

---

## Physical Problem

When a thin liquid film ruptures, surface tension forces cause the film edge to retract rapidly.

The classical Taylor–Culick velocity is given by:

V = sqrt( 2 σ / (ρ h) )

Where:

σ — surface tension  
ρ — fluid density  
h — film thickness  

This project extends the classical framework to include:

- elastic effects
- viscoelastic behavior
- dynamic film response
- non-Newtonian properties

---

## Objectives

The main objectives of this project are:

- model thin film retraction dynamics
- analyze the Taylor–Culick velocity
- investigate elastic effects on film rupture
- simulate interface motion
- study stability and energy balance
- support physical understanding through numerical simulation

---

## Core Physics

The model involves:

- surface tension forces
- inertia
- elasticity
- momentum conservation
- interface dynamics
- energy balance

---

## Repository Structure

```text
elastic-taylor-culick-retraction/

README.md

src/
    film_model.py
    retraction_solver.py
    interface_dynamics.py

examples/
    basic_retraction_simulation.py
    parameter_study.py

results/
    retraction_velocity_plot.png
    film_profile_evolution.png

docs/
    physical_model.md
    governing_equations.md
    assumptions.md
