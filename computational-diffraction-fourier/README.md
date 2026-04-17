# Tensor-Train Accelerated Diffraction Solver  
### Rigorous Fourier-Space Solution to the 1D Grating Diffraction Problem

This repository provides an implementation supporting the scientific work:

**"Logarithmically complex rigorous Fourier-space solution to the 1D grating diffraction problem"**

The numerical method implemented is the:

**Tensor-Train Accelerated Generalised Source Method (TTGSM)**

---

## Scientific Context

Rigorous diffraction modeling is a fundamental problem in computational electromagnetics and wave physics.  
Accurate simulation of diffraction by periodic structures is critical in:

- photonics
- optical engineering
- metamaterials
- wave propagation
- electromagnetic scattering
- spectral numerical methods

This repository focuses on the numerical resolution of the:

1D periodic grating diffraction problem

using a Fourier-space formulation accelerated via tensor decomposition.

---

## Method Implemented

### Tensor-Train Accelerated Generalised Source Method

The TTGSM combines:

- Generalised Source Method
- Fourier spectral formulation
- Tensor-Train compression
- logarithmic computational complexity scaling

Key advantages:

- reduced memory footprint
- improved numerical scalability
- rigorous spectral accuracy
- efficient handling of large harmonic expansions

---

## Mathematical Framework

The solver operates in Fourier space and solves the diffraction problem using:

Maxwell-type wave equations

spectral decomposition

tensor-train representation

iterative linear solvers

The method achieves:

logarithmic complexity growth

instead of quadratic or cubic scaling typical of classical solvers.

---

## Repository Structure
