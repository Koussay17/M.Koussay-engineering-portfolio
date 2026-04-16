# Koussay Mansouri — Aerospace & Theoretical Engineering Portfolio

Academic portfolio oriented toward aerospace engineering, theoretical physics, advanced fluid mechanics, thermal sciences, structural mechanics, and physics-based numerical simulation.

## Scientific Positioning

This repository gathers academic work developed through a strong foundation in mathematics, physics, continuum mechanics, fluid mechanics, thermal sciences, structural dynamics, and aerospace-oriented engineering analysis.

My academic profile is built around:

- Aerospace-oriented fluid mechanics
- Theoretical and applied physics
- Compressible, rotational, irrotational and viscous flows
- Thermodynamics, thermotechnics and heat transfer
- Structural mechanics and dynamics
- FEM, CFD and scientific simulation
- Instrumentation, physical measurement and systems interpretation

## Portfolio Architecture

### Theoretical Physics & Cosmology
- [Cosmological Bounce in String Quantum Gravity](./theoretical-physics-cosmology/cosmological-bounce-in-string-quantum-gravity.pdf.pdf)

Advanced theoretical work on cosmological singularity resolution, string-inspired quantum gravity, T-duality, effective string actions, and non-singular bounce cosmology.

### Aerospace Flow Physics
- [Cylinder Wake – Aerodynamic Coefficients](./aerospace-flow-physics/cylinder-wake-aerodynamic-coefficients.pdf.pdf)
- [Nozzle Study](./aerospace-flow-physics/nozzle-study.pdf.pdf)
- [Singular Pressure Losses and Flow Measurement](./advanced-fluid-thermal-sciences/singular-pressure-losses-and-flow-measurement.pdf.pdf)
- [High-Lift Devices / Flaps](./aerospace-flow-physics/high-lift-devices-flaps.pdf.pdf)

### Advanced Fluid & Thermal Sciences
- [Transition to Turbulence and Friction Factor](./advanced-fluid-thermal-sciences/transition-to-turbulence-and-friction-factor.pdf.pdf)
- [Crossflow Heat Exchangers](./advanced-fluid-thermal-sciences/crossflow-heat-exchangers.pdf.pdf)
- [Savonius Rotor – Dimensional Analysis and CFD Reproduction](./advanced-fluid-thermal-sciences/dimensional-analysis-similitude-savonius-rotor.pdf.pdf)

### Dynamics, Control & Physical Systems
- [Automatic Control – Home I/O / Simulink](./dynamics-control-physical-systems/automatic-control-homeio-simulink.pdf.pdf)
- [Interdistance Regulation](./dynamics-control-physical-systems/interdistance-regulation.pdf.pdf)

### Energetics, Fields & Measurement Systems
- [Three-Phase Measurements and Energy Assessment](./energetics-fields-measurement-systems/three-phase-measurements-energy-assessment.pdf.pdf)

### Signals, Communication & Applied Physical Investigation
- [Telecommunication ESE Project](./signals-communication-physical-investigation/telecommunication-ese-project.pdf.pdf)

### Solid Mechanics & Geometric Reasoning
- [Cone Revolution](./solid-mechanics-geometric-reasoning/cone-revolution.pdf.pdf)
- [Rolling Mobile](./solid-mechanics-geometric-reasoning/rolling-mobile.pdf.pdf)
- [Bicone Experiment](./solid-mechanics-geometric-reasoning/bicone-experiment.pdf.pdf)

### Scientific Computing & Numerical Methods
- [Polynomial Interpolation](./scientific-computing-numerical-methods/polynomial-interpolation.ipynb)
- [Iterative Methods](./scientific-computing-numerical-methods/iterative-methods.ipynb)
- [Linear Systems Resolution](./scientific-computing-numerical-methods/linear-direct-systems-methods.ipynb)

---

## ✦ Featured Theoretical Physics Work

# Cosmological Bounce in String Quantum Gravity

> Research-oriented theoretical work on non-singular early-universe cosmology, string-inspired quantum gravity, T-duality, effective string actions, and high-curvature dynamics.

This paper studies a central problem in theoretical cosmology: the classical Big Bang picture leads to an initial singularity, meaning that curvature and density diverge and the classical description breaks down.  
The objective of this work is to investigate whether this singular origin can be replaced by a **non-singular cosmological bounce** in a **string-inspired quantum-gravity framework**, where the Universe passes from contraction to expansion through a regulated high-curvature phase rather than an infinite breakdown. 

The work is structured around the following theoretical ingredients:
- effective string cosmology,
- T-duality,
- dilaton dynamics,
- α′ high-curvature corrections,
- modified cosmological equations,
- proper-time deformation through the bounce. 

### ✦ What is actually studied

This work does **not** merely restate general cosmology.  
It focuses on a precise theoretical question:

- how to encode string-inspired corrections into cosmological dynamics,
- how the dilaton modifies the evolution,
- how the corrected equations avoid a singular behavior,
- and how the notion of proper time changes across the bounce regime.

In other words, the paper studies a cosmological model in which the earliest phase of the Universe is governed not only by classical expansion laws, but by a corrected effective dynamics capable of regularizing the transition itself. 
---

## ✦ Core Equations of the Work

### 1. Effective string action with high-curvature corrections

$$
S_{\mathrm{eff}}=
\frac{1}{2\kappa^2}
\int d^4x\,\sqrt{-g}\,e^{-2\Phi}
\left[
R
+4(\nabla\Phi)^2
-\frac{1}{12}H^2
+\frac{\alpha'}{4}\mathcal{R}^2
+\mathcal{O}(\alpha'^2)
\right]
$$

This is the central theoretical starting point of the work.  
It describes an effective cosmological model where curvature, the dilaton field, and string-scale corrections all contribute to the dynamics. The crucial point is that the \(\alpha'\) sector introduces high-curvature corrections that become relevant near the bounce. 

---

### 2. Cosmological background used in the model

$$
ds^2 = -dt^2 + a(t)^2\,d\vec{x}^{\,2}
$$

The model is then specialized to a flat homogeneous and isotropic cosmological background.  
The function \(a(t)\) is the scale factor, and it controls the global contraction/expansion behavior of the Universe. 
---

### 3. Variables used to rewrite the dynamics

$$
\lambda(t)=\ln a(t),
\qquad
H=\dot\lambda=\frac{\dot a}{a},
\qquad
\bar{\Phi}(t)=\Phi(t)-3\lambda(t)
$$

These variables are introduced to express the cosmological system in a compact and physically interpretable way.  
The Hubble parameter \(H\) tracks the expansion rate, while the shifted dilaton \(\bar{\Phi}\) packages the coupling between geometry and dilaton evolution.
---

### 4. Corrected cosmological equations

$$
-\dot{\bar{\Phi}}^{\,2}+6H^2+\alpha'c_1H^4+\alpha'c_2H^2\dot{\bar{\Phi}}^{\,2}=0
$$

$$
\ddot{\bar{\Phi}}+\dot{\bar{\Phi}}^{\,2}-3H^2+\alpha'c_3H^2\left(3H^2+\dot{\bar{\Phi}}^{\,2}\right)=0
$$

These are the main corrected evolution equations of the work.  
The first is the effective constraint equation.  
The second governs the evolution of the shifted dilaton.  
The key feature is that the \(\alpha'\)-dependent nonlinear terms modify the classical high-curvature behavior and open the door to a non-singular bounce.
---

### 5. Tree-level limit of the system

$$
\dot{\bar{\Phi}}^{\,2}=6H^2
$$

$$
\ddot{\bar{\Phi}}+\dot{\bar{\Phi}}^{\,2}=3H^2
$$

These equations correspond to the lower-order regime before the full high-curvature corrections dominate.  
They help show what the uncorrected string cosmology would look like, and therefore what is genuinely changed by the corrected dynamics.

---

### 6. Proper-time dilation relation

$$
d\tau = e^{-\Phi/2}\,dt
$$

This is one of the most interesting physical relations in the work.  
It expresses the difference between the cosmological time parameter \(t\) and the Einstein-frame proper time \(\tau\).  
Near the bounce, the evolution of the dilaton modifies the rate at which physical time flows, which gives a physically meaningful interpretation of the bounce beyond geometry alone.
---

### 7. Compact form of the bounce constraint

$$
-\dot\Phi^2 + 6H\dot\Phi - 3H^2
+\frac{\alpha'}{2}\left(a_1 H^4 + a_2 H^2 \dot\Phi^2\right)=0
$$

This compact form makes the competition between expansion, dilaton evolution, and high-curvature corrections especially clear.  
It condenses the structure of the bounce mechanism into a single expression: the corrected dynamics are no longer purely Einsteinian, and the singular regime can be regularized.
---

## ✦ Physical Interpretation

The logic of the paper is not to add corrections arbitrarily, but to show that a coherent string-inspired effective theory naturally modifies the cosmological equations at high curvature.  
The physical picture is the following:

- the classical singular regime is not accepted as final,
- the effective dynamics are corrected by string-scale contributions,
- the dilaton becomes dynamically relevant,
- time itself is no longer trivially identified across frames,
- and the contraction-to-expansion transition becomes mathematically describable without infinite breakdown. 

This makes the work a genuinely research-oriented study in:
- theoretical cosmology,
- gravitational physics,
- effective string dynamics,
- and mathematically grounded early-universe modeling. 

---

## ✦ Why this work stands out in the portfolio

This paper marks a clear shift from advanced engineering coursework toward **theoretical physics research-style reasoning**.  
It combines:
- abstract theoretical structure,
- action-based modeling,
- nontrivial cosmological interpretation,
- string-inspired corrections,
- and a physically meaningful treatment of the singularity problem.
[Read the paper](./theoretical-physics-cosmology/cosmological-bounce-in-string-quantum-gravity.pdf.pdf)

---

## Main Theoretical Axes

- Continuum mechanics
- Solid and structural mechanics
- Ideal and viscous fluid mechanics
- Compressible-flow reasoning
- Boundary layer and turbulence
- Thermodynamics and heat transfer
- Turbomachinery and propulsion-oriented analysis
- Finite Element Method (FEM)
- Computational Fluid Dynamics (CFD)
- Measurement, modeling and physical interpretation

## Academic Intent

This portfolio is not limited to software usage or laboratory execution.
It reflects a physics-grounded engineering mindset focused on understanding the mathematical and physical structure of complex systems, especially in aerospace and multiphysics environments.

## Author

**Koussay Mansouri**
