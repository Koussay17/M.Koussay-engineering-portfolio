<div align="center">

# Koussay Mansouri

### Theoretical Physics · Continuum Mechanics · Scientific Computing

*A research portfolio built on one rule: if it can be proved, prove it.
If it can be measured, measure it. If it fails, say where.*

<br>

![Papers](https://img.shields.io/badge/research_papers-10-1f3b73?style=for-the-badge)
![Pages](https://img.shields.io/badge/pages-138-2f7d4f?style=for-the-badge)
![Theorems](https://img.shields.io/badge/theorems_proved-40+-b3452c?style=for-the-badge)
![Language](https://img.shields.io/badge/papers-français-7a5aa8?style=for-the-badge)

![Python](https://img.shields.io/badge/Python-numpy_·_scipy-3776AB?style=flat-square&logo=python&logoColor=white)
![LaTeX](https://img.shields.io/badge/LaTeX-source_included-008080?style=flat-square&logo=latex&logoColor=white)
![Reproducible](https://img.shields.io/badge/reproducible-every_figure-success?style=flat-square)

</div>

---

## What this repository is

Ten self-contained research papers in theoretical physics, each pairing an
analytical derivation with an **independent numerical verification**. Every
paper follows the same discipline:

| | |
|---|---|
| **Prove what is provable** | theorems with complete proofs, not citations |
| **Verify, don't illustrate** | a simulation that *shows* an effect proves nothing — the growth rate must match the exact dispersion relation |
| **Validate without external references** | energy conservation, constraint propagation, jump conditions — quantities that test the code against *itself* |
| **State the limits** | every paper ends with what does **not** work, and why |

The papers are written in French. Every number quoted below is produced by the
accompanying Python script and is reproducible from a clean checkout.

---

## The papers

| # | Paper | Field | Headline result | pp. |
|:-:|---|---|---|:-:|
| 01 | [Rebond cosmologique en gravité des cordes](theoretical-physics-cosmology/) | String cosmology | Duality threshold: **no bounce below the string scale** | 23 |
| 01b | [Carnet de calculs](theoretical-physics-cosmology/) | Companion | Every derivation of #01, from scratch | 31 |
| 02 | [Instabilité double faisceau](plasma-two-stream-instability-pic/) | Kinetic plasma theory | PIC matches Vlasov theory to **0.70 %** | 17 |
| 03 | [Mécanique non holonome du roulement](solid-mechanics-geometric-reasoning/) | Geometric mechanics | Lie bracket **predicts** the measured holonomy | 14 |
| 04 | [De Boltzmann à Navier–Stokes](lattice-boltzmann-fluid-simulation/) | Kinetic theory | The $-\tfrac12$ is worth **271 %** | 10 |
| 05 | [Le piège de Fourier des réseaux](computational-diffraction-fourier/) | Periodic EM | Li's rules: $N^{+0.07}\!\to\!N^{-1.15}$ | 9 |
| 06 | [Instabilité de cisaillement](Kelvin-Helmholtz-Instability/) | Hydrodynamic stability | Exact neutral mode, residual $5.6\!\times\!10^{-17}$ | 9 |
| 07 | [Intégration symplectique et chaos](celestial-dynamics-model%20multi-body-gravity-simulation/) | Dynamical systems | Order-2 **beats** order-4 after 700 periods | 9 |
| 08 | [Rétraction viscoélastique](elastic-taylor-culick-retraction/) | Polymer rheology | Oldroyd-B is a **category error** here | 9 |
| 09 | [Turbulence : ce qui est démontré](advanced-fluid-thermal-sciences/) | Turbulence | $\zeta_3 = 0.98397$ against exactly 1 | 8 |
| 10 | [Lois de conservation hyperboliques](aerospace-flow-physics/) | Gas dynamics | Weak solutions are **not unique** — exhibited | 9 |

<br>

<details>
<summary><b>01 — Rebond cosmologique non singulier en gravité des cordes</b> · <i>seuil de dualité, théorèmes d'obstruction, solution globale</i></summary>

<br>

**The question.** Can the Big Bang singularity be replaced by a bounce, using
the $\mathrm{O}(d,d)$-invariant $\alpha'$ corrections of string theory?

**What is proved.** Working in the Hohm–Zwiebach framework, where all
$\alpha'$ corrections collapse into a single even function $F(H)$:

- **Consistency** — the conservation law $\dot{\bar\rho}+dH\bar p=0$ holds for
  *any* resummation $F$.
- **No string-frame vacuum bounce** — if $H(t_0)=0$ in vacuum then $H\equiv 0$.
- **The threshold theorem** — an Einstein-frame bounce occurs iff
  $\dot{\bar\Phi}=-H$, which forces $W(H)\equiv HF'-F \le H^2$. This is
  *identically violated* at leading order, where $W=dH^2$. **A bounce cannot
  happen below the string scale.**
- **Graceful-exit obstruction** — while $F\ge0$, $\dot{\bar\Phi}$ is monotone,
  so the exit to a decelerating FRW branch is impossible without either
  $W<0$ or negative potential energy.
- **Branch singularity** — the system blows up at any simple zero of $F''$.

**The solution.** With $F''=6/(1+H^2/H_*^2)$, an $\mathrm{O}(d,d)$ string gas
and a duality-invariant potential:

```
threshold          H_thr    = 2.389862 H*     (root of 3 ln(1+x) = x)
bounce             |H_b|    = 1.0134 × H_thr
constraint         |C|max   = 3.9 × 10⁻¹²     (never imposed — diagnostic only)
coupling           g_s      ≤ 8.4 × 10⁻³      (perturbative throughout)
past asymptotics   a_E ∝ |τ|^0.4990           radiation-dominated contraction
future asymptotics a_E ∝  τ ^0.4997           radiation-dominated expansion
α′ share of the NEC-violation budget          99.83 %
```

Radiation → bounce → radiation, geodesically complete both ways.
The exotic energy needed at tree level ($V=11.84$) drops to $0.0203$: the
bounce is a **string effect**, quantifiably.

**What is not done.** $F$ is postulated, not derived from amplitudes.
Perturbations are not computed. A radiation-dominated contraction is unstable
to anisotropy — quantified in the companion notebook: the initial shear must
satisfy $(\sigma^2/\kappa^2\rho)_{\rm init}\lesssim10^{-9}$.

**Companion: `Carnet de calculs`** — 31 pages rebuilding everything from
manifolds and connections. Levi-Civita existence *and* uniqueness, Palatini,
the lapse as a Lagrange multiplier, the full conformal transformation law
derived by expanding nine terms, Gauss–Hermite, the exact
$\mathrm{tr}(\dot{\mathcal S}^{2k})=(-1)^k 2d\,4^k H^{2k}$, and the exact
propagation law $\mathcal C(t)=\mathcal C(t_0)e^{\bar\Phi(t)-\bar\Phi(t_0)}$
that *predicts* the numerical error growth. Plus eight worked exercises.

</details>

<details>
<summary><b>02 — L'instabilité double faisceau dans Vlasov–Poisson</b> · <i>réponse de Landau, critère de Penrose, saturation par piégeage</i></summary>

<br>

**The question.** How does a collisionless plasma go unstable, and can a
particle code reproduce kinetic theory quantitatively?

**What is proved.**

- The dielectric function with the **Landau contour** — including *why*
  principal value is wrong and why it is Laplace, not Fourier, in time.
- **Penrose's criterion** by the Nyquist argument, in the strong form that also
  gives the marginal wavenumber $k_c=\omega_p\sqrt{\mathcal P}$ with
  $\mathcal P=\int[f_0(v)-f_0(v_m)](v-v_m)^{-2}\dd v$.
- **Cold beams solved exactly** — the dispersion relation reduces to a
  biquadratic: unstable for $kv_0<\omega_p$, with
  $\gamma_{\max}=\omega_p/(2\sqrt2)$ at $kv_0/\omega_p=\sqrt{3/8}$.
- **New**: $\varepsilon(k,i\gamma)$ is *real* for a symmetric distribution, so
  the warm problem reduces to a real bisection — no complex root finding.

```
cold exact           γ_max = ω_p/(2√2) = 0.353553   at kv₀ = √(3/8) = 0.612372
Penrose threshold    v₀/v_th = 1.277307   (below it, NO wavenumber is unstable)
cross-check          Penrose k_c vs Faddeeva root agree to 0.3 %
cold limit of Penrose → k_c v₀ = ω_p   exactly the algebraic threshold
PIC vs kinetic theory                    0.70 % mean over 6 modes
saturation           ω_b/γ = 1.35        trapping estimate confirmed
energy conservation  2.5 × 10⁻³
```

**Method note.** The PIC code is written from scratch with a **tensor quiet
start** ($2\times2048\times256$ particles, exactly uniform deposited density,
zero shot noise), giving eight growth times of clean linear phase. The
measurement protocol declares three modes **unmeasurable** rather than quoting
a contaminated number.

**What is not done.** 1D–1V has no analogue of Squire's theorem, so the
reduction is a modelling assumption, not a controlled simplification —
contrast with paper 06, where it *is* proved.

</details>

<details>
<summary><b>03 — Mécanique non holonome du roulement</b> · <i>non-intégrabilité, holonomie, trois systèmes exacts</i></summary>

<br>

**The question.** Rolling without slipping constrains velocities but not
configurations. What is the geometry of that, and can the defect be *measured*?

**What is proved.**

- Explicit Lie brackets for the rolling sphere: $[X_1,X_2]=a^{-2}E_z\notin\mathcal D$,
  so by **Frobenius** the constraint is non-integrable. Then $[X_1,E_z]$ and
  $[X_2,E_z]$ generate the rest — bracket-generating, hence **completely
  controllable** by Chow–Rashevskii (the "parallel parking theorem").
- Lagrange–d'Alembert equations are **not variational**: the vakonomic
  equations differ, and coincide only when the constraint is integrable.
- The rolling cone with fixed apex is **holonomic** — its constraint integrates
  to $\psi\sin\alpha+\varphi=\mathrm{const}$. *Rolling does not imply
  non-holonomic.*

**The measurement.** The bracket *is* the curvature of the rolling connection,
so it predicts $\Theta\simeq\mathcal A/a^2$. Exact transport by ordered
products of $\mathfrak{so}(3)$ exponentials:

```
𝒜/a² = 0.0025   →   Θ/𝒜 = 0.9997918      the area law, confirmed
second order:  Θ = 𝒜(1 − 𝒜/12 + …)  for a square
shape dependence at equal area 𝒜 = 0.36:
    square 0.6×0.6      Θ = 0.349774
    rectangle 2.4×0.15  Θ = 0.279577      ← non-abelian correction
reversal check   ‖R_direct · R_inverse − I‖ = 7 × 10⁻¹³
```

First order is universal (the curvature); second order depends on loop shape —
the signature of a **non-abelian** structure group. For an abelian connection
(Aharonov–Bohm) $\Theta$ would be exactly proportional to area.

**Three systems, three lessons.**

| system | result |
|---|---|
| rolling cone | $I_c=\tfrac{3}{20}mh^2\sin^2\alpha(6+\tan^2\alpha)$, uniform precession, **holonomic** |
| double cone paradox | $\tan\theta<\tan\beta\tan\alpha$ — no mass, no inertia, no $g$ |
| Chaplygin sleigh | **an attractor with energy conserved to $1.8\times10^{-12}$** |

The sleigh is the gem: $\dot v=a\omega^2$, $\dot\omega=-\frac{ma}{J}v\omega$,
integrated in closed form via $\dot\phi=-\lambda\sin\phi$. Numerics agree with
the exact solution to $1.9\times10^{-12}$, $\omega(40)=10^{-14}$,
$v/v_\infty=1.0000000000$. Divergence of the vector field is $-\frac{ma}{J}v\ne0$:
non-holonomic systems are neither Hamiltonian nor dissipative — a third class.

</details>

<details>
<summary><b>04 — De Boltzmann à Navier–Stokes</b> · <i>Chapman–Enskog, Gauss–Hermite, la viscosité au pour-mille</i></summary>

<br>

**The question.** Why is Navier–Stokes universal, when it contains almost none
of the molecular physics?

**What is proved.**

- Chapman–Enskog in full, with the **solvability conditions** that actually
  produce the hydrodynamic equations: order $\mathrm{Kn}^0\to$ Euler,
  order $\mathrm{Kn}^1\to$ Navier–Stokes with $\nu=c_s^2\tau$.
- Why one stops: Burnett (order $\mathrm{Kn}^2$) is linearly **unstable**
  (Bobylev). Going further in an asymptotic expansion does not help.
- The D2Q9 weights are **not arbitrary** — they are 3-point Gauss–Hermite
  quadrature. $\tfrac23\!\cdot\!\tfrac23=\tfrac49$,
  $\tfrac16\!\cdot\!\tfrac23=\tfrac19$, $\tfrac16\!\cdot\!\tfrac16=\tfrac1{36}$,
  and $c_s^2=1/3$ is the abscissa rescaling.
- The 4th-order isotropy condition
  $\sum_iw_ic_{i\alpha}c_{i\beta}c_{i\gamma}c_{i\delta}=c_s^4(\delta\delta+\text{perm})$
  is *exactly* what makes the stress tensor isotropic.
- The $-\tfrac12$ comes from **trapezoidal integration** of the collision term
  plus the change of variable $\bar f_i=f_i+\frac{\Delta t}{2\tau}(f_i-f^{\rm eq}_i)$.

```
lattice checks    Σw = 1.000000000000000     Σw c_x² = 0.333333333333333
isotropy          Σw c_x⁴ = 3 Σw c_x²c_y²   (0.333333 vs 0.333333)   exact

ν measured vs c_s²(τ−½)                  mean deviation  0.156 %
ν measured vs naive c_s²τ                mean deviation    271 %
  at τ = 0.55                            naive is wrong by a factor 10
Taylor–Green (independent geometry)      0.002 % – 0.32 %
spatial convergence   error ratios       4.04 · 4.04 · 3.98   → order 2
compressibility error                    ∝ Ma^1.966
```

**Method note kept in the paper.** My first compressibility test was *badly
designed*: for the shear wave $(\mathbf u\!\cdot\!\nabla)\mathbf u=0$
identically, so it is an exactly linear solution — perfect for measuring $\nu$,
blind to compressibility. Hence the switch to Taylor–Green. Saying this is
worth more than a smooth table.

</details>

<details>
<summary><b>05 — Le piège de Fourier des réseaux de diffraction</b> · <i>règles de factorisation de Li</i></summary>

<br>

**The question.** Why did the Fourier modal method converge beautifully in TE
and catastrophically in TM for twenty years?

**The answer is not physics and not numerics.** It is Fourier analysis: the
Fourier series of a *product* is not always the truncated convolution.

**What is proved.** Li's three rules, and the electromagnetic diagnosis that
follows in three lines from the interface conditions ($E_y,E_z$ continuous;
$D_x=\varepsilon E_x$ continuous, hence $E_x$ **discontinuous**):

- $D_y=\varepsilon E_y$ — no concurrent discontinuity → **Laurent**.
  *This is why TE never had a problem.*
- $D_x=\varepsilon E_x$ — concurrent, complementary → **inverse rule**
  $[\![1/\varepsilon]\!]^{-1}$.

The subtlety that makes the error so easy: in the TM equation both products
involve the **same** factor $1/\varepsilon$ and require **different** rules —
because $\partial_xH_y\propto\varepsilon E_z$ is discontinuous while
$\partial_zH_y\propto D_x$ is continuous.

```
validation (lossless metal → energy exactly conserved)
    R + T = 1.000000000000        in TM and TE, at every truncation

convergence of DE₀ (TM, metallic lamellar grating)
    N = 32    naive 1.33e-3    Li 1.58e-4    gain  8.4
    N = 48    naive 1.77e-2    Li 3.11e-4    gain 56.9   ← naive error GREW ×13
    N = 64    naive 8.20e-3    Li 1.28e-4    gain 64.2

fitted order      naive N^(+0.07)     Li N^(−1.15)     TE N^(−2.86)
same limit?       0.301692 ± 2.9e-5   vs   0.301969 ± 4.9e-4      ✓
```

**Honest ceiling.** Even corrected, convergence caps at $N^{-1.15}$, not
$N^{-3}$ — that is the **Meixner edge singularity** at the four metal corners,
$|\mathbf E|\sim\rho^{\sigma-1}$, which bounds Fourier coefficients regardless
of factorization. The diagnosis names the right remedy: not a better
convolution, but adaptive spatial resolution near the edges.

</details>

<details>
<summary><b>06 — Instabilité de cisaillement non visqueuse</b> · <i>quatre théorèmes, la couche en tanh résolue</i></summary>

<br>

**The question.** Shear instability needs no viscosity, no stratification, no
rotation. What exactly can be proved about it?

**Four theorems, one identity.** Everything follows from multiplying the
Rayleigh equation by $\varphi^*/(U-c)$ and integrating:

$$\int\big(|\varphi'|^2+k^2|\varphi|^2\big)\dd y+\int\frac{U''}{U-c}|\varphi|^2\dd y=0$$

| theorem | content |
|---|---|
| **Squire** | 2D disturbances dominate — the reduction is *proved*, not assumed |
| **Rayleigh** | imaginary part ⟹ $U''$ must change sign |
| **Fjørtoft** | real part + Rayleigh's null term ⟹ vorticity max must be interior |
| **Howard** | $F=\varphi/(U-c)$ and $(U-U_{\min})(U-U_{\max})\le0$ ⟹ semicircle |

**The exact solution.** For $U=\tanh y$: $\varphi=\operatorname{sech} y$ at
$k=1$, $c=0$ solves the Rayleigh equation exactly (three lines, via
$\tanh^2=1-\operatorname{sech}^2$).

```
analytic residual of the exact neutral mode        5.6 × 10⁻¹⁷
numerical kc_i at k = 1                            3.56 × 10⁻⁵    (→ 0)
maximum growth   kδ = 0.444939   kc_i = 0.189701
reference values        0.4446           0.1897     → 0.08 % / 0.00 %
⟹ selected wavelength   λ ≈ 14.1 δ    (parameter-free prediction)

Howard's semicircle   max |c−c₀|/r = 0.827 ≤ 1     no violations
Rayleigh test  U = 1−y² (no inflection)  →  max c_i = 0 exactly, all k
vortex sheet   γ = kΔU/2  saturates Howard's bound but diverges with k
regularization  k_max δ = 0.4477 invariant over a decade of δ
```

**Limit worth stating.** The Rayleigh equation predicts plane Poiseuille is
stable. It is not — it transitions at $\mathrm{Re}\simeq5772$. Viscosity can
**destabilize** (Tollmien–Schlichting), by a viscous critical-layer mechanism
entirely absent here.

</details>

<details>
<summary><b>07 — Intégration symplectique et chaos hamiltonien</b> · <i>hamiltonien modifié, KAM, Lyapunov</i></summary>

<br>

**The question.** Why does nobody integrate the Solar System with Runge–Kutta?

**What is proved.** By symmetric BCH, leapfrog integrates *exactly* a nearby
Hamiltonian:

$$\widetilde H=H+\Delta t^2\Big[-\tfrac1{24}|\nabla V|^2+\tfrac1{12}p_ip_j\partial_i\partial_jV\Big]+O(\Delta t^4)$$

Brackets computed in full: $\{\{T,V\},V\}=|\nabla V|^2$,
$\{\{V,T\},T\}=p_ip_j\partial_i\partial_jV$. **Coefficients cross-checked by an
entirely independent route** — leapfrog on the harmonic oscillator is a
$2\times2$ symplectic matrix whose exact quadratic invariant is
$\tfrac12p^2+\tfrac12\omega^2(1-\Delta t^2\omega^2/4)q^2$, matching BCH after
normalization.

```
Kepler, 2500 orbits, Δt = 0.05
    leapfrog  |ΔH/H₀|max = 6.84e-3    envelope growth  t^(−0.000)   BOUNDED
    RK4       |ΔH/H₀|max = 2.38e-2    envelope growth  t^(+0.942)   SECULAR
    crossover at ~700 periods — order 2 becomes 3× more accurate than order 4

the modified Hamiltonian, measured
    Δt      max|ΔH|     order      max|ΔH̃|     order      gain
    0.020   5.44e-4      2.00      4.91e-7      4.01      1 107
    0.005   3.40e-5      2.00      1.92e-9      4.00     17 726

angular momentum (central force ⟹ ∇V ∥ q)        8.6 × 10⁻¹⁴   exactly
```

If $\widetilde H$ were merely formal, nothing would force its variation to be
*two orders* higher. It is, exactly.

**Chaos.** Hénon–Heiles Poincaré sections and Benettin exponents:

| $E$ | orbit 1 | orbit 2 | orbit 3 |
|---|---|---|---|
| $E_{\rm esc}/2$ | $-0.00000$ | $+0.00003$ | $+0.00051$ |
| $3E_{\rm esc}/4$ | $-0.00000$ | $+0.00002$ | $\mathbf{+0.03495}$ |
| $E_{\rm esc}$ | $+0.06759$ | $+0.10835$ | $+0.13647$ |

**There is no critical energy.** At the same energy some orbits sit on tori and
others in the stochastic sea — phase space is *mixed*, and the chaotic fraction
grows continuously. Exactly what KAM predicts.

**Cross-reference.** A Hamiltonian Lyapunov spectrum is symmetric
($\pm\lambda_i$, sum zero), so a Hamiltonian system cannot have an attractor —
compare the Chaplygin sleigh in paper 03, which *does*, precisely because it is
not Hamiltonian.

</details>

<details>
<summary><b>08 — Rétraction capillaire d'un film viscoélastique</b> · <i>Taylor–Culick, la singularité à Wi = 1/2</i></summary>

<br>

**The question.** A soap film retracts at constant speed. Add polymer and it
slows by a factor of five. What breaks?

**Taylor–Culick, and the mistake everyone makes.** $V=\sqrt{2\sigma/\rho h}$ is
a **momentum** balance. The energy balance gives $\sqrt{4\sigma/\rho h}$ —
wrong by exactly $\sqrt2$, measured ratio $1.414214$. Half the surface energy
is dissipated at the film/rim junction **even in a perfect fluid**: it is a
perfectly inelastic collision, and viscosity sets the thickness of the joint,
not the amount dissipated. The transient integrates to $V$ constant *from
$t=0$* (0.0000 % deviation).

**The exact singularity.** In uniaxial extension, Oldroyd-B has

$$A_{xx}=\frac{1}{1-2\mathrm{Wi}},\qquad \eta_E=3\eta_s+\frac{3\eta_p}{(1-2\mathrm{Wi})(1+\mathrm{Wi})}$$

which **diverges at $\mathrm{Wi}=1/2$**. Beyond it, exponential growth at rate
exactly $(2\mathrm{Wi}-1)/\lambda$.

```
steady state       A_xx at Wi = 0.45     10.000000  vs  10.000000 exact
extensional visc.  η_E at Wi = 0.48      45.908103  vs  45.908108 exact
blow-up rate       Wi = 1.0               1.00000000 vs  1.00000000 exact
FENE-P saturation  A_xx → L²(1−1/2Wi)     499.833   vs  500 predicted
```

**The real numerical failure mode.** The conformation tensor is a covariance
matrix — it *must* stay positive definite. Explicit Euler, relaxing from
$A_0=50$:

| $\Delta t/\lambda$ | direct | log-conformation |
|---|---|---|
| 0.90 | $5.900$ | $20.698$ |
| 1.05 | $\mathbf{-1.450}$ | $17.868$ |
| 2.00 | $\mathbf{-48.000}$ | $7.043$ |

Negative **in a single step** once $\Delta t/\lambda>1$. This is the
high-Weissenberg-number problem that blocked computational rheology for twenty
years. Fattal–Kupferman fixes it *by construction* — but is **not more
accurate**, a nuance the paper states.

**The conclusion, and it is negative.** The rim strain rate is
$\dot\varepsilon\sim V/h$, so the threshold is crossed once
$\lambda>\lambda_c=h/(2V)\propto h^{3/2}$. For a $10\,\mu$m water film,
$\lambda_c=1.3\,\mu$s. Real relaxation times are **three to seven orders of
magnitude larger**. The elastocapillary regime of interest lies *entirely*
beyond the divergence. Using Oldroyd-B there is not an approximation — it is a
category error.

</details>

<details>
<summary><b>09 — Turbulence : ce qui est démontré et ce qui ne l'est pas</b></summary>

<br>

**Three levels, stated separately.**

**Proved.** Buckingham's theorem (by a rank argument: $\dim\ker D=n-r$), which
gives K41 in one line. Then the only exact non-trivial result in the field —
Kármán–Howarth–Monin → the $4/3$ law → **Kolmogorov's $4/5$ law**:

$$\langle\delta u_L^3\rangle=-\tfrac45\,\varepsilon\,r\qquad\Longrightarrow\qquad \zeta_3=1\ \text{exactly}$$

No similarity hypothesis — only Navier–Stokes, homogeneity, isotropy,
stationarity. The **negative sign** is the physical content: it encodes the
direction of the cascade.

**Assumed.** The dissipative anomaly: $\varepsilon$ finite as $\nu\to0$,
requiring $\langle|\nabla\mathbf u|^2\rangle\to\infty$ like $1/\nu$. Developed
turbulence is a **singular** object in the inviscid limit.

**Wrong in detail.** K41 for $p\ne3$, by Landau's objection:
$\langle\varepsilon^{p/3}\rangle\ne\langle\varepsilon\rangle^{p/3}$. The real
exponents are **anomalous dimensions** — not by analogy, in the exact sense of
field theory. Nobody can compute them.

```
validation first    nonlinear term conserves energy    6.3 × 10⁻¹³
                    and generalized helicity           1.8 × 10⁻¹³

diagnosed artifact  GOY period-3 oscillation, 35 % modulation
                    uncorrected  ζ₃ = 0.896   (10 % error on an EXACT exponent)
                    corrected    ζ₃ = 0.98397 (1.60 %)  ← triad geometric mean

spectrum            E(k) ∝ k^(−1.702)     not −5/3, and that is CORRECT:
                    ζ₂ > 2/3, and experiment also gives ≈ −1.70
```

| $p$ | $\zeta_p/\zeta_3$ | K41 | She–Levêque |
|---|---|---|---|
| 3 | 1.0000 | 1.0000 | 1.0000 |
| 6 | **1.7240** | 2.0000 | 1.7778 |
| 8 | 2.1485 | 2.6667 | 2.2105 |

The anomaly vanishes at $p=3$ — imposed by the theorem, not fitted. And the
shell model is *more* intermittent than real turbulence: reducing to one
amplitude per octave removes geometry. It reproduces the **existence** of the
anomaly, not its value. Stated in the paper.

</details>

<details>
<summary><b>10 — Lois de conservation hyperboliques et écoulements compressibles</b></summary>

<br>

**The question.** A shock is a discontinuity solving a PDE where derivatives do
not exist. How is that legitimate?

**What is proved.**

- **Classical solutions die.** $t_b=-1/\min u_0'$, from the vanishing of the
  characteristic Jacobian $1+t\,u_0'$. For $u_0=-\sin\pi x$: $t_b=1/\pi$,
  verified to $3\times10^{-6}$ %. The initial datum is analytic — the
  singularity is made by the equation.
- **Weak solutions are not unique** — *exhibited*, not invoked. For Burgers
  with $u_L=0$, $u_R=1$, both the rarefaction fan and the "expansion shock" of
  speed $s=1/2$ are weak solutions, and the second satisfies Rankine–Hugoniot
  exactly (residual zero, checked).
- **The second law selects.** Along the Hugoniot curve,
  $\Delta s/c_v=\ln P+\gamma\ln[(\beta P+1)/(P+\beta)]$ is strictly negative
  for $P<1$: expansion shocks are forbidden.
- **A derived corollary.** Expanding in $P=1+\varepsilon$, the $\varepsilon$
  *and* $\varepsilon^2$ terms cancel **exactly** (because $a-b=-1/\gamma$ and
  $a+b=1$), leaving

$$\frac{\Delta s}{c_v}=\frac{\gamma^2-1}{12\gamma^2}\,\varepsilon^3+O(\varepsilon^4)$$

  Coefficient $0.04081633$ for $\gamma=1.4$; measured ratio $0.9955$ at
  $\varepsilon=0.003$. **A weak shock is almost isentropic** — dissipation
  appears only at third order in amplitude. That is why linear acoustics works
  so well.

```
exact Riemann solver, Sod         p* = 0.303130178051   u* = 0.927452620049
Rankine–Hugoniot in shock frame
    mass       ρv                 discrepancy  5.6 × 10⁻¹⁷
    momentum   p + ρv²            discrepancy  1.1 × 10⁻¹⁶
    enthalpy   h + v²/2           discrepancy  8.9 × 10⁻¹⁶
entropy        +0.054 at the shock,  exactly 0 across the rarefaction

nozzle         A/A*(0.999) = A/A*(1.001) = 1.00000083   → M = 1 at a minimum
normal shock at M₁ = 2   0.577350 / 4.500000 / 0.720874
tables                   0.57735  / 4.50000  / 0.72087
shock position   p_b/p₀ = 0.5  →  A_s/A_t = 2.460494,  M₁ = 2.4256
```

**The throat is a saddle point.** Two branches pass through $M=1$ for each
$A/A^*$; a generic flow does not transit, only a precise back pressure selects
the separatrix. *That* is sonic choking.

**Open problem flagged.** In several space dimensions, uniqueness of admissible
weak solutions is **not** proved — De Lellis and Székelyhidi showed persistent
non-uniqueness for Euler.

</details>

---

## Verified results at a glance

Every number below is a *measurement against an independent prediction*, not a
fitted parameter.

| paper | quantity | prediction | measured | agreement |
|---|---|---|---|---|
| 02 | growth rate $\gamma(k)$ | kinetic dispersion relation | 6 modes | **0.70 %** |
| 02 | Penrose vs Faddeeva $k_c$ | two independent methods | — | **0.3 %** |
| 03 | rolling holonomy | $\Theta=\mathcal A/a^2$ (Lie bracket) | $\Theta/\mathcal A$ | **0.9997918** |
| 04 | LBM viscosity | $c_s^2(\tau-\tfrac12)$ | 8 values of $\tau$ | **0.156 %** |
| 04 | spatial convergence | order 2 | ratios | **4.04 · 4.04 · 3.98** |
| 05 | energy balance | $R+T=1$ | lossless metal | **1.000000000000** |
| 06 | exact neutral mode | $\operatorname{sech}y$ at $k=1$ | residual | **5.6 × 10⁻¹⁷** |
| 06 | max growth rate | $0.1897$ at $0.4446$ | $0.189701$ at $0.444939$ | **0.00 % / 0.08 %** |
| 07 | modified Hamiltonian | order $\Delta t^4$ | fitted order | **4.00** |
| 07 | angular momentum | exact | drift | **8.6 × 10⁻¹⁴** |
| 08 | FENE-P saturation | $L^2(1-1/2\mathrm{Wi})=500$ | $499.833$ | **0.03 %** |
| 09 | $\zeta_3$ | exactly 1 (the 4/5 law) | $0.98397$ | **1.60 %** |
| 10 | Rankine–Hugoniot | exact jump conditions | 3 relations | **10⁻¹⁶** |
| 01 | Hamiltonian constraint | exactly zero | never imposed | **3.9 × 10⁻¹²** |

---

## How the papers talk to each other

Not a collection — a corpus. Some threads run across several papers:

```
Squire's theorem
    06  PROVED — 2D disturbances dominate, so the reduction is legitimate
    02  ABSENT — no analogue in plasma, so 1D–1V remains an assumption
        → the same reduction is rigorous in one paper and a hypothesis in the other

Attractors and phase-space volume
    07  Hamiltonian ⟹ symmetric Lyapunov spectrum ±λᵢ ⟹ NO attractor
    03  Chaplygin sleigh HAS an attractor with energy exactly conserved
        → because non-holonomic systems are a third class, neither of the two

Dissipation independent of viscosity
    10  shock: ν sets the thickness, not the entropy produced
    08  Taylor–Culick: half the surface energy lost in a perfect fluid
    09  dissipative anomaly: ε finite as ν → 0
        → the same structural fact in three different physical settings

Anomalous dimensions
    09  ζₚ ≠ p/3 from fluctuations at all scales
    01  α′ corrections resummed to all orders, non-perturbative regime
        → the same field-theoretic notion, in turbulence and in string theory

Kinetic theory → macroscopic laws
    02  Vlasov: what survives when collisions are removed
    04  Boltzmann: how Navier–Stokes emerges when they dominate
        → the two opposite limits of the same hierarchy
```

---

## Reference implementations by other authors

Several directories contain **well-known reference codes written by others**,
kept here as the numerical baselines the papers are built against. All original
headers are intact. Full credit to the authors:

| directory | author | source |
|---|---|---|
| `plasma-two-stream-instability-pic/pic.py` | **Philip Mocz** (Princeton) | *Create Your Own Plasma PIC Simulation* |
| `Kelvin-Helmholtz-Instability/finitevolume.py` | **Philip Mocz** | *Create Your Own Finite Volume Fluid Simulation* |
| `lattice-boltzmann-fluid-simulation/latticeboltzmann.py` | **Philip Mocz** | *Create Your Own Lattice Boltzmann Simulation* |
| `celestial-dynamics-model .../nbody.py` | **Philip Mocz** | *Create Your Own N-body Simulation* |
| `elastic-taylor-culick-retraction/*.c` | **Vatsal Sanjay** (Physics of Fluids, Twente) | Basilisk C viscoelastic solvers |
| `droplet-fluid-dynamics-flow-simulation/*` | **Vatsal Sanjay** | Basilisk C two-phase solvers |
| `computational-physics-methods/` | **Hans Petter Langtangen** | *Programming for Computations* (book website) |

**My contribution in this repository is the ten papers and their code.** Every
simulation used in a paper — the PIC code, the D2Q9 solver, the RCWA
implementation, the Rayleigh eigensolver, the shell model, the exact Riemann
solver, the symplectic integrators, the rheological ODE solvers, the
cosmological integrator — is written independently and shipped with the paper.
Where a paper discusses the same physics as one of the reference codes above,
the reference code is the *baseline*, and the paper's contribution is the theory
and the quantitative validation.

---

## Reproducibility

Every figure and every number in every paper comes from a single script,
shipped alongside. No plasma, CFD, electromagnetic or rheology library is used —
only `numpy`, `scipy` and `matplotlib`.

```bash
git clone https://github.com/Koussay17/M.Koussay-engineering-portfolio.git
cd M.Koussay-engineering-portfolio
pip install numpy scipy matplotlib

# each paper directory contains its own script, e.g.
python plasma-two-stream-instability-pic/twostream.py   # theory + PIC
python plasma-two-stream-instability-pic/analyse.py     # fits + figures
```

Random seeds are fixed. LaTeX sources are included, so every paper recompiles:

```bash
pdflatex paper.tex && pdflatex paper.tex && pdflatex paper.tex
```

---

## Repository map

```
theoretical-physics-cosmology/       01  string cosmology  +  computational notebook
plasma-two-stream-instability-pic/   02  Vlasov–Poisson         [+ ref: Mocz]
solid-mechanics-geometric-reasoning/ 03  non-holonomic mechanics
lattice-boltzmann-fluid-simulation/  04  Chapman–Enskog         [+ ref: Mocz]
computational-diffraction-fourier/   05  Li's factorization rules
Kelvin-Helmholtz-Instability/        06  shear instability      [+ ref: Mocz]
celestial-dynamics-model .../        07  symplectic & chaos     [+ ref: Mocz]
elastic-taylor-culick-retraction/    08  viscoelasticity        [+ ref: Sanjay]
droplet-fluid-dynamics-flow-.../     08  two-phase flow         [+ ref: Sanjay]
advanced-fluid-thermal-sciences/     09  turbulence  +  experimental reports
aerospace-flow-physics/              10  hyperbolic conservation laws
dynamics-control-physical-systems/       control & system dynamics
scientific-computing-numerical-.../      numerical methods
energetics-fields-measurement-.../       energy & measurement systems
signals-communication-physical-.../      signals & communication
computational-physics-methods/           [ref: Langtangen textbook site]
```

---

## Toolchain

<div align="center">

![Python](https://img.shields.io/badge/Python-3776AB?style=flat-square&logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?style=flat-square&logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?style=flat-square&logo=scipy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557c?style=flat-square)
![LaTeX](https://img.shields.io/badge/LaTeX-008080?style=flat-square&logo=latex&logoColor=white)
![C](https://img.shields.io/badge/C-A8B9CC?style=flat-square&logo=c&logoColor=black)
![MATLAB](https://img.shields.io/badge/MATLAB-0076A8?style=flat-square)
![Basilisk](https://img.shields.io/badge/Basilisk_C-adaptive_multiphase-444?style=flat-square)

</div>

Numerical methods actually used across the papers: Dormand–Prince 8(5,3) and
Radau for stiff ODEs, Runge–Kutta 4, symplectic Störmer–Verlet, generalized
eigenvalue problems, spectral Poisson solvers, cloud-in-cell particle
deposition, Faddeeva-function root finding, Toeplitz operator inversion,
Benettin's variational algorithm.

---

<div align="center">

### Contact

**Koussay Mansouri** · Sousse, Tunisia

[![GitHub](https://img.shields.io/badge/GitHub-Koussay17-181717?style=for-the-badge&logo=github)](https://github.com/Koussay17)

<br>

*"Ce que ce travail ne fait pas est tout aussi important."*
Every paper ends with its own limits — stated, quantified, and left open.

</div>
