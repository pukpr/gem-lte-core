# Theoretical Framework: Mapping Landreman 3D MHD Toroidal Equilibria to Geophysical Fluid Dynamics (GEM-LTE)

## 1. Executive Summary & Core Premise

This document establishes the mathematical and physical equivalence between Matt Landreman's closed-form, analytic 3D Magnetohydrodynamic (MHD) toroidal equilibria (*arXiv:2609.26742*) and non-autonomous Geophysical Fluid Dynamics (GFD) models implemented in the `pukpr/gem-lte-core` framework.

The primary premise is that large-scale, long-term climate oscillations—specifically the **El Niño–Southern Oscillation (ENSO)** and the **Quasi-Biennial Oscillation (QBO)**—can be modeled as topologically protected, 3D toroidal standing waves on invariant manifolds. Rather than being chaotic, autonomous internal fluid instabilities, these structures are **non-autonomous, paced systems** whose spatial shapes are dictated by Landreman-like Euler/MHD equilibrium surfaces and whose temporal phase trajectories are driven by astronomical tidal forcing (draconic, tropical, sidereal, and solar-annual cycles).

---

## 2. Isomorphism Between Magnetohydrodynamics (MHD) and Geophysical Fluid Dynamics (GFD)

The mathematical bridge connecting Landreman's plasma equilibria to GFD relies on the exact structural isomorphism between steady, incompressible magnetohydrodynamics and inviscid, rotational geostrophic fluid flows.

### 2.1 The Vector Field Mapping
In steady-state magnetostatic equilibrium, force balance is governed by:
$$\nabla p = (\nabla \times \mathbf{B}) \times \mathbf{B}$$

In a rotating frame of reference (under the geostrophic/Euler limit), fluid momentum balance is governed by:
$$(\mathbf{u} \cdot \nabla)\mathbf{u} + 2\mathbf{\Omega} \times \mathbf{u} = -\nabla \Phi - \frac{1}{\rho}\nabla P_{\text{fluid}}$$

By defining a **Beltrami flow condition** where the fluid velocity field $\mathbf{u}$ is aligned with the magnetic field vector $\mathbf{B}$ ($\mathbf{u} = \mathbf{B} / \sqrt{\mu_0 \rho}$), the advective tension $(\mathbf{u} \cdot \nabla)\mathbf{u}$ and Lorentz force $(\mathbf{B} \cdot \nabla)\mathbf{B}$ map identically.

### 2.2 Physical Field Analogs

| Physical Concept | Landreman 3D MHD Model | GEM-LTE Geophysical Model |
| :--- | :--- | :--- |
| **Primary Vector Field** | Magnetic Field Vector $\mathbf{B}$ | Fluid Velocity Vector $\mathbf{u}$ |
| **Field Curl / Density** | Current Density $\mathbf{J} = \nabla \times \mathbf{B}$ | Vorticity Vector $\boldsymbol{\omega} = \nabla \times \mathbf{u}$ |
| **Invariant Boundary** | Iso-Pressure Surface $p(\psi) = \text{const}$ | Potential Vorticity / Thermocline Surface $\Psi = \text{const}$ |
| **Field Pitch / Winding** | Rotational Transform $\iota = \frac{d\psi_{\text{pol}}}{d\psi_{\text{tor}}}$ | Phase Winding Ratio $\iota = \frac{\omega_{\text{poloidal}}}{\omega_{\text{toroidal}}}$ |
| **Conserved Quantity** | Magnetic Helicity $H_M = \int \mathbf{A} \cdot \mathbf{B} \, dV$ | Kinematic/Potential Helicity $H_K = \int \mathbf{u} \cdot \boldsymbol{\omega} \, dV$ |

---

## 3. Reconciling Temporal Phase Space (GEM-LTE) and Spatial Toroidal Geometry (Landreman)

A standing wave is a space-time separable field:
$$\Psi(x, y, z, t) = S(x, y, z) \cdot T(t)$$

1. **Spatial Container $S(x, y, z)$ (Landreman):** The physical geometry of the ocean basin or global atmospheric envelope acts as a boundary condition fixing the 3D topological manifold (the "static geometric track" or invariant torus).
2. **Temporal Phase Trajectory $T(t)$ (GEM-LTE):** The astronomical tidal forcing drives the movement along that static track over time via non-linear phase modulation:
   $$\mathbf{r}(t) = \begin{bmatrix} A \cos(M \cdot \omega_0 \sin(t)) \\ A \sin(M \cdot \omega_0 \sin(t)) \end{bmatrix}$$
3. **Equivalence:** Increasing the temporal modulation index $M$ in GEM-LTE increases the spatial pitch/twist (rotational transform $\iota$) of the Landreman streamtubes in 3D physical space. High temporal phase velocity modulation manifests as tight spatial coils (vorticity concentration).

---

## 4. Dual Classification of Climate Modes

The framework divides large-scale geophysical oscillations into two distinct spatial wavenumber classes governed by domain geometry.

```
                         Geophysical Domain Geometry
                                      │
              ┌───────────────────────┴───────────────────────┐
              ▼                                               ▼
   Unbounded Global Ring                                Bounded Ocean Basin
      (k_ϕ = 0 Mode)                                     (k_x ≠ 0 Mode)
              │                                               │
              ▼                                               ▼
 QBO (Stratospheric Ring)                           ENSO (Pacific Basin)
• Axisymmetric (Grad-Shafranov limit)              • Non-axisymmetric 3D Helical Torus
• Driven by Draconic Lunar Cycle                   • Driven by Tropical/Sidereal/Solar Tides
• ~28-month pacing beat                            • ~3-7 year pacing beat
```

### 4.1 The QBO: Axisymmetric $k_\phi = 0$ Complete Toroidal Mode
* **Domain & Geometry:** Global equatorial stratospheric band wrapping around the globe ($0 \le \phi < 2\pi$) without land barriers.
* **Spatial Symmetry:** Strict zonal axisymmetry ($k_\phi = 0$). Landreman's 3D MHD equations reduce to the axisymmetric **Grad-Shafranov Equation** (or the GFD Eliassen-Palm flux streamfunction relation).
* **Forcing & Pacing:** Because $k_\phi = 0$, horizontal spatial gradients vanish. The system selects equatorially symmetric, $k=0$ global impulses—specifically the **draconic lunar month** ($T \approx 27.2122\text{ days}$, nodal crossing cycle).
* **Mechanism:** The non-autonomous beating between high-frequency $k=0$ draconic impulses and annual/seasonal harmonics ($T = 1.0\text{ year}$) paces the strict **~28-month (~2.33 year)** downward phase descent and velocity flipping of the stratospheric zonal winds.

### 4.2 The ENSO Mode: Non-Axisymmetric $k_x \neq 0$ Bounded Basin Mode
* **Domain & Geometry:** Tropical Pacific Ocean bounded by physical continental margins at $x = \pm L_x$.
* **Spatial Symmetry:** Non-axisymmetric standing mode ($k_x \neq 0$). Closed circulation forms a Landreman-style 3D non-axisymmetric helical streamtube where eastward equatorial Kelvin waves (depth $z \approx -150\text{ m}$) reflect off boundaries into westward off-equatorial Rossby waves ($z \approx -50\text{ m}$ to $0\text{ m}$).
* **Forcing & Pacing:** The fixed spatial boundary condition ($u_x = 0$ at $x = \pm L_x$) breaks zonal axisymmetry, unlocking coupling to **tropical lunar declination cycles** ($T \approx 27.3216\text{ days}$) and **sidereal lunar orbital cycles** ($T \approx 27.3217\text{ days}$).
* **Mechanism:** The multi-harmonic interaction between tropical, sidereal, draconic, and annual solar components paces the characteristic **3-to-7-year quasi-periodic ENSO recurrence band**.

---

## 5. Non-Autonomous Stability Enforcement

In ideal MHD (Landreman), stability is **unforced and topological**, relying on exact mathematical integrability ($\nabla p = \mathbf{J} \times \mathbf{B}$) to prevent field line chaos. 

In real, dissipative geophysical fluids, stability is **externally enforced and non-autonomous**:
1. The system behaves as a **driven non-autonomous attractor**.
2. External, deterministic astronomical tidal forcing continuously "re-anchors" the phase trajectory on the 3D torus before viscous dissipation or turbulent shear can collapse the surface topology.
3. The rotational transform $\iota$ is directly set by the ratio of astronomical forcing frequencies ($\iota = \omega_{\text{tidal}} / \omega_{\text{annual}}$).

---

## 6. Mathematical Pipeline: GEM-LTE to Landreman Surface Construction

To programmatically convert GEM-LTE winding parameters into 3D physical fields:

1. **Extract Winding Ratio:** Obtain $\omega_{\text{poloidal}}$ and $\omega_{\text{toroidal}}$ from `gem-lte-core` harmonic analysis.
2. **Set Rotational Transform:** Compute $\iota = \frac{\omega_{\text{poloidal}}}{\omega_{\text{toroidal}}}$.
3. **Map 3D Coordinates:**
   $$X(\phi, \theta) = \left( R_0 + r \cos(\theta) \right) \cos(\phi)$$
   $$Y(\phi, \theta) = \left( R_0 + r \cos(\theta) \right) \sin(\phi)$$
   $$Z(\phi, \theta) = r \sin(\theta)$$
   where poloidal angle is phase-locked via $\theta = \iota \cdot \phi + M \sin(\phi)$.
4. **Evaluate Equilibrium Vector Field:**
   $$u_\phi = U_0 \left( \frac{R_0}{R} \right), \quad u_\theta = u_\phi \left( \frac{r}{R_0} \right) \iota$$
   $$\mathbf{u}(x,y,z) = -u_\phi \sin(\phi)\mathbf{e}_x + u_\phi \cos(\phi)\mathbf{e}_y + u_\theta \cos(\theta)\mathbf{e}_z$$

---

## 7. Conclusion & Methodological Takeaway

By unifying Landreman's 3D MHD solutions with GEM-LTE's tidal winding framework:
* **Landreman provides the spatial geometry:** Invariant 3D toroidal surfaces that prevent turbulent breakdown.
* **GEM-LTE provides the temporal driver:** Non-autonomous astronomical tidal pacing (draconic for $k=0$ QBO; tropical/sidereal/solar for $k \neq 0$ ENSO).

Together, they prove that complex, low-frequency climate oscillations are mathematically structured, non-autonomous 3D standing waves anchored to invariant toroidal manifolds.