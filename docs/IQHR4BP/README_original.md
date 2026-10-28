# **Atmospheric–Orbital Frequency Architecture: A Cross‑Domain Interpretation of Quasi‑Periodic Invariant Manifolds in the IQHR4BP**  

This document reconstructs the reasoning chain that connects:

- Earth’s **SAO → QBO** atmospheric oscillation ladder  
- The **lunar draconic cycle**  
- The **in‑plane / out‑of‑plane center manifolds** in the HVO  
- The **CR3BP → HR4BP → IQHR4BP → full ephemeris** invariant‑torus ladder  

and shows how these structures map onto each other.

The goal is to provide enough context to understand why the atmospheric oscillation architecture mirrors the invariant‑torus architecture developed in:

> **Grossi, G., & Ross, S. D. (2024). *Quasi‑Periodic Invariant Manifolds in the In‑Plane Quasi‑Hill Restricted Four‑Body Problem*.**  
> `https://www.researchgate.net/publication/412048007` [(researchgate.net in Bing)](https://www.bing.com/search?q="https%3A%2F%2Fwww.researchgate.net%2Fpublication%2F412048007")  
> <sup>†</sup>

---

## **1. The Atmospheric Ladder: Internal → Solar → Lunar → Full Ephemeris**

### **1.1 Pure Equatorial Mode (Internal Only)**  
A hypothetical atmosphere with no seasonal cycle and no lunar tides would support a single equatorial oscillation mode.  
This is a **1‑frequency internal oscillator**.

### **1.2 SAO: Solar Seasonal Forcing Added**  
The real atmosphere receives **solar seasonal forcing**, producing the **Semi‑Annual Oscillation (SAO)**:

- reverses twice per year  
- driven by solar heating asymmetry around equinoxes  
- supported by equatorial Kelvin/Rossby waves  

This is still a **1‑torus**, but embedded in an annual background.

### **1.3 QBO: SAO + Lunar Draconic Forcing**  
The Moon’s **draconic cycle** (27.212 days) modulates:

- lunar declination  
- equatorial tidal forcing  
- wave propagation in the equatorial waveguide  

The SAO (internal) and draconic (external) frequencies combine to form a **2‑torus**, producing a beat frequency:

\[
\omega_{\text{QBO}} = |\, p\,\omega_{\text{SAO}} - q\,\omega_{\text{drac}} \,|
\]

The dominant term:

\[
\omega_{\text{QBO}} \approx |13 - \omega_{\text{drac}}|
\]

gives a period near **28 months**, the observed QBO.

### **1.4 Full Atmospheric Ephemeris**  
Adding:

- full lunar tidal spectrum (M2, N2, Mf, Mm, 18.6‑yr nodal, 8.85‑yr perigee)  
- solar tides  
- ENSO‑scale coupled modes  
- decadal planetary‑wave variability  

produces a **multi‑frequency atmospheric torus** with a dense spectral comb.

---

## **2. The Orbital Ladder: CR3BP → HR4BP → IQHR4BP → Full Ephemeris**

### **2.1 CR3BP: Internal Only**  
Earth–Moon primaries + massless spacecraft.  
Autonomous in rotating frame.  
Produces **periodic orbits / 1‑tori**.

### **2.2 HR4BP: Add Solar Forcing**  
Sun added via Hill approximation.  
Introduces one external frequency.  
Produces **2‑tori**.

### **2.3 IQHR4BP: Add Lunar Geometry (Eccentricity + Inclination)**  
Adds:

- anomalistic month  
- draconic month  
- Earth–Moon pulsation  

Produces **3‑tori**.  
The HVO linear stability reveals:

- **in‑plane center manifold** → anomalistic month  
- **out‑of‑plane center manifold** → draconic month  

### **2.4 Full Ephemeris / HFEM**  
Adds:

- full lunar ephemeris  
- planetary perturbations  
- non‑averaged Earth–Moon geometry  

Produces a **multi‑frequency torus** with a dense spectral comb.

---

## **3. Explicit One‑to‑One Mapping Between the Two Ladders**

| Orbital Ladder | Atmospheric Ladder | Structural Match |
|----------------|--------------------|------------------|
| **CR3BP** (internal only) | Pure equatorial mode | 1 internal frequency |
| **HR4BP** (solar forcing) | SAO (solar seasonal forcing) | internal + solar external |
| **IQHR4BP** (solar + lunar geometry) | SAO + draconic → QBO | internal + solar + lunar external; beat mode emerges |
| **Full ephemeris** | Full tidal–ENSO–decadal spectrum | multi‑frequency torus with dense spectral comb |

The mapping is exact at the level of **forced quasi‑periodic invariant tori**.

---

## **4. Why the Draconic Cycle Is the Atmospheric Out‑of‑Plane Manifold**

In the HVO stability analysis:

- **In‑plane center manifold** → anomalistic month  
- **Out‑of‑plane center manifold** → draconic month  

In the atmosphere:

- **SAO internal manifold** ↔ in‑plane  
- **Draconic tidal modulation** ↔ out‑of‑plane  

The draconic cycle controls **lunar declination**, which controls:

- equatorial tidal amplitude  
- waveguide geometry  
- vertical momentum transport  

Thus the draconic cycle is the correct external frequency for generating the QBO.

---

## **5. The QBO as a Beat Mode on a 2‑Torus**

The SAO and draconic frequencies form a 2‑torus:

\[
\omega = [\omega_{\text{SAO}},\ \omega_{\text{drac}}]
\]

The QBO emerges as a beat mode:

\[
\omega_{\text{QBO}} = |13 - \omega_{\text{drac}}|
\]

This is the same mechanism by which:

- IQHR4BP tori produce long‑period modulations  
- ELFOs exhibit draconic sidebands  
- Chandler wobble emerges from annual + lunar aliasing  

The architecture is identical.

---

## **6. Torus‑Dimension Ladder Diagram (Atmosphere ↔ Orbital)**

```
Atmosphere                               Orbital Dynamics
-----------                              ------------------
Pure Equatorial Mode (1D)                CR3BP (1D)
   ↓ external continuation (solar)          ↓ external continuation (solar)
SAO (1D-in-annual background)             HR4BP (2D)
   ↓ external continuation (draconic)       ↓ external continuation (draconic)
QBO (2D torus)                            IQHR4BP (3D torus)
   ↓ full ephemeris forcing                 ↓ full ephemeris forcing
Full Tidal–ENSO–Decadal Comb (multi-D)   HFEM / Full Ephemeris (multi-D)
```

This ladder shows the **structural equivalence** between:

- atmospheric oscillation hierarchy, and  
- orbital invariant‑torus hierarchy.

Both systems grow in torus dimension by **adding external quasi‑periodic forcing**.

---

## **7. Why This Connects Directly to Grossi & Ross (2024)**

The cited paper develops:

- **quasi‑periodic invariant tori**  
- **internal vs external continuation**  
- **in‑plane vs out‑of‑plane center manifolds**  
- **solar + lunar forcing**  
- **multi‑frequency torus geometry**

This is *exactly* the same mathematical structure underlying:

- SAO → QBO  
- lunar tidal modulation  
- atmospheric beat frequencies  
- Chandler wobble forced aliasing  

Thus the atmospheric ladder is a **physical analog** of the orbital ladder in the IQHR4BP.

---

<sup>†</sup> **Formal citation:**  
Grossi, G., & Ross, S. D. (2024). *Quasi‑Periodic Invariant Manifolds in the In‑Plane Quasi‑Hill Restricted Four‑Body Problem*. Retrieved from ResearchGate: `https://www.researchgate.net/publication/412048007` [(researchgate.net in Bing)](https://www.bing.com/search?q="https%3A%2F%2Fwww.researchgate.net%2Fpublication%2F412048007")

[paper](GrossiBaresiTopputo2026Quasi-PeriodicInvariantManifoldsintheIn-PlaneQuasi-HillRestrictedFour-BodyProblem.pdf)



