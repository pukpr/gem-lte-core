--  GEM.Zonal: closed-form zonal (long-period) tidal forcing (2026-10-06).
--
--  Replaces the 42-constituent table with the degree-2 zonal tide potential
--     V(t) = sum over Moon, Sun of w_b * (abar_b / r_b)**3 * (3 sin**2 dec_b - 1) / 2,
--  Moon from the full Meeus ch. 47 series (lteMod gem-ephemeris.adb, fixed) in the
--  fundamental arguments D, M, M', F; Sun from its Keplerian orbit. Every tidal
--  constituent (and its nodal sidebands and cross terms) follows from the
--  expansion, so no per-constituent amplitude or phase is free.
--  Enabled by TIDES=ZONAL (default TABLE = old behaviour). The forcing is the
--  RATE dV/dt (dLOD is a rate), scaled and aligned to the daily dLOD record:
--     forcing(t) = ZONAL_KAPPA * dV/dt (JD_Of(t) - ZONAL_TAU).
--  Tier-1 generating-factor multipliers (default 1.0): ZONAL_WSUN (Sun/Moon
--  weight), ZONAL_ECC (lunar eccentricity term), ZONAL_EVECT (evection),
--  ZONAL_VAR (variation), ZONAL_ANNEQ (annual equation), ZONAL_INCL (inclination).
--  Reference implementation: experiments/Oct2026_generator/zonal_generator.py.
package GEM.Zonal is
   function On return Boolean;
   --  TIDES=HYBRID: On is True as well, and Tide_Sum adds the constituent table
   --  (lpap) to the generator as correction lines. The corrections are
   --  regularized in the search metric: score * (1 - HYBRID_LAMBDA * (sum of
   --  A**2 / 2) / Forcing_Variance), i.e. HYBRID_LAMBDA times the corrections'
   --  share of the generator's variance (default 1.0; 0 = no penalty).
   function Hybrid return Boolean;
   --  TIDES=BLEND: On is True as well, and Tide_Sum returns
   --     Rho * generator + (1 - Rho) * s * table,   Rho = the searched .p parameter
   --  "rho" (1.0 when the file has none),
   --  with the table (lpap) scaled by s = sqrt (Forcing_Variance / (sum of A**2 / 2))
   --  so both parts are normalized to the generator's variance. Rho = 1 is
   --  exactly TIDES=ZONAL. Use LOCKT=TRUE to hold the table at its dLOD values.
   --  The table must be calibrated at the base model year (.p year = 0); its
   --  clock is anchored at 1990.5 like the generator's.
   function Blend return Boolean;
   function Forcing_Variance return Long_Float;
   --  Julian date of a decimal-year date: JD = 2448074.761995 + (T - 1990.5) * (365.2422484 + YEAR),
   --  shared by dLOD and climate files; YEAR sets the comb's model-year length as in TABLE mode.
   --  Year_Len > 0: days per calendar year for this call (lt.exe passes
   --  Year_Length (D.B.Year), so the .p "year" term retunes the clock);
   --  0 = the startup clock (365.2422484 + resp YEAR).
   function JD_Of (T : Long_Float; Year_Len : Long_Float := 0.0) return Long_Float;
   function Potential (JD : Long_Float) return Long_Float;
   --  dV/dt per day, central difference with h = 0.5 day
   function Rate (JD : Long_Float) return Long_Float;
   --  the forcing value at decimal year T
   --  Integ: jerk factor (shfT), adds (Integ / 2 pi) * d(forcing)/dt, t in years
   --  Lod: LOD factor (resp LOD_TERM), adds Lod * 2 pi * (time integral of the
   --  forcing, t in years) = Lod * 2 pi * Kappa * (V - mean V) / (days per year):
   --  the forcing is the dLOD rate, this mixes in a little of LOD itself. Per
   --  line it is Lod / Freq * Amplitude in quadrature (the jerk is Integ * Freq).
   function Forcing_At
     (T : Long_Float; Integ : Long_Float := 0.0; Year_Len : Long_Float := 0.0;
      Lod : Long_Float := 0.0) return Long_Float;
   --  Forcing_At for every date in T, bit-identical to calling it per date.
   --  The generator's own values (rate, its derivative, potential) depend only
   --  on the dates and Year_Len, so they are cached per task (two date/year
   --  sets): a search step that leaves the year alone costs three multiplies
   --  per date instead of up to seven evaluations of the lunar series.
   type LF_Vec is array (Integer range <>) of Long_Float;
   function Forcing_Series
     (T : LF_Vec; Integ : Long_Float := 0.0; Year_Len : Long_Float := 0.0;
      Lod : Long_Float := 0.0) return LF_Vec;
end GEM.Zonal;
