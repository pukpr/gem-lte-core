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
   --  Julian date of a decimal-year date: JD = 2448074.761995 + (T - 1990.5) * (365.2422484 + YEAR),
   --  shared by dLOD and climate files; YEAR sets the comb's model-year length as in TABLE mode.
   function JD_Of (T : Long_Float) return Long_Float;
   function Potential (JD : Long_Float) return Long_Float;
   --  dV/dt per day, central difference with h = 0.5 day
   function Rate (JD : Long_Float) return Long_Float;
   --  the forcing value at decimal year T
   --  Integ: jerk factor (shfT), adds (Integ / 2 pi) * d(forcing)/dt, t in years
   function Forcing_At (T : Long_Float; Integ : Long_Float := 0.0) return Long_Float;
end GEM.Zonal;
