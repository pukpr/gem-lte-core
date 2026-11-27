--  GEM.Tidal_Poly: Meeus-free compact gravitational tidal polynomial forcing.
--
--  Replaces the 52-line empirical Doodson table (TIDES=TABLE) and the
--  Meeus ch. 47 ephemeris generator (TIDES=ZONAL) with a compact algebraic
--  polynomial basis in the 4 fundamental linear orbital mean angles:
--     L  (Moon mean longitude),   M' (Moon mean anomaly),
--     F  (Moon arg. of latitude), D  (Moon mean elongation),
--  plus derived Ls = L - D (Sun mean longitude), Omega = L - F (node), and
--  mean solar day phase theta_d for optional diurnal/semidiurnal tiers.
--
--  Enabled by TIDES=POLY (default TABLE = old behavior).
--  Tier selection via POLY_TIER (default Z10):
--    Z9          :  9 core zonal gravitational monomials (18 params)
--    Z10         : 10 zonal monomials (Z9 + 2D+M' [9.6137d], 20 params)
--    Z12         : 12 zonal monomials (Z10 + Mqm [6.859d] + Msqm [7.096d], 24 params)
--    Z15         : 15 zonal monomials (Z10 + L, F, M'-Om, 2L-M'-Om, Mtm, 30 params)
--    Z10_D6      : 16 monomials (Z10 + K1, O1, P1, M2, N2, K2, 32 params)
--    Z15_D6      : 21 monomials (Z15 + K1, O1, P1, M2, N2, K2, 42 params)
--    Z10_N5      : 15 monomials (Z10 + 5 multi-year nodal/perigee terms N5, 30 params)
--    Z10_D6_N5   : 21 monomials (Z10 + D6 + N5, 42 params)
--    Z10_S5_N5   : 20 monomials (Z10 + Mtm + 4 nodal sidebands S4 + N5, 40 params)
--    Z12_S5_N5   : 22 monomials (Z12 + Mtm + S4 + N5, 44 params)
--    Z10_S5_D6_N5: 26 monomials (Z10 + Mtm + S4 + D6 + N5, 52 params)
--    Z12_S5_D6_N5: 28 monomials (Z12 + Mtm + S4 + D6 + N5, 56 params)
--
--  Sub-annual terms (|period| <= 1 yr) are calibrated & anchored to dlod3.dat;
--  multi-year terms (Is_Long_Term, |period| > 1 yr) are left free of the dLOD
--  prior (with a light L2 ridge) to capture long-term decadal LOD / ocean drift.
--
--  Active coefficients are stored in slots 1 .. Active_Count of the standard
--  Long_Periods_Amp_Phase array (with slots Active_Count+1 .. 52 zeroed),
--  preserving full compatibility with .p files, TEST_ONLY, and Refine.

with GEM.LTE;

package GEM.Tidal_Poly is

   Max_Poly_Terms : constant Positive := 32;

   type LF_Vec is array (Integer range <>) of Long_Float;

   --  True when TIDES=POLY (case-insensitive).
   function On return Boolean;

   --  Active tier name (uppercase: Z9, Z10, Z15, Z10_D6, Z15_D6, Z10_N5, Z10_D6_N5).
   function Tier_Name return String;

   --  Number of active polynomial basis pairs (each pair has Amp, Phase -> 2*N params).
   function Active_Count return Positive;

   --  True if term K (1 .. Active_Count) is a multi-year term (|period| > 1 yr)
   --  that should be allowed to vary freely of the sub-annual dLOD prior.
   function Is_Long_Term (K : Positive) return Boolean;

   --  Short human-readable label for term K (1 .. Active_Count).
   function Term_Name (K : Positive) return String;

   --  Approximate period in days for term K (1 .. Active_Count).
   function Term_Period_Days (K : Positive) return Long_Float;

   --  Exact linear frequency in cycles/year for term K at Year_Len days/yr.
   function Term_Freq_CPY
     (K : Positive; Year_Len : Long_Float := GEM.LTE.Year_Length)
      return Long_Float;

   --  Evaluate the polynomial tidal forcing on decimal-year dates T using
   --  Constituents (Constituents'First .. Constituents'First + Active_Count - 1).
   --  Applies analytic jerk (Integ) and LOD (Lod) corrections per basis term.
   function Forcing_Series
     (T            : LF_Vec;
      Constituents : GEM.LTE.Long_Periods_Amp_Phase;
      Year_Len     : Long_Float := GEM.LTE.Year_Length;
      Integ        : Long_Float := 0.0;
      Lod          : Long_Float := 0.0) return LF_Vec;

   --  Calibrate the active polynomial basis on a dLOD file (e.g. dlod3.dat)
   --  via exact linear least squares. Returns a full Long_Periods_Amp_Phase
   --  array (1 .. GEM.LTE.Doodson_Args'Length) with slots 1 .. Active_Count
   --  filled and remaining slots zeroed.
   function Calibrate_dLOD
     (File_Name       : String;
      Year_Correction : Long_Float := 0.0)
      return GEM.LTE.Long_Periods_Amp_Phase;

end GEM.Tidal_Poly;
