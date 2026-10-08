with Ada.Task_Attributes;
with Ada.Unchecked_Deallocation;
with Ada.Numerics.Long_Elementary_Functions;
package body GEM.Zonal is
   --  v2 (2026-10-06): full Meeus ch. 47 tables from lteMod gem-ephemeris.adb, with three fixes
   --  (a duplicated latitude argument corrected to 2D-M-M'+F; L in radians inside the L-F, L,
   --  L+-M' terms; obliquity T^2 coefficient -1.638889e-7 deg), the eccentricity factor E on
   --  M terms, polynomial fundamental arguments; clock: see Year_Days below.
   --  Generated from experiments/Oct2026_generator/meeus_full.json (reference: zonal_full.py).
   use Ada.Numerics.Long_Elementary_Functions;
   D2R : constant Long_Float := Ada.Numerics.Pi / 180.0;

   --  TIDES=HYBRID: the generator plus the constituent table as correction lines
   --  (see Tide_Sum); everything keyed on Zonal_On behaves as in TIDES=ZONAL.
   Hybrid_On : constant Boolean := GEM.Getenv ("TIDES", "TABLE") = "HYBRID";
   --  TIDES=BLEND: rho * generator + (1 - rho) * the constituent table (lpap, e.g.
   --  the dLOD regression), the table first scaled to the generator's variance
   --  (see Tide_Sum); rho is the searched .p parameter "rho".
   Blend_On : constant Boolean := GEM.Getenv ("TIDES", "TABLE") = "BLEND";
   Zonal_On : constant Boolean :=
     Hybrid_On or else Blend_On or else GEM.Getenv ("TIDES", "TABLE") = "ZONAL";
   Kappa : constant Long_Float := GEM.Getenv ("ZONAL_KAPPA", -5.866344);
   Tau   : constant Long_Float := GEM.Getenv ("ZONAL_TAU", 0.92);
   W_Sun : constant Long_Float := 0.4600 * GEM.Getenv ("ZONAL_WSUN", 1.0);
   M_Ecc   : constant Long_Float := GEM.Getenv ("ZONAL_ECC", 1.0);
   M_Evect : constant Long_Float := GEM.Getenv ("ZONAL_EVECT", 1.0);
   M_Var   : constant Long_Float := GEM.Getenv ("ZONAL_VAR", 1.0);
   M_AnnEq : constant Long_Float := GEM.Getenv ("ZONAL_ANNEQ", 1.0);
   M_Incl  : constant Long_Float := GEM.Getenv ("ZONAL_INCL", 1.0);

   --  Clock: the tides are evaluated at true Julian dates; the model (and its annual
   --  impulse comb) advances one MODEL year per decimal year, of length
   --  365.2422484 + YEAR days, exactly as in TIDES=TABLE. YEAR therefore keeps its
   --  production role of setting the comb aliases, e.g. the Mt (9.13295 d) beat that
   --  scales the 120 -> 60 yr cycle (production YEAR 0.00405305 -> 127.3 yr). The
   --  mapping is anchored at 1990.5 (mid-dLOD record) so the dLOD calibration (tau)
   --  is stable as YEAR changes.
   T_Anchor  : constant Long_Float := 1990.5;
   JD_Anchor : constant Long_Float := 2448074.761995;   -- 1880-01-01 + 110.5 tropical years
   Year_Days : constant Long_Float := 365.2422484 + GEM.Getenv ("YEAR", 0.0);

   --  argument multipliers in the order D, M, M', F, A1, A2, A3, L
   type Mults is array (1 .. 8) of Integer;
   type Vals is array (1 .. 8) of Long_Float;
   type Term is record
      Mu : Mults;
      C : Long_Float;
      Is_Sin : Boolean;
   end record;
   type Terms is array (Positive range <>) of Term;
   Lon : constant Terms := (
      ((0, 0, 1, 0, 0, 0, 0, 0), 6.288774, True),
      ((2, 0, -1, 0, 0, 0, 0, 0), 1.274027, True),
      ((2, 0, 0, 0, 0, 0, 0, 0), 0.658314, True),
      ((0, 0, 2, 0, 0, 0, 0, 0), 0.213618, True),
      ((0, 1, 0, 0, 0, 0, 0, 0), -0.185116, True),
      ((0, 0, 0, 2, 0, 0, 0, 0), -0.114332, True),
      ((2, 0, -2, 0, 0, 0, 0, 0), 0.058793, True),
      ((2, -1, -1, 0, 0, 0, 0, 0), 0.057066, True),
      ((2, 0, 1, 0, 0, 0, 0, 0), 0.053322, True),
      ((2, -1, 0, 0, 0, 0, 0, 0), 0.045758, True),
      ((0, 1, -1, 0, 0, 0, 0, 0), -0.040923, True),
      ((1, 0, 0, 0, 0, 0, 0, 0), -0.03472, True),
      ((0, 1, 1, 0, 0, 0, 0, 0), -0.030383, True),
      ((2, 0, 0, -2, 0, 0, 0, 0), 0.015327, True),
      ((0, 0, 1, 2, 0, 0, 0, 0), -0.012528, True),
      ((0, 0, 1, -2, 0, 0, 0, 0), 0.01098, True),
      ((4, 0, -1, 0, 0, 0, 0, 0), 0.010675, True),
      ((0, 0, 3, 0, 0, 0, 0, 0), 0.010034, True),
      ((4, 0, -2, 0, 0, 0, 0, 0), 0.008548, True),
      ((2, 1, -1, 0, 0, 0, 0, 0), -0.007888, True),
      ((2, 1, 0, 0, 0, 0, 0, 0), -0.006766, True),
      ((1, 0, -1, 0, 0, 0, 0, 0), -0.005163, True),
      ((1, 1, 0, 0, 0, 0, 0, 0), 0.004987, True),
      ((2, -1, 1, 0, 0, 0, 0, 0), 0.004036, True),
      ((2, 0, 2, 0, 0, 0, 0, 0), 0.003994, True),
      ((4, 0, 0, 0, 0, 0, 0, 0), 0.003841, True),
      ((2, 0, -3, 0, 0, 0, 0, 0), 0.003665, True),
      ((0, 1, -2, 0, 0, 0, 0, 0), -0.002689, True),
      ((2, 0, -1, 2, 0, 0, 0, 0), -0.002602, True),
      ((2, -1, -2, 0, 0, 0, 0, 0), 0.00239, True),
      ((1, 0, 1, 0, 0, 0, 0, 0), -0.002348, True),
      ((2, -2, 0, 0, 0, 0, 0, 0), 0.002236, True),
      ((0, 1, 2, 0, 0, 0, 0, 0), -0.00212, True),
      ((0, 2, 0, 0, 0, 0, 0, 0), -0.002069, True),
      ((2, -2, -1, 0, 0, 0, 0, 0), 0.002048, True),
      ((2, 0, 1, -2, 0, 0, 0, 0), -0.001773, True),
      ((2, 0, 0, 2, 0, 0, 0, 0), -0.001595, True),
      ((4, -1, -1, 0, 0, 0, 0, 0), 0.001215, True),
      ((0, 0, 2, 2, 0, 0, 0, 0), -0.00111, True),
      ((3, 0, -1, 0, 0, 0, 0, 0), -0.000892, True),
      ((2, 1, 1, 0, 0, 0, 0, 0), -0.00081, True),
      ((4, -1, -2, 0, 0, 0, 0, 0), 0.000759, True),
      ((0, 2, -1, 0, 0, 0, 0, 0), -0.000713, True),
      ((2, 2, -1, 0, 0, 0, 0, 0), -0.0007, True),
      ((2, 1, 2, 0, 0, 0, 0, 0), 0.000691, True),
      ((2, -1, 0, -2, 0, 0, 0, 0), 0.000596, True),
      ((4, 0, 1, 0, 0, 0, 0, 0), 0.000549, True),
      ((0, 0, 4, 0, 0, 0, 0, 0), 0.000537, True),
      ((4, -1, 0, 0, 0, 0, 0, 0), 0.00052, True),
      ((1, 0, -2, 0, 0, 0, 0, 0), -0.000487, True),
      ((2, 1, 0, -2, 0, 0, 0, 0), -0.000399, True),
      ((0, 0, 2, -2, 0, 0, 0, 0), -0.000381, True),
      ((1, 1, 1, 0, 0, 0, 0, 0), 0.000351, True),
      ((3, 0, -2, 0, 0, 0, 0, 0), -0.00034, True),
      ((4, 0, -3, 0, 0, 0, 0, 0), 0.00033, True),
      ((2, -1, -2, 0, 0, 0, 0, 0), 0.000327, True),
      ((0, 2, 1, 0, 0, 0, 0, 0), -0.000323, True),
      ((1, 1, -1, 0, 0, 0, 0, 0), 0.000299, True),
      ((2, 0, 3, 0, 0, 0, 0, 0), 0.000294, True),
      ((0, 0, 0, 0, 1, 0, 0, 0), 0.003958, True),
      ((0, 0, 0, -1, 0, 0, 0, 1), 0.001962, True),
      ((0, 0, 0, 0, 0, 1, 0, 0), 0.000318, True));
   Lat : constant Terms := (
      ((0, 0, 0, 1, 0, 0, 0, 0), 5.128122, True),
      ((0, 0, 1, 1, 0, 0, 0, 0), 0.280602, True),
      ((0, 0, 1, -1, 0, 0, 0, 0), 0.277693, True),
      ((2, 0, 0, -1, 0, 0, 0, 0), 0.173237, True),
      ((2, 0, -1, 1, 0, 0, 0, 0), 0.055413, True),
      ((2, 0, -1, -1, 0, 0, 0, 0), 0.046271, True),
      ((2, 0, 0, 1, 0, 0, 0, 0), 0.032573, True),
      ((0, 0, 2, 1, 0, 0, 0, 0), 0.017198, True),
      ((2, 0, 1, -1, 0, 0, 0, 0), 0.009266, True),
      ((0, 0, 2, -1, 0, 0, 0, 0), 0.008822, True),
      ((2, -1, 0, -1, 0, 0, 0, 0), 0.008216, True),
      ((2, 0, -2, -1, 0, 0, 0, 0), 0.004324, True),
      ((2, 0, 1, 1, 0, 0, 0, 0), 0.0042, True),
      ((2, 1, 0, -1, 0, 0, 0, 0), -0.003359, True),
      ((2, -1, -1, 1, 0, 0, 0, 0), 0.002463, True),
      ((2, -1, 0, 1, 0, 0, 0, 0), 0.002211, True),
      ((2, -1, -1, -1, 0, 0, 0, 0), 0.002065, True),
      ((0, 1, -1, -1, 0, 0, 0, 0), -0.00187, True),
      ((4, 0, -1, -1, 0, 0, 0, 0), 0.001828, True),
      ((0, 1, 0, 1, 0, 0, 0, 0), -0.001794, True),
      ((0, 0, 0, 3, 0, 0, 0, 0), -0.001749, True),
      ((0, 1, -1, 1, 0, 0, 0, 0), -0.001565, True),
      ((1, 0, 0, 1, 0, 0, 0, 0), -0.001491, True),
      ((0, 1, 1, 1, 0, 0, 0, 0), -0.001475, True),
      ((0, 1, 1, -1, 0, 0, 0, 0), -0.00141, True),
      ((0, 1, 0, -1, 0, 0, 0, 0), -0.001344, True),
      ((1, 0, 0, -1, 0, 0, 0, 0), -0.001335, True),
      ((0, 0, 3, 1, 0, 0, 0, 0), 0.001107, True),
      ((4, 0, 0, -1, 0, 0, 0, 0), 0.001021, True),
      ((4, 0, -1, 1, 0, 0, 0, 0), 0.000833, True),
      ((0, 0, 1, -3, 0, 0, 0, 0), 0.000777, True),
      ((4, 0, -2, 1, 0, 0, 0, 0), 0.000671, True),
      ((2, 0, 0, -3, 0, 0, 0, 0), 0.000607, True),
      ((2, 0, 2, -1, 0, 0, 0, 0), 0.000596, True),
      ((2, -1, 1, -1, 0, 0, 0, 0), 0.000491, True),
      ((2, 0, -2, 1, 0, 0, 0, 0), -0.000451, True),
      ((0, 0, 3, -1, 0, 0, 0, 0), 0.000439, True),
      ((2, 0, 2, 1, 0, 0, 0, 0), 0.000422, True),
      ((2, 0, -3, -1, 0, 0, 0, 0), 0.000421, True),
      ((2, 1, -1, 1, 0, 0, 0, 0), -0.000366, True),
      ((2, 1, 0, 1, 0, 0, 0, 0), -0.000351, True),
      ((4, 0, 0, 1, 0, 0, 0, 0), 0.000331, True),
      ((2, -1, 1, 1, 0, 0, 0, 0), 0.000315, True),
      ((2, -2, 0, -1, 0, 0, 0, 0), 0.000302, True),
      ((0, 0, 1, 3, 0, 0, 0, 0), -0.000283, True),
      ((2, 1, 1, -1, 0, 0, 0, 0), -0.000229, True),
      ((1, 1, 0, -1, 0, 0, 0, 0), 0.000223, True),
      ((1, 1, 0, 1, 0, 0, 0, 0), 0.000223, True),
      ((0, 1, -2, -1, 0, 0, 0, 0), -0.00022, True),
      ((2, 1, -1, -1, 0, 0, 0, 0), -0.00022, True),
      ((1, 0, 1, 1, 0, 0, 0, 0), -0.000185, True),
      ((2, -1, -2, -1, 0, 0, 0, 0), 0.000181, True),
      ((0, 1, 2, 1, 0, 0, 0, 0), -0.000177, True),
      ((4, 0, -2, -1, 0, 0, 0, 0), 0.000176, True),
      ((4, -1, -1, -1, 0, 0, 0, 0), 0.000166, True),
      ((1, 0, 1, -1, 0, 0, 0, 0), -0.000164, True),
      ((4, 0, 1, -1, 0, 0, 0, 0), 0.000132, True),
      ((1, 0, -1, -1, 0, 0, 0, 0), -0.000119, True),
      ((4, -1, 0, -1, 0, 0, 0, 0), 0.000115, True),
      ((2, -2, 0, 1, 0, 0, 0, 0), 0.000107, True),
      ((0, 0, 0, 0, 0, 0, 0, 1), -0.002235, True),
      ((0, 0, 0, 0, 0, 0, 1, 0), 0.000382, True),
      ((0, 0, 0, -1, 1, 0, 0, 0), 0.000175, True),
      ((0, 0, 0, 1, 1, 0, 0, 0), 0.000175, True),
      ((0, 0, -1, 0, 0, 0, 0, 1), 0.000127, True),
      ((0, 0, 1, 0, 0, 0, 0, 1), -0.000115, True));
   Dist : constant Terms := (
      ((0, 0, 1, 0, 0, 0, 0, 0), -20905.355, False),
      ((2, 0, -1, 0, 0, 0, 0, 0), -3699.111, False),
      ((2, 0, 0, 0, 0, 0, 0, 0), -2955.968, False),
      ((0, 0, 2, 0, 0, 0, 0, 0), -569.925, False),
      ((0, 1, 0, 0, 0, 0, 0, 0), 48.888, False),
      ((0, 0, 0, 2, 0, 0, 0, 0), -3.149, False),
      ((2, 0, -2, 0, 0, 0, 0, 0), 246.158, False),
      ((2, -1, -1, 0, 0, 0, 0, 0), -152.138, False),
      ((2, 0, 1, 0, 0, 0, 0, 0), -170.733, False),
      ((2, -1, 0, 0, 0, 0, 0, 0), -204.586, False),
      ((0, 1, -1, 0, 0, 0, 0, 0), -129.62, False),
      ((1, 0, 0, 0, 0, 0, 0, 0), 108.743, False),
      ((0, 1, 1, 0, 0, 0, 0, 0), 104.755, False),
      ((2, 0, 0, -2, 0, 0, 0, 0), 10.321, False),
      ((0, 0, 1, -2, 0, 0, 0, 0), 79.661, False),
      ((4, 0, -1, 0, 0, 0, 0, 0), -34.782, False),
      ((0, 0, 3, 0, 0, 0, 0, 0), -23.21, False),
      ((4, 0, -2, 0, 0, 0, 0, 0), -21.636, False),
      ((2, 1, -1, 0, 0, 0, 0, 0), 24.208, False),
      ((2, 1, 0, 0, 0, 0, 0, 0), 30.824, False),
      ((1, 0, -1, 0, 0, 0, 0, 0), -8.379, False),
      ((1, 1, 0, 0, 0, 0, 0, 0), -16.675, False),
      ((2, -1, 1, 0, 0, 0, 0, 0), -12.831, False),
      ((2, 0, 2, 0, 0, 0, 0, 0), -10.445, False),
      ((4, 0, 0, 0, 0, 0, 0, 0), -11.65, False),
      ((2, 0, -3, 0, 0, 0, 0, 0), 14.403, False),
      ((0, 1, -2, 0, 0, 0, 0, 0), -7.003, False),
      ((2, -1, -2, 0, 0, 0, 0, 0), 10.056, False),
      ((1, 0, 1, 0, 0, 0, 0, 0), 6.322, False),
      ((2, -2, 0, 0, 0, 0, 0, 0), -9.884, False),
      ((0, 1, 2, 0, 0, 0, 0, 0), 5.751, False),
      ((2, -2, -1, 0, 0, 0, 0, 0), -4.95, False),
      ((2, 0, 1, -2, 0, 0, 0, 0), 4.13, False),
      ((4, -1, -1, 0, 0, 0, 0, 0), -3.958, False),
      ((3, 0, -1, 0, 0, 0, 0, 0), 3.258, False),
      ((2, 1, 1, 0, 0, 0, 0, 0), 2.616, False),
      ((4, -1, -2, 0, 0, 0, 0, 0), -1.897, False),
      ((0, 2, -1, 0, 0, 0, 0, 0), -2.117, False),
      ((2, 2, -1, 0, 0, 0, 0, 0), 2.354, False),
      ((4, 0, 1, 0, 0, 0, 0, 0), -1.423, False),
      ((0, 0, 4, 0, 0, 0, 0, 0), -1.117, False),
      ((4, -1, 0, 0, 0, 0, 0, 0), -1.571, False),
      ((1, 0, -2, 0, 0, 0, 0, 0), -1.739, False),
      ((0, 0, 2, -2, 0, 0, 0, 0), -4.421, False),
      ((0, 2, 1, 0, 0, 0, 0, 0), 1.165, False),
      ((2, 0, -1, -2, 0, 0, 0, 0), 8.752, False));

   function On return Boolean is (Zonal_On);
   function Hybrid return Boolean is (Hybrid_On);
   function Blend return Boolean is (Blend_On);
   --  0.0669 = variance of the daily forcing over the dLOD record at the default Kappa
   function Forcing_Variance return Long_Float is (0.0669 * (Kappa / (-5.866344))**2);
   function JD_Of (T : Long_Float; Year_Len : Long_Float := 0.0) return Long_Float is
     (JD_Anchor + (T - T_Anchor) * (if Year_Len > 0.0 then Year_Len else Year_Days));

   function Is_Arg (X : Term; A, B, C, D : Integer) return Boolean is
     (X.Mu (1) = A and X.Mu (2) = B and X.Mu (3) = C and X.Mu (4) = D
      and X.Mu (5) = 0 and X.Mu (6) = 0 and X.Mu (7) = 0 and X.Mu (8) = 0);
   --  tier-1 multiplier for a longitude/distance term
   function LR_Mult (X : Term) return Long_Float is
     (if Is_Arg (X, 0, 0, 1, 0) then M_Ecc elsif Is_Arg (X, 2, 0, -1, 0) then M_Evect
      elsif Is_Arg (X, 2, 0, 0, 0) then M_Var elsif Is_Arg (X, 0, 1, 0, 0) then M_AnnEq else 1.0);

   function Series (S : Terms; V : Vals; E : Long_Float; Latitude : Boolean) return Long_Float is
      Acc, Ang, Ef : Long_Float := 0.0;
   begin
      for X of S loop
         Ang := 0.0;
         for K in 1 .. 8 loop
            Ang := Ang + Long_Float (X.Mu (K)) * V (K);
         end loop;
         Ef := E ** (abs X.Mu (2)) * (if Latitude then M_Incl else LR_Mult (X));
         Acc := Acc + Ef * X.C * (if X.Is_Sin then Sin (Ang) else Cos (Ang));
      end loop;
      return Acc;
   end Series;

   function Potential (JD : Long_Float) return Long_Float is
      T  : constant Long_Float := (JD - 2451545.0) / 36525.0;
      L  : constant Long_Float := 218.3164477 + (481267.88123421 + (-0.0015786 + (1.85584E-6 - 1.53388E-8 * T) * T) * T) * T;
      M  : constant Long_Float := 357.5291092 + (35999.0502909 + (-0.0001536 - 4.0833E-8 * T) * T) * T;
      Mm : constant Long_Float := 134.9633964 + (477198.8675055 + (0.0087414 + (1.43741E-6 - 6.7972E-8 * T) * T) * T) * T;
      DD : constant Long_Float := 297.8501921 + (445267.1114034 + (-1.8819E-3 + (1.831945E-6 - 8.84447E-9 * T) * T) * T) * T;
      FF : constant Long_Float := 93.272095 + (483202.0175233 + (-0.0036539 + (-2.8361E-7 + 1.15833E-9 * T) * T) * T) * T;
      A1 : constant Long_Float := 119.75 + 131.849 * T;
      A2 : constant Long_Float := 53.09 + 479264.290 * T;
      A3 : constant Long_Float := 313.45 + 481266.484 * T;
      Eps : constant Long_Float := (23.43929111 + (-0.0130041667 + (-1.638889E-7 + 5.036111E-7 * T) * T) * T) * D2R;
      E  : constant Long_Float := 1.0 - (0.002516 + 0.0000074 * T) * T;
      V  : constant Vals := (DD * D2R, M * D2R, Mm * D2R, FF * D2R, A1 * D2R, A2 * D2R, A3 * D2R, L * D2R);
      Lam : constant Long_Float := (L + Series (Lon, V, E, False)) * D2R;
      Bet : constant Long_Float := Series (Lat, V, E, True) * D2R;
      Rm  : constant Long_Float := (385000.56 + Series (Dist, V, E, False)) / 385000.56;
      SDec : constant Long_Float := Sin (Bet) * Cos (Eps) + Cos (Bet) * Sin (Eps) * Sin (Lam);
      Res : Long_Float := Rm ** (-3) * (3.0 * SDec * SDec - 1.0) / 2.0;
      --  Sun (Keplerian, equation of the centre)
      L0 : constant Long_Float := 280.46646 + 36000.76983 * T;
      Mr : constant Long_Float := M * D2R;
      C  : constant Long_Float := (1.914602 - 0.004817 * T) * Sin (Mr) + 0.019993 * Sin (2.0 * Mr) + 0.000289 * Sin (3.0 * Mr);
      Ec : constant Long_Float := 0.016708634 - 0.000042037 * T;
      Nu : constant Long_Float := Mr + C * D2R;
      Rs : constant Long_Float := (1.000001018 * (1.0 - Ec * Ec)) / (1.0 + Ec * Cos (Nu));
      Sd : constant Long_Float := Sin (Eps) * Sin ((L0 + C) * D2R);
   begin
      Res := Res + W_Sun * Rs ** (-3) * (3.0 * Sd * Sd - 1.0) / 2.0;
      return Res;
   end Potential;

   function Rate (JD : Long_Float) return Long_Float is
     ((Potential (JD + 0.5) - Potential (JD - 0.5)) / 1.0);

   --  Integ = the jerk factor (shfT): forcing + (Integ / 2 pi) * d(forcing)/dt,
   --  t in years -- the same term Tide_Sum_Diff adds per line
   --  (Integ * Freq * Amplitude in quadrature).
   --  mean of V over 1880-2021 (set at elaboration), removed from the LOD term
   V_Mean : Long_Float := 0.0;
   function Forcing_At
     (T : Long_Float; Integ : Long_Float := 0.0; Year_Len : Long_Float := 0.0;
      Lod : Long_Float := 0.0) return Long_Float is
      YD : constant Long_Float := (if Year_Len > 0.0 then Year_Len else Year_Days);
      JD : constant Long_Float := JD_Of (T, YD) - Tau;
      F  : Long_Float := Kappa * Rate (JD);
   begin
      if Integ /= 0.0 then
         F := F + Integ / (2.0 * Ada.Numerics.Pi) * Kappa
           * (Rate (JD + 0.25) - Rate (JD - 0.25)) / 0.5 * YD;
      end if;
      if Lod /= 0.0 then
         F := F + Lod * 2.0 * Ada.Numerics.Pi * Kappa * (Potential (JD) - V_Mean) / YD;
      end if;
      return F;
   end Forcing_At;

   type Vec_Access is access LF_Vec;
   procedure Free is new Ada.Unchecked_Deallocation (LF_Vec, Vec_Access);
   type Slot is record
      N : Integer := 0;
      YD, T_First, T_Mid, T_Last : Long_Float := 0.0;
      R0, DR, PV : Vec_Access;      -- Rate, Rate (+0.25) - Rate (-0.25), Potential - V_Mean
      Stamp : Natural := 0;
   end record;
   type Slot_Array is array (1 .. 2) of Slot;
   type Cache is record
      S : Slot_Array;
      Clock : Natural := 0;
   end record;
   type Cache_Access is access Cache;
   package Cache_Attr is new Ada.Task_Attributes (Cache_Access, null);

   function Forcing_Series
     (T : LF_Vec; Integ : Long_Float := 0.0; Year_Len : Long_Float := 0.0;
      Lod : Long_Float := 0.0) return LF_Vec
   is
      YD : constant Long_Float := (if Year_Len > 0.0 then Year_Len else Year_Days);
      Res : LF_Vec (T'Range);
      C : Cache_Access := Cache_Attr.Value;
      K : Integer := 0;
      Mid : constant Integer := (if T'Length = 0 then T'First else T'First + T'Length / 2);
   begin
      if T'Length = 0 then
         return Res;
      end if;
      if C = null then
         C := new Cache;
         Cache_Attr.Set_Value (C);
      end if;
      for I in C.S'Range loop
         if C.S (I).N = T'Length and then C.S (I).YD = YD
           and then C.S (I).T_First = T (T'First) and then C.S (I).T_Mid = T (Mid)
           and then C.S (I).T_Last = T (T'Last)
         then
            K := I;
         end if;
      end loop;
      if K = 0 then
         K := (if C.S (1).Stamp <= C.S (2).Stamp then 1 else 2);   -- least recently used
         Free (C.S (K).R0); Free (C.S (K).DR); Free (C.S (K).PV);
         C.S (K).N := T'Length; C.S (K).YD := YD;
         C.S (K).T_First := T (T'First); C.S (K).T_Mid := T (Mid); C.S (K).T_Last := T (T'Last);
         C.S (K).R0 := new LF_Vec (T'Range);
         for I in T'Range loop
            C.S (K).R0 (I) := Rate (JD_Of (T (I), YD) - Tau);
         end loop;
      end if;
      C.Clock := C.Clock + 1;
      C.S (K).Stamp := C.Clock;
      if Integ /= 0.0 and then C.S (K).DR = null then
         C.S (K).DR := new LF_Vec (T'Range);
         for I in T'Range loop
            C.S (K).DR (I) := Rate (JD_Of (T (I), YD) - Tau + 0.25) - Rate (JD_Of (T (I), YD) - Tau - 0.25);
         end loop;
      end if;
      if Lod /= 0.0 and then C.S (K).PV = null then
         C.S (K).PV := new LF_Vec (T'Range);
         for I in T'Range loop
            C.S (K).PV (I) := Potential (JD_Of (T (I), YD) - Tau) - V_Mean;
         end loop;
      end if;
      --  same expressions, in the same order, as Forcing_At
      for I in T'Range loop
         declare
            F : Long_Float := Kappa * C.S (K).R0 (C.S (K).R0'First + (I - T'First));
         begin
            if Integ /= 0.0 then
               F := F + Integ / (2.0 * Ada.Numerics.Pi) * Kappa
                 * C.S (K).DR (C.S (K).DR'First + (I - T'First)) / 0.5 * YD;
            end if;
            if Lod /= 0.0 then
               F := F + Lod * 2.0 * Ada.Numerics.Pi * Kappa * C.S (K).PV (C.S (K).PV'First + (I - T'First)) / YD;
            end if;
            Res (I) := F;
         end;
      end loop;
      return Res;
   end Forcing_Series;
begin
   declare
      N : constant Integer := 17_166;   -- 1880-01-01 to 2021, every 3 days
      S : Long_Float := 0.0;
   begin
      for I in 0 .. N loop
         S := S + Potential (2_407_715.5 + 3.0 * Long_Float (I));
      end loop;
      V_Mean := S / Long_Float (N + 1);
   end;
end GEM.Zonal;
