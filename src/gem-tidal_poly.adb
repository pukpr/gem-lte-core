with Ada.Numerics;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Long_Elementary_Functions;
with Ada.Text_IO;
with Ada.Long_Float_Text_IO;
with GEM.LTE.Primitives;

package body GEM.Tidal_Poly is

   package RA is new Ada.Numerics.Generic_Real_Arrays (Real => Long_Float);
   use RA;
   use Ada.Numerics.Long_Elementary_Functions;

   Pi     : constant Long_Float := Ada.Numerics.Pi;
   Two_Pi : constant Long_Float := 2.0 * Pi;
   D2R    : constant Long_Float := Pi / 180.0;

   --  Julian date anchor shared with GEM.Zonal and dlod3.dat
   JD_Anchor   : constant Long_Float := 2_448_074.761_995;
   Year_Anchor : constant Long_Float := 1990.5;
   JD_J2000    : constant Long_Float := 2_451_545.0;

   --  IAU / Chapront linear mean orbital arguments at J2000.0 (degrees)
   --  and linear rates (degrees per Julian century = 36525 days).
   --  ZERO Meeus perturbation tables!
   L0_Deg  : constant Long_Float := 218.316_447_7;
   L1_Deg  : constant Long_Float := 481_267.881_234_21;
   Mp0_Deg : constant Long_Float := 134.963_396_4;
   Mp1_Deg : constant Long_Float := 477_198.867_505_5;
   F0_Deg  : constant Long_Float := 93.272_095_0;
   F1_Deg  : constant Long_Float := 483_202.017_523_3;
   D0_Deg  : constant Long_Float := 297.850_192_1;
   D1_Deg  : constant Long_Float := 445_267.111_403_4;

   --  Derived linear rates in cycles/day
   R_L  : constant Long_Float := L1_Deg / (360.0 * 36_525.0);
   R_Mp : constant Long_Float := Mp1_Deg / (360.0 * 36_525.0);
   R_F  : constant Long_Float := F1_Deg / (360.0 * 36_525.0);
   R_D  : constant Long_Float := D1_Deg / (360.0 * 36_525.0);
   R_Om : constant Long_Float := R_L - R_F;  --  Omega = L - F (retrograde)

   function Upper (S : String) return String is
      R : String := S;
   begin
      for I in R'Range loop
         if R (I) in 'a' .. 'z' then
            R (I) := Character'Val (Character'Pos (R (I)) - 32);
         end if;
      end loop;
      return R;
   end Upper;

   Poly_On : constant Boolean := Upper (GEM.Getenv ("TIDES", "TABLE")) = "POLY";
   Tier_Raw : constant String := Upper (GEM.Getenv ("POLY_TIER", "Z10"));
   Lod_Floor : constant Long_Float := GEM.Getenv ("LOD_FLOOR", 0.0);

   type Tier_Kind is
     (Tier_Z9, Tier_Z10, Tier_Z12, Tier_Z15, Tier_Z10_D6, Tier_Z15_D6,
      Tier_Z10_N5, Tier_Z10_D6_N5, Tier_Z10_S5_N5, Tier_Z12_S5_N5,
      Tier_Z10_S5_D6_N5, Tier_Z12_S5_D6_N5);

   function Resolve_Tier return Tier_Kind is
   begin
      if Tier_Raw = "Z9" or else Tier_Raw = "Z18" then
         return Tier_Z9;
      elsif Tier_Raw = "Z10" or else Tier_Raw = "Z20" then
         return Tier_Z10;
      elsif Tier_Raw = "Z12" then
         return Tier_Z12;
      elsif Tier_Raw = "Z15" or else Tier_Raw = "Z30" then
         return Tier_Z15;
      elsif Tier_Raw = "Z10_D6" or else Tier_Raw = "F32" then
         return Tier_Z10_D6;
      elsif Tier_Raw = "Z15_D6" then
         return Tier_Z15_D6;
      elsif Tier_Raw = "Z10_N5" then
         return Tier_Z10_N5;
      elsif Tier_Raw = "Z10_D6_N5" then
         return Tier_Z10_D6_N5;
      elsif Tier_Raw = "Z10_S5_N5" then
         return Tier_Z10_S5_N5;
      elsif Tier_Raw = "Z12_S5_N5" then
         return Tier_Z12_S5_N5;
      elsif Tier_Raw = "Z10_S5_D6_N5" then
         return Tier_Z10_S5_D6_N5;
      elsif Tier_Raw = "Z12_S5_D6_N5" then
         return Tier_Z12_S5_D6_N5;
      else
         return Tier_Z10;
      end if;
   end Resolve_Tier;

   Active_Tier : constant Tier_Kind := Resolve_Tier;

   --  Each gravitational polynomial term is an integer combination
   --  n_L * L + n_Mp * M' + n_Om * Omega + n_D * D + n_Day * theta_d,
   --  evaluated algebraically via phasor products of the 5 elementary angles.
   type Term_Spec is record
      N_L   : Integer;
      N_Mp  : Integer;
      N_Om  : Integer;
      N_D   : Integer;
      N_Day : Integer;
      Label : String (1 .. 10);
   end record;

   --  Master ordered catalog of 30 Meeus-free gravitational monomials:
   --   1..9   : Z9  (core zonal: Mf, Mf_nodal, 2F, Mm, MSm, Mt, Mt_nodal, MSf, Ssa)
   --   10     : Z10 (+ 2D+M' [9.6137d])
   --   11..15 : Z15 (+ L, F, M'-Om, 2L-M'-Om, Mtm)
   --   16..21 : D6  (K1, O1, P1, M2, N2, K2)
   --   22..26 : N5  (multi-year nodal/perigee: Om 18.61y, 2Om 9.31y, p 8.85y, 2p 4.43y, p+N 6.00y)
   --   27..30 : S4  (nodal sidebands of MSm, MSf, Mtm identified by 30-yr Baltic diagnostic)
   Catalog : constant array (1 .. Max_Poly_Terms) of Term_Spec :=
     (1  => ( 2,  0,  0,  0, 0, "Mf(2L)    "),  -- 13.6608 d
      2  => ( 2,  0, -1,  0, 0, "Mf_nod    "),  -- 13.6334 d
      3  => ( 2,  0, -2,  0, 0, "2F        "),  -- 13.6061 d
      4  => ( 0,  1,  0,  0, 0, "Mm(Mp)    "),  -- 27.5545 d
      5  => ( 0, -1,  0,  2, 0, "MSm(ev)   "),  -- 31.8119 d
      6  => ( 2,  1,  0,  0, 0, "Mt        "),  --  9.1329 d
      7  => ( 2,  1, -1,  0, 0, "Mt_nod    "),  --  9.1207 d
      8  => ( 0,  0,  0,  2, 0, "MSf(2D)   "),  -- 14.7653 d
      9  => ( 2,  0,  0, -2, 0, "Ssa(2Ls)  "),  -- 182.621 d
      10 => ( 0,  1,  0,  2, 0, "2D+Mp     "),  --  9.6137 d
      11 => ( 1,  0,  0,  0, 0, "T(L)      "),  -- 27.3216 d
      12 => ( 1,  0, -1,  0, 0, "D(F)      "),  -- 27.2122 d
      13 => ( 0,  1, -1,  0, 0, "Mp-Om     "),  -- 27.4432 d
      14 => ( 2, -1, -1,  0, 0, "2L-Mp-Om  "),  -- 26.9850 d
      15 => ( 2, -1,  0,  2, 0, "Mtm       "),  --  9.5569 d
      16 => ( 1,  0,  0, -1, 1, "K1        "),  --  0.9973 d
      17 => (-1,  0,  0, -1, 1, "O1        "),  --  1.0758 d
      18 => (-1,  0,  0,  1, 1, "P1        "),  --  1.0027 d
      19 => ( 0,  0,  0, -2, 2, "M2        "),  --  0.5175 d
      20 => ( 0, -1,  0, -2, 2, "N2        "),  --  0.5274 d
      21 => ( 2,  0,  0, -2, 2, "K2        "),  --  0.4986 d
      22 => ( 0,  0, -1,  0, 0, "N(18.6y)  "),  -- 6798.38 d
      23 => ( 0,  0, -2,  0, 0, "2N(9.3y)  "),  -- 3399.19 d
      24 => ( 1, -1,  0,  0, 0, "p(8.85y)  "),  -- 3231.50 d
      25 => ( 2, -2,  0,  0, 0, "2p(4.43y) "),  -- 1615.75 d
      26 => ( 1, -1, -1,  0, 0, "p+N(6.0y) "),  -- 2190.24 d
      27 => ( 0, -1, -1,  2, 0, "MSm_nod   "),  -- 31.6639 d
      28 => ( 0,  0, -1,  2, 0, "MSf_nod+  "),  -- 14.7333 d
      29 => ( 0,  0,  1,  2, 0, "MSf_nod-  "),  -- 14.7975 d
      30 => ( 2, -1, -1,  2, 0, "Mtm_nod   "),  --  9.5435 d
      31 => ( 2,  2,  0,  0, 0, "Mqm(6.86) "),  --  6.8594 d
      32 => ( 2,  0,  0,  2, 0, "Msqm(7.1) ")); --  7.0958 d

   --  Map active slot K (1 .. Active_Count) to its Catalog index
   function Catalog_Index (K : Positive) return Positive is
   begin
      case Active_Tier is
         when Tier_Z9     => return K;
         when Tier_Z10    => return K;
         when Tier_Z12    =>
            if K <= 10 then
               return K;
            else
               return K + 20; --  11..12 -> 31..32 (Mqm, Msqm)
            end if;
         when Tier_Z15    => return K;
         when Tier_Z10_D6 =>
            if K <= 10 then
               return K;
            else
               return K + 5;  --  11..16 -> 16..21 (D6)
            end if;
         when Tier_Z15_D6 => return K;
         when Tier_Z10_N5 =>
            if K <= 10 then
               return K;
            else
               return K + 11; --  11..15 -> 22..26 (N5)
            end if;
         when Tier_Z10_D6_N5 =>
            if K <= 10 then
               return K;
            else
               return K + 5;  --  11..16 -> 16..21 (D6), 17..21 -> 22..26 (N5)
            end if;
         when Tier_Z10_S5_N5 =>
            if K <= 10 then
               return K;
            elsif K = 11 then
               return 15;     --  Mtm (9.557d)
            elsif K <= 15 then
               return K + 15; --  12..15 -> 27..30 (S4)
            else
               return K + 6;  --  16..20 -> 22..26 (N5)
            end if;
         when Tier_Z12_S5_N5 =>
            if K <= 10 then
               return K;
            elsif K <= 12 then
               return K + 20; --  11..12 -> 31..32 (Mqm, Msqm)
            elsif K = 13 then
               return 15;     --  Mtm (9.557d)
            elsif K <= 17 then
               return K + 13; --  14..17 -> 27..30 (S4)
            else
               return K + 4;  --  18..22 -> 22..26 (N5)
            end if;
         when Tier_Z10_S5_D6_N5 =>
            if K <= 10 then
               return K;
            elsif K = 11 then
               return 15;     --  Mtm (9.557d)
            elsif K <= 15 then
               return K + 15; --  12..15 -> 27..30 (S4)
            else
               return K;      --  16..21 -> 16..21 (D6), 22..26 -> 22..26 (N5)
            end if;
         when Tier_Z12_S5_D6_N5 =>
            if K <= 10 then
               return K;
            elsif K <= 12 then
               return K + 20; --  11..12 -> 31..32 (Mqm, Msqm)
            elsif K = 13 then
               return 15;     --  Mtm (9.557d)
            elsif K <= 17 then
               return K + 13; --  14..17 -> 27..30 (S4)
            else
               return K - 2;  --  18..23 -> 16..21 (D6), 24..28 -> 22..26 (N5)
            end if;
      end case;
   end Catalog_Index;

   function On return Boolean is (Poly_On);

   function Tier_Name return String is
   begin
      case Active_Tier is
         when Tier_Z9           => return "Z9";
         when Tier_Z10          => return "Z10";
         when Tier_Z12          => return "Z12";
         when Tier_Z15          => return "Z15";
         when Tier_Z10_D6       => return "Z10_D6";
         when Tier_Z15_D6       => return "Z15_D6";
         when Tier_Z10_N5       => return "Z10_N5";
         when Tier_Z10_D6_N5    => return "Z10_D6_N5";
         when Tier_Z10_S5_N5    => return "Z10_S5_N5";
         when Tier_Z12_S5_N5    => return "Z12_S5_N5";
         when Tier_Z10_S5_D6_N5 => return "Z10_S5_D6_N5";
         when Tier_Z12_S5_D6_N5 => return "Z12_S5_D6_N5";
      end case;
   end Tier_Name;

   function Active_Count return Positive is
   begin
      case Active_Tier is
         when Tier_Z9           => return 9;
         when Tier_Z10          => return 10;
         when Tier_Z12          => return 12;
         when Tier_Z15          => return 15;
         when Tier_Z10_D6       => return 16;
         when Tier_Z15_D6       => return 21;
         when Tier_Z10_N5       => return 15;
         when Tier_Z10_D6_N5    => return 21;
         when Tier_Z10_S5_N5    => return 20;
         when Tier_Z12_S5_N5    => return 22;
         when Tier_Z10_S5_D6_N5 => return 26;
         when Tier_Z12_S5_D6_N5 => return 28;
      end case;
   end Active_Count;

   function Is_Long_Term (K : Positive) return Boolean is
   begin
      return abs Term_Period_Days (K) > 365.25;
   end Is_Long_Term;

   function Term_Name (K : Positive) return String is
   begin
      return Catalog (Catalog_Index (K)).Label;
   end Term_Name;

   function Rate_CPD (Idx : Positive) return Long_Float is
      S : Term_Spec renames Catalog (Idx);
   begin
      return Long_Float (S.N_L) * R_L
        + Long_Float (S.N_Mp) * R_Mp
        + Long_Float (S.N_Om) * R_Om
        + Long_Float (S.N_D) * R_D
        + Long_Float (S.N_Day);
   end Rate_CPD;

   function Term_Period_Days (K : Positive) return Long_Float is
      R : constant Long_Float := Rate_CPD (Catalog_Index (K));
   begin
      if abs R > 1.0e-15 then
         return 1.0 / R;
      else
         return 0.0;
      end if;
   end Term_Period_Days;

   function Term_Freq_CPY
     (K : Positive; Year_Len : Long_Float := GEM.LTE.Year_Length)
      return Long_Float
   is
   begin
      return Year_Len * Rate_CPD (Catalog_Index (K));
   end Term_Freq_CPY;

   --  Complex phasor helpers (pure algebraic products, no extra trig calls!)
   type Phasor is record
      C, S : Long_Float;
   end record;

   function Mul (A, B : Phasor) return Phasor is
     ((C => A.C * B.C - A.S * B.S,
       S => A.S * B.C + A.C * B.S));

   function Conj (A : Phasor) return Phasor is
     ((C => A.C, S => -A.S));

   function Pow (A : Phasor; N : Integer) return Phasor is
   begin
      case N is
         when 0  => return (1.0, 0.0);
         when 1  => return A;
         when -1 => return Conj (A);
         when 2  => return (A.C * A.C - A.S * A.S, 2.0 * A.C * A.S);
         when -2 => return (A.C * A.C - A.S * A.S, -2.0 * A.C * A.S);
         when others =>
            declare
               Base : constant Phasor := (if N > 0 then A else Conj (A));
               R    : Phasor := (1.0, 0.0);
            begin
               for I in 1 .. abs N loop
                  R := Mul (R, Base);
               end loop;
               return R;
            end;
      end case;
   end Pow;

   type Basis_Vec is array (1 .. Max_Poly_Terms) of Phasor;

   procedure Compute_Basis
     (Date_Yr  : in Long_Float;
      Year_Len : in Long_Float;
      N_Act    : in Positive;
      Need_Day : in Boolean;
      B        : out Basis_Vec)
   is
      JD : constant Long_Float :=
        JD_Anchor + (Date_Yr - Year_Anchor) * Year_Len;
      T  : constant Long_Float := (JD - JD_J2000) / 36_525.0;

      --  4 fundamental linear mean angles (radians)
      L_Rad  : constant Long_Float := (L0_Deg  + L1_Deg  * T) * D2R;
      Mp_Rad : constant Long_Float := (Mp0_Deg + Mp1_Deg * T) * D2R;
      F_Rad  : constant Long_Float := (F0_Deg  + F1_Deg  * T) * D2R;
      D_Rad  : constant Long_Float := (D0_Deg  + D1_Deg  * T) * D2R;
      Om_Rad : constant Long_Float := L_Rad - F_Rad;

      P_L  : constant Phasor := (Cos (L_Rad),  Sin (L_Rad));
      P_Mp : constant Phasor := (Cos (Mp_Rad), Sin (Mp_Rad));
      P_Om : constant Phasor := (Cos (Om_Rad), Sin (Om_Rad));
      P_D  : constant Phasor := (Cos (D_Rad),  Sin (D_Rad));
      P_2L : constant Phasor := (P_L.C * P_L.C - P_L.S * P_L.S,
                                 2.0 * P_L.C * P_L.S);
      P_2D : constant Phasor := (P_D.C * P_D.C - P_D.S * P_D.S,
                                 2.0 * P_D.C * P_D.S);
      P_Day  : Phasor := (1.0, 0.0);
      P_2Day : Phasor := (1.0, 0.0);
   begin
      if Need_Day then
         declare
            Day_Frac : constant Long_Float :=
              (JD - JD_J2000) - Long_Float'Floor (JD - JD_J2000);
            Th_Day : constant Long_Float := Two_Pi * Day_Frac;
         begin
            P_Day  := (Cos (Th_Day), Sin (Th_Day));
            P_2Day := (P_Day.C * P_Day.C - P_Day.S * P_Day.S,
                       2.0 * P_Day.C * P_Day.S);
         end;
      end if;

      for K in 1 .. N_Act loop
         declare
            S : Term_Spec renames Catalog (Catalog_Index (K));
            Z : Phasor := (1.0, 0.0);
         begin
            if S.N_L = 2 then
               Z := P_2L;
            elsif S.N_L /= 0 then
               Z := Pow (P_L, S.N_L);
            end if;

            if S.N_Mp /= 0 then
               Z := Mul (Z, Pow (P_Mp, S.N_Mp));
            end if;

            if S.N_Om /= 0 then
               Z := Mul (Z, Pow (P_Om, S.N_Om));
            end if;

            if S.N_D = 2 then
               Z := Mul (Z, P_2D);
            elsif S.N_D = -2 then
               Z := Mul (Z, Conj (P_2D));
            elsif S.N_D /= 0 then
               Z := Mul (Z, Pow (P_D, S.N_D));
            end if;

            if S.N_Day = 1 then
               Z := Mul (Z, P_Day);
            elsif S.N_Day = 2 then
               Z := Mul (Z, P_2Day);
            elsif S.N_Day /= 0 then
               Z := Mul (Z, Pow (P_Day, S.N_Day));
            end if;

            B (K) := Z;
         end;
      end loop;
   end Compute_Basis;

   function Lod_Factor (Freq : Long_Float; Lod : Long_Float) return Long_Float is
   begin
      if Lod = 0.0 or else Freq = 0.0 then
         return 0.0;
      elsif Lod_Floor = 0.0 then
         return Lod / Freq;
      else
         return Lod / (Freq * Sqrt (1.0 + (Lod_Floor / Freq)**2));
      end if;
   end Lod_Factor;

   function Forcing_Series
     (T            : LF_Vec;
      Constituents : GEM.LTE.Long_Periods_Amp_Phase;
      Year_Len     : Long_Float := GEM.LTE.Year_Length;
      Integ        : Long_Float := 0.0;
      Lod          : Long_Float := 0.0) return LF_Vec
   is
      N_Act    : constant Positive := Active_Count;
      Need_Day : constant Boolean :=
        Active_Tier = Tier_Z10_D6 or else Active_Tier = Tier_Z15_D6
        or else Active_Tier = Tier_Z10_D6_N5
        or else Active_Tier = Tier_Z10_S5_D6_N5
        or else Active_Tier = Tier_Z12_S5_D6_N5;
      Res      : LF_Vec (T'Range);
      First_C  : constant Integer := Constituents'First;

      --  Precompute linear sin/cos coefficients and jerk/LOD weights per term:
      --    Amp * sin(theta + phase) = (Amp*cos(phase))*sin(theta) + (Amp*sin(phase))*cos(theta)
      --    d/dt sin(theta + phase) / (2*pi) = Freq * cos(theta + phase)
      --      = Freq * [(Amp*cos(phase))*cos(theta) - (Amp*sin(phase))*sin(theta)]
      --    Match Tide_Sum_Diff (Cos_Phase = False):
      --      Val = Amp*sin(theta+ph) + (Integ*Freq - Lod_F)*Amp*cos(theta+ph)
      --          = A_sin * sin(theta) + A_cos * cos(theta)
      A_Sin : array (1 .. Max_Poly_Terms) of Long_Float := (others => 0.0);
      A_Cos : array (1 .. Max_Poly_Terms) of Long_Float := (others => 0.0);
      B     : Basis_Vec;
   begin
      for K in 1 .. N_Act loop
         declare
            AP   : GEM.LTE.Amp_Phase renames Constituents (First_C + K - 1);
            Cp   : constant Long_Float := AP.Amplitude * Cos (AP.Phase);
            Sp   : constant Long_Float := AP.Amplitude * Sin (AP.Phase);
            Freq : constant Long_Float := Term_Freq_CPY (K, Year_Len);
            Quad : constant Long_Float := Integ * Freq - Lod_Factor (Freq, Lod);
         begin
            --  Amp*sin(th+ph) = Cp*sin(th) + Sp*cos(th)
            --  Quad*Amp*cos(th+ph) = Quad*(Cp*cos(th) - Sp*sin(th))
            A_Sin (K) := Cp - Quad * Sp;
            A_Cos (K) := Sp + Quad * Cp;
         end;
      end loop;

      for I in T'Range loop
         Compute_Basis (T (I), Year_Len, N_Act, Need_Day, B);
         declare
            Sum : Long_Float := 0.0;
         begin
            for K in 1 .. N_Act loop
               Sum := Sum + A_Sin (K) * B (K).S + A_Cos (K) * B (K).C;
            end loop;
            Res (I) := Sum;
         end;
      end loop;
      return Res;
   end Forcing_Series;

   function Calibrate_dLOD
     (File_Name       : String;
      Year_Correction : Long_Float := 0.0)
      return GEM.LTE.Long_Periods_Amp_Phase
   is
      use GEM.LTE.Primitives;
      D        : constant Data_Pairs := Make_Data (File_Name);
      N_Obs    : constant Integer := D'Length;
      N_Act    : constant Positive := Active_Count;
      Need_Day : constant Boolean :=
        Active_Tier = Tier_Z10_D6 or else Active_Tier = Tier_Z15_D6
        or else Active_Tier = Tier_Z10_D6_N5
        or else Active_Tier = Tier_Z10_S5_D6_N5
        or else Active_Tier = Tier_Z12_S5_D6_N5;
      --  Calibrate only sub-annual terms (|period| <= 1 yr) against dlod3.dat;
      --  multi-year nodal/perigee terms stay at 0.0 (free of sub-annual dLOD prior).
      Sub_Idx  : array (1 .. Max_Poly_Terms) of Positive := (others => 1);
      N_Sub    : Natural := 0;
      Y_Len    : constant Long_Float := GEM.LTE.Year_Length (Year_Correction);
      B        : Basis_Vec;

      AP : GEM.LTE.Long_Periods_Amp_Phase (GEM.LTE.Doodson_Args'Range) :=
        (others => (Amplitude => 0.0, Phase => 0.0));
      Model : Data_Pairs := D;
   begin
      for K in 1 .. N_Act loop
         if not Is_Long_Term (K) then
            N_Sub := N_Sub + 1;
            Sub_Idx (N_Sub) := K;
         end if;
      end loop;

      declare
         N_Cols : constant Positive := 2 * N_Sub + 1;  --  [1, sin_1, cos_1, ...]
         XtX    : Real_Matrix (1 .. N_Cols, 1 .. N_Cols) := (others => (others => 0.0));
         Xty    : Real_Vector (1 .. N_Cols) := (others => 0.0);
         Row    : Real_Vector (1 .. N_Cols);
         Beta   : Real_Vector (1 .. N_Cols);
      begin
         for I in D'Range loop
            Compute_Basis (D (I).Date, Y_Len, N_Act, Need_Day, B);
            Row (1) := 1.0;
            for S in 1 .. N_Sub loop
               Row (2 * S)     := B (Sub_Idx (S)).S;
               Row (2 * S + 1) := B (Sub_Idx (S)).C;
            end loop;
            declare
               Y : constant Long_Float := D (I).Value;
            begin
               for C1 in 1 .. N_Cols loop
                  Xty (C1) := Xty (C1) + Row (C1) * Y;
                  for C2 in C1 .. N_Cols loop
                     XtX (C1, C2) := XtX (C1, C2) + Row (C1) * Row (C2);
                  end loop;
               end loop;
            end;
         end loop;

         for C1 in 1 .. N_Cols loop
            XtX (C1, C1) := XtX (C1, C1) + 1.0e-12;
            for C2 in C1 + 1 .. N_Cols loop
               XtX (C2, C1) := XtX (C1, C2);
            end loop;
         end loop;

         Beta := Solve (XtX, Xty);

         for S in 1 .. N_Sub loop
            declare
               K     : constant Positive := Sub_Idx (S);
               C_Sin : constant Long_Float := Beta (2 * S);
               C_Cos : constant Long_Float := Beta (2 * S + 1);
               Amp   : constant Long_Float := Sqrt (C_Sin * C_Sin + C_Cos * C_Cos);
            begin
               AP (K).Amplitude := Amp;
               AP (K).Phase := (if Amp > 0.0 then Arctan (C_Cos, C_Sin) else 0.0);
            end;
         end loop;

         declare
            T_Vec : LF_Vec (D'Range);
         begin
            for I in D'Range loop
               T_Vec (I) := D (I).Date;
            end loop;
            declare
               F_Vec : constant LF_Vec := Forcing_Series (T_Vec, AP, Y_Len);
            begin
               for I in D'Range loop
                  Model (I).Value := Beta (1) + F_Vec (I);
               end loop;
            end;
         end;
      end;

      Ada.Text_IO.Put ("TIDES=POLY tier=" & Tier_Name & " (" & N_Sub'Img
                       & " sub-annual +" & Natural'Image (N_Act - N_Sub)
                       & " multi-yr, N=" & N_Obs'Img & ") dLOD CC=");
      Ada.Long_Float_Text_IO.Put (CC (D, Model), Fore => 1, Aft => 5, Exp => 0);
      Ada.Text_IO.New_Line;

      return AP;
   end Calibrate_dLOD;

end GEM.Tidal_Poly;
