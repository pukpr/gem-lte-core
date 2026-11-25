with Ada.Calendar;
with Ada.Directories;
with Ada.Long_Float_Text_IO;
with Ada.Numerics;
with Ada.Numerics.Generic_Real_Arrays;
with Ada.Numerics.Long_Elementary_Functions;
with Ada.Strings.Fixed;
with Ada.Text_IO;
with Ada.Unchecked_Deallocation;
with System.Task_Info;

package body GEM.LTE.Primitives.Refine is

   package RA is new Ada.Numerics.Generic_Real_Arrays (Real => Long_Float);
   use RA;
   use Ada.Numerics.Long_Elementary_Functions;
   use type Ada.Calendar.Time;
   use type Shared.Leg_Own_Set;

   Pi : constant Long_Float := Ada.Numerics.Pi;
   Two_Pi : constant Long_Float := 2.0 * Pi;
   Min_LT : constant Long_Float := 0.01;

   function Fmt (X : Long_Float; Aft : Natural := 4) return String is
      Buf : String (1 .. 32);
   begin
      Ada.Long_Float_Text_IO.Put (Buf, X, Aft => Aft, Exp => 0);
      return Ada.Strings.Fixed.Trim (Buf, Ada.Strings.Both);
   end Fmt;

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

   function Leg_P_Mode return Shared.Leg_Own_Set is
      V : constant String := Upper (GEM.Getenv ("LEG_P", "NONE"));
   begin
      if V = "PARAMS" then
         return Shared.Own_Params;
      elsif V = "BESSEL" then
         return Shared.Own_Bessel;
      else
         return (if GEM.Getenv ("LEG_BESSEL", False)
                 then Shared.Own_Bessel else Shared.Own_None);
      end if;
   end Leg_P_Mode;

   --  Convert LPAP (Amp, Phase) <-> linear vector Theta (2*NLP):
   --    Theta(2*J-1) = Amp * Cos(Phase),  Theta(2*J) = Amp * Sin(Phase)
   procedure LPAP_To_Theta
     (AP : in Long_Periods_Amp_Phase; Th : out Real_Vector)
   is
      K : Integer := 0;
   begin
      for J in AP'Range loop
         K := K + 1;
         Th (2 * K - 1) := AP (J).Amplitude * Cos (AP (J).Phase);
         Th (2 * K)     := AP (J).Amplitude * Sin (AP (J).Phase);
      end loop;
   end LPAP_To_Theta;

   procedure Theta_To_LPAP
     (Th : in Real_Vector; AP : in out Long_Periods_Amp_Phase)
   is
      K : Integer := 0;
      C, S : Long_Float;
   begin
      for J in AP'Range loop
         K := K + 1;
         C := Th (2 * K - 1);
         S := Th (2 * K);
         AP (J).Amplitude := Sqrt (C * C + S * S);
         if AP (J).Amplitude > 0.0 then
            AP (J).Phase := Arctan (S, C);
         else
            AP (J).Phase := 0.0;
         end if;
      end loop;
   end Theta_To_LPAP;

   --  Compute pre-Bessel forcing F(t) for a given leg parameter set D_Leg
   --  and LPAP constituents, matching GEM.LTE.Primitives.Solution.Calc_Forcing.
   function Calc_Pre_Bessel
     (Template     : in Data_Pairs;
      LP           : in Long_Periods;
      Constituents : in Long_Periods_Amp_Phase;
      B            : in Shared.Param_B;
      Init_Val     : in Long_Float;
      Strict_IDate : in Boolean;
      IDate        : in Long_Float;
      SPY          : in Long_Float;
      Impulse_Only : in Boolean;
      Imp_Mod12    : in Boolean;
      Use_LagC     : in Boolean := True) return Real_Vector
   is
      SPY_Int : constant Integer := Integer (SPY);

      function Impulse_Delta (Time : Long_Float) return Long_Float is
         Trunc : constant Integer :=
           Integer ((Time - Long_Float'Floor (Time)) * SPY);
         DPos  : constant Integer :=
           (if Imp_Mod12
            then Integer (abs B.DelB * SPY) mod SPY_Int
            else Integer (abs B.DelB * SPY));
      begin
         if Trunc = DPos then
            return B.DelA;
         elsif Trunc = (DPos + SPY_Int / 2) mod SPY_Int then
            return B.Asym;
         else
            return 0.0;
         end if;
      end Impulse_Delta;

      function Impulse_Smear (Time : Long_Float) return Long_Float is
         Trunc : constant Integer :=
           Integer ((Time - Long_Float'Floor (Time)) * SPY);
         DPos  : constant Integer :=
           (if Imp_Mod12
            then Integer (abs B.DelB * SPY) mod SPY_Int
            else Integer (abs B.DelB * SPY));
         DPos2 : constant Integer := (DPos + SPY_Int / 2) mod SPY_Int;
      begin
         if Trunc = DPos2 - 1 then
            return (if B.DelA < 0.0 then -B.DelA else 0.0);
         elsif Trunc = DPos2 then
            return B.Asym;
         elsif Trunc = DPos2 + 1 then
            return (if B.DelA > 0.0 then B.DelA else 0.0);
         else
            return 0.0;
         end if;
      end Impulse_Smear;

      function Impulse_Base (Time : Long_Float) return Long_Float is
      begin
         if Impulse_Only then
            return Impulse_Delta (Time);
         else
            return Impulse_Smear (Time);
         end if;
      end Impulse_Base;

      function Impulse (Time : Long_Float) return Long_Float is
      begin
         if B.delC = 0.0 then
            return Impulse_Base (Time);
         else
            return Impulse_Base (Time)
              - abs B.delC
                * Impulse_Base (Time - Long_Float'Copy_Sign (0.25, B.delC));
         end if;
      end Impulse;

      function Impulse_Amplify is new Amplify (Impulse => Impulse);

      function Shape (Raw : Data_Pairs) return Data_Pairs is
      begin
         if B.mB = 0.0 then
            return Raw;
         end if;
         declare
            R2 : constant Long_Float := Long_Float'Min (abs B.mB, 0.9999);
            G  : Data_Pairs := Raw;
         begin
            G (G'First).Value := (1.0 - R2) * Raw (Raw'First).Value;
            for I in G'First + 1 .. G'Last loop
               G (I).Value := (1.0 - R2) * Raw (I).Value + R2 * G (I - 1).Value;
            end loop;
            return G;
         end;
      end Shape;

      Der   : constant Long_Float := 1.0 - B.mA;
      Lag_C : constant Long_Float := (if Use_LagC then B.mP else 0.0);
      Out_F : Real_Vector (1 .. Template'Length);
   begin
      if Strict_IDate then
         declare
            Ext : constant Data_Pairs :=
              Extend_Backward (Template, IDate, SPY);
            Tide : constant Data_Pairs :=
              Tide_Sum
                (Template => Ext, Constituents => Constituents, Periods => LP,
                 Ref_Time => 0.0, Scaling => 0.0,
                 Year_Len => Year_Length (B.Year), Integ => B.shiftT,
                 Lod => 0.0, Rho => B.rho, ZTau => B.ztau, ZDist => B.zdist,
                 ZTilt => B.ztilt, SsaA => B.ssaA, SsaP => B.ssaP);
            Imp : constant Data_Pairs :=
              Impulse_Amplify
                (Raw => Tide, Offset => 0.0, Ramp => 0.0,
                 Start => Ext (Ext'First).Date);
            Ext_F : constant Data_Pairs :=
              IIR (Raw => Shape (Imp), lagA => Der, lagC => Lag_C,
                   iA => Init_Val, Start => IDate - 0.1 / SPY);
         begin
            for I in Template'Range loop
               Out_F (I - Template'First + 1) := Ext_F (I).Value;
            end loop;
         end;
      else
         declare
            Tide : constant Data_Pairs :=
              Tide_Sum
                (Template => Template, Constituents => Constituents,
                 Periods => LP, Ref_Time => 0.0, Scaling => 0.0,
                 Year_Len => Year_Length (B.Year), Integ => B.shiftT,
                 Lod => 0.0, Rho => B.rho, ZTau => B.ztau, ZDist => B.zdist,
                 ZTilt => B.ztilt, SsaA => B.ssaA, SsaP => B.ssaP);
            Imp : constant Data_Pairs :=
              Impulse_Amplify
                (Raw => Tide, Offset => 0.0, Ramp => 0.0,
                 Start => Template (Template'First).Date);
            F_Res : constant Data_Pairs :=
              IIR (Raw => Shape (Imp), lagA => Der, lagC => Lag_C,
                   iA => Init_Val, Start => IDate);
         begin
            for I in Template'Range loop
               Out_F (I - Template'First + 1) := F_Res (I).Value;
            end loop;
         end;
      end if;
      return Out_F;
   end Calc_Pre_Bessel;

   --  3-point binomial filter matching GEM.LTE.Primitives.Filter9Point exactly
   procedure Apply_F9 (V : in out Real_Vector; N_Pass : in Integer) is
      N : constant Integer := V'Length;
      Tmp : Real_Vector (1 .. N);
   begin
      if N_Pass <= 0 or else N < 3 then
         return;
      end if;
      for Pass in 1 .. N_Pass loop
         Tmp := V;
         for I in 2 .. N - 1 loop
            V (I) := 0.25 * Tmp (I - 1) + 0.5 * Tmp (I) + 0.25 * Tmp (I + 1);
         end loop;
      end loop;
   end Apply_F9;

   --  Exact 9-point filter using Primitives.Filter9Point on a Data_Pairs buffer
   procedure Apply_F9_Exact
     (Dates : in Real_Vector; V : in out Real_Vector; N_Pass : in Integer)
   is
      N : constant Integer := V'Length;
      DP : Data_Pairs (1 .. N);
   begin
      if N_Pass <= 0 then
         return;
      end if;
      for I in 1 .. N loop
         DP (I) := (Date => Dates (I), Value => V (I));
      end loop;
      for Pass in 1 .. N_Pass loop
         DP := Filter9Point (DP);
      end loop;
      for I in 1 .. N loop
         V (I) := DP (I).Value;
      end loop;
   end Apply_F9_Exact;

   function Pearson_CC
     (Y, D : in Real_Vector; Mask : in Real_Vector) return Long_Float
   is
      Cnt : Long_Float := 0.0;
      Sy, Sd, My, Md, Vyy, Vdd, Vyd : Long_Float := 0.0;
   begin
      for I in Y'Range loop
         if Mask (I) > 0.5 then
            Cnt := Cnt + 1.0;
            Sy := Sy + Y (I);
            Sd := Sd + D (I);
         end if;
      end loop;
      if Cnt < 3.0 then
         return 0.0;
      end if;
      My := Sy / Cnt;
      Md := Sd / Cnt;
      for I in Y'Range loop
         if Mask (I) > 0.5 then
            declare
               Dy : constant Long_Float := Y (I) - My;
               Dd : constant Long_Float := D (I) - Md;
            begin
               Vyy := Vyy + Dy * Dy;
               Vdd := Vdd + Dd * Dd;
               Vyd := Vyd + Dy * Dd;
            end;
         end if;
      end loop;
      if Vyy <= 1.0e-24 or else Vdd <= 1.0e-24 then
         return 0.0;
      end if;
      return Vyd / Sqrt (Vyy * Vdd);
   end Pearson_CC;

   --  Solve (A^T A + ridge*I) x = A^T b using RA.Solve with fallback
   function Solve_Normal
     (G : in Real_Matrix; Rhs : in Real_Vector; Ridge : in Long_Float := 1.0e-10)
      return Real_Vector
   is
      N : constant Integer := Rhs'Length;
      A : Real_Matrix := G;
      Eps : Long_Float := Ridge;
   begin
      for Attempt in 1 .. 6 loop
         for I in 1 .. N loop
            A (I, I) := G (I, I) + Eps;
         end loop;
         begin
            return Solve (A, Rhs);
         exception
            when others =>
               Eps := Eps * 100.0;
         end;
      end loop;
      return (1 .. N => 0.0);
   end Solve_Normal;

   procedure Run (D : in out Shared.Param_S) is
      T_Start : constant Ada.Calendar.Time := Ada.Calendar.Clock;
      NLP : constant Integer := D.NLP;
      NP  : constant Integer := 2 * NLP;
      NLT : constant Integer := D.NLT;

      --  Configuration
      N_Legs_Raw : constant Integer := GEM.Getenv ("LEGS", 0);
      N_Legs     : constant Positive :=
        (if N_Legs_Raw >= 2 then N_Legs_Raw else 1);
      Multi_Leg  : constant Boolean := N_Legs_Raw >= 2;
      Leg_P      : constant Shared.Leg_Own_Set := Leg_P_Mode;
      Has_Bessel : constant Boolean :=
        Multi_Leg and then (Leg_P = Shared.Own_Bessel or else Leg_P = Shared.Own_Params);

      Lam          : constant Long_Float := GEM.Getenv ("REFINE_LAM", 1.0);
      Do_Homotopy  : constant Boolean    := GEM.Getenv ("REFINE_HOMOTOPY", False);
      H_Lam_Raw    : constant Long_Float := GEM.Getenv ("REFINE_HOMOTOPY_LAM", -1.0);
      H_Lam        : constant Long_Float :=
        (if H_Lam_Raw >= 0.0 then H_Lam_Raw else Lam);
      Rounds       : constant Positive   := GEM.Getenv ("REFINE_ROUNDS", 4);
      Iters        : constant Positive   := GEM.Getenv ("REFINE_ITERS", 8);
      Init_Step    : constant Long_Float := GEM.Getenv ("REFINE_STEP", 0.3);
      Do_Per_Leg   : constant Boolean    := GEM.Getenv ("REFINE_PER_LEG", True);
      Init_LStep   : constant Long_Float := GEM.Getenv ("REFINE_LEG_STEP", 0.05);
      Leg_Lam      : constant Long_Float := GEM.Getenv ("REFINE_LEG_LAM", 1.0e-4);
      Do_Harm      : constant Boolean    := GEM.Getenv ("REFINE_HARM_SCREEN", True);
      Max_Swaps    : constant Natural    := GEM.Getenv ("REFINE_HARM_MAX_SWAPS", 3);
      Min_Gain     : constant Long_Float := GEM.Getenv ("REFINE_HARM_MIN_GAIN", 1.0e-4);
      Do_Apply     : constant Boolean    := GEM.Getenv ("REFINE_APPLY", True);

      F9_N         : constant Integer    := GEM.Getenv ("F9", 0);
      Exclude_On   : constant Boolean    := GEM.Getenv ("EXCLUDE", False);
      Train_Start  : constant Long_Float := GEM.Getenv ("TRAIN_START", 0.0);
      Train_End    : constant Long_Float := GEM.Getenv ("TRAIN_END", 9999.0);
      Strict_IDate : constant Boolean    := GEM.Getenv ("STRICT_IDATE", False);
      IDate        : constant Long_Float := GEM.Getenv ("IDATE", 1880.0);
      SPY          : constant Long_Float := GEM.Getenv ("SAMPLING", 12.0);
      Impulse_Only : constant Boolean    := GEM.Getenv ("IMPULSE", True);
      Imp_Mod12    : constant Boolean    := GEM.Getenv ("IMPULSE_MOD12", False);
      Trend_On     : constant Boolean    := GEM.Getenv ("TREND", True);
      Accel_On     : constant Boolean    := GEM.Getenv ("ACCEL", True);
      Max_H        : constant Integer    := GEM.Getenv ("MAXH", 30);
      IR_Lag       : constant Integer    := IR_Lag_Samples;
      NM_Cfg       : constant Integer    := GEM.Getenv ("NM", 3);

      --  Count NH from response file's NH list (e.g. "1 1 1 1 1 1 1")
      function Count_NH return Integer is
         S : constant String := GEM.Getenv ("NH", "1 1 1 1 1");
         C : Integer := 0;
         In_Tok : Boolean := False;
      begin
         for Ch of S loop
            if Ch = ' ' or else Ch = ASCII.HT or else Ch = '"' then
               In_Tok := False;
            elsif not In_Tok then
               In_Tok := True;
               C := C + 1;
            end if;
         end loop;
         return Integer'Max (1, Integer'Min (C, Shared.Max_Harmonics));
      end Count_NH;

      NH_Cfg : constant Integer := Count_NH;

      --  Per-leg workspace
      type Vec_Acc is access Real_Vector;
      type Mat_Acc is access Real_Matrix;
      procedure Free_Mat is new Ada.Unchecked_Deallocation (Real_Matrix, Mat_Acc);

      type Leg_State is record
         N       : Integer := 0;
         N_Train : Integer := 0;
         P       : access Shared.Param_S;
         T       : Vec_Acc;   -- dates (1..N)
         Data    : Vec_Acc;   -- F9-filtered climate data (1..N)
         Mask    : Vec_Acc;   -- 1.0 on training rows, 0.0 on held-out / zero
         F_Zero  : Vec_Acc;   -- homogeneous (lpap=0) response (1..N)
         F_Cur   : Vec_Acc;   -- current pre-Bessel F (1..N)
         J       : Mat_Acc;   -- (1..N, 1..NP) exact analytic pre-Bessel Jacobian
         CW_On   : Boolean := False;
         CW_Freq : Long_Float := 0.0;
         NM      : Integer := 3;
         NH      : Integer := 5;
         Acc_Ref : Long_Float := 1880.0;
         --  Per-leg continuous state (Block B)
         ImpA, ImpB, BesK, Offs, BG, BG2, IR : Long_Float := 0.0;
         LT      : Real_Vector (1 .. 16) := (others => 0.0);
         Harm    : Ns (1 .. Shared.Max_Harmonics) := (others => 0);
      end record;

      Legs : array (1 .. N_Legs) of Leg_State;

      --  Shared dLOD prior Gram matrix Phi_G (NP x NP) and reference Th_Ref
      Phi_G  : Real_Matrix (1 .. NP, 1 .. NP) := (others => (others => 0.0));
      Th_Ref : Real_Vector (1 .. NP) := (others => 0.0);
      Th_0   : Real_Vector (1 .. NP) := (others => 0.0);
      X_Scl  : Real_Vector (1 .. NP) := (others => 1.0e-3);
      Free_M : array (1 .. NP) of Boolean := (others => True);
      N_Free : Integer := 0;

      function DLOD_CC (Th : in Real_Vector) return Long_Float is
         P_Th  : constant Real_Vector := Phi_G * Th;
         P_Ref : constant Real_Vector := Phi_G * Th_Ref;
         Num   : constant Long_Float  := Th * P_Ref;
         D1    : constant Long_Float  := Th * P_Th;
         D2    : constant Long_Float  := Th_Ref * P_Ref;
      begin
         if D1 <= 1.0e-24 or else D2 <= 1.0e-24 then
            return 0.0;
         end if;
         return Num / Sqrt (D1 * D2);
      end DLOD_CC;

      --  Apply Bessel correction to pre-Bessel F for Leg L
      procedure Post_Bessel
        (L : in Leg_State; F_In : in Real_Vector; F_Out : out Real_Vector)
      is
         K_B  : constant Long_Float :=
           (if L.NM = 1 then L.BesK * (1.0 - L.BG) else L.BesK);
         E2S  : constant Long_Float := L.Offs;
         E2C  : constant Long_Float := (if L.NM = 1 then L.BG2 else L.BG);
         Arg  : Long_Float;
      begin
         for I in F_In'Range loop
            Arg := Two_Pi * K_B * F_In (I);
            F_Out (I) :=
              F_In (I)
              + L.ImpA * Sin (Arg) + L.ImpB * Cos (Arg)
              + E2S * L.ImpA * Sin (2.0 * Arg) + E2C * L.ImpB * Cos (2.0 * Arg);
         end loop;
      end Post_Bessel;

      --  Evaluate Leg L for a given pre-Bessel F_In, optionally capping |k| <= Kmax
      procedure Eval_Leg
        (L      : in Leg_State;
         F_In   : in Real_Vector;
         Model  : out Real_Vector;
         CC_Tr  : out Long_Float;
         Kmax   : in Long_Float := 1.0e9)
      is
         N : constant Integer := L.N;
         Fp : Real_Vector (1 .. N);
         --  Active wavenumbers: LT(1..NM) + LT(NM)*Harm(1..NH)
         Ks : Real_Vector (1 .. 64);
         NK : Integer := 0;
         K_Base : constant Long_Float := L.LT (L.NM);
      begin
         Post_Bessel (L, F_In, Fp);
         for J in 1 .. L.NM loop
            if abs L.LT (J) > 0.0 and then abs L.LT (J) <= Kmax then
               NK := NK + 1;
               Ks (NK) := L.LT (J);
            end if;
         end loop;
         for J in 1 .. L.NH loop
            if L.Harm (J) /= 0 then
               declare
                  Kh : constant Long_Float := K_Base * Long_Float (L.Harm (J));
               begin
                  if abs Kh <= Kmax then
                     NK := NK + 1;
                     Ks (NK) := Kh;
                  end if;
               end;
            end if;
         end loop;

         --  Build design matrix X (N x NC):
         --    1, Fp, sin/cos(2pi*k*Fp), [t, (t-Acc_Ref)^2, ann1/2, sem1/2], [cw1/2]
         declare
            NC : Integer := 2 + 2 * NK;
            Col : Integer;
         begin
            if Trend_On then
               NC := NC + (if Accel_On then 6 else 5);
            end if;
            if L.CW_On then
               NC := NC + 2;
            end if;

            declare
               X   : Real_Matrix (1 .. N, 1 .. NC);
               XtX : Real_Matrix (1 .. NC, 1 .. NC) := (others => (others => 0.0));
               XtY : Real_Vector (1 .. NC) := (others => 0.0);
               Y_Tgt : Real_Vector (1 .. N);
               Beta : Real_Vector (1 .. NC);
               Lag : constant Integer := Integer'Max (1, IR_Lag);
            begin
               --  Target y_ir = Data - IR * Data(t - lag) on training rows?
               --  Wait: in lt.exe, Model = LTE(Fp) + IR * Data(t - lag), so
               --  OLS regresses (Data - IR * Data_lag) on X!
               for I in 1 .. N loop
                  if L.IR /= 0.0 and then I > Lag then
                     Y_Tgt (I) := L.Data (I) + L.IR * L.Data (I - Lag);
                  else
                     Y_Tgt (I) := L.Data (I);
                  end if;
                  X (I, 1) := 1.0;
                  X (I, 2) := Fp (I);
                  Col := 2;
                  for J in 1 .. NK loop
                     declare
                        A : constant Long_Float := Two_Pi * Ks (J) * Fp (I);
                     begin
                        Col := Col + 1; X (I, Col) := Sin (A);
                        Col := Col + 1; X (I, Col) := Cos (A);
                     end;
                  end loop;
                  if Trend_On then
                     declare
                        Ti : constant Long_Float := L.T (I);
                        Dt : constant Long_Float := Ti - L.Acc_Ref;
                     begin
                        Col := Col + 1; X (I, Col) := Ti - L.Acc_Ref;
                        if Accel_On then
                           Col := Col + 1; X (I, Col) := Dt * Dt;
                        end if;
                        Col := Col + 1; X (I, Col) := Sin (Two_Pi * Ti);
                        Col := Col + 1; X (I, Col) := Cos (Two_Pi * Ti);
                        Col := Col + 1; X (I, Col) := Sin (4.0 * Pi * Ti);
                        Col := Col + 1; X (I, Col) := Cos (4.0 * Pi * Ti);
                     end;
                  end if;
                  if L.CW_On then
                     declare
                        Acw : constant Long_Float := Two_Pi * L.CW_Freq * L.T (I);
                     begin
                        Col := Col + 1; X (I, Col) := Sin (Acw);
                        Col := Col + 1; X (I, Col) := Cos (Acw);
                     end;
                  end if;
               end loop;

               --  Accumulate normal equations on training mask
               for I in 1 .. N loop
                  if L.Mask (I) > 0.5 then
                     declare
                        Yi : constant Long_Float := Y_Tgt (I);
                     begin
                        for C1 in 1 .. NC loop
                           declare
                              X1 : constant Long_Float := X (I, C1);
                           begin
                              XtY (C1) := XtY (C1) + X1 * Yi;
                              for C2 in C1 .. NC loop
                                 XtX (C1, C2) := XtX (C1, C2) + X1 * X (I, C2);
                              end loop;
                           end;
                        end loop;
                     end;
                  end if;
               end loop;
               for C1 in 1 .. NC loop
                  for C2 in 1 .. C1 - 1 loop
                     XtX (C1, C2) := XtX (C2, C1);
                  end loop;
               end loop;

               Beta := Solve_Normal (XtX, XtY, Ridge => 1.0e-10);

               declare
                  M_Raw : Real_Vector (1 .. N);
               begin
                  for I in 1 .. N loop
                     declare
                        S : Long_Float := 0.0;
                     begin
                        for C1 in 1 .. NC loop
                           S := S + X (I, C1) * Beta (C1);
                        end loop;
                        M_Raw (I) := S;
                     end;
                  end loop;
                  for I in 1 .. N loop
                     if L.IR /= 0.0 and then I > Lag then
                        Model (I) := M_Raw (I) - L.IR * M_Raw (I - Lag);
                     else
                        Model (I) := M_Raw (I);
                     end if;
                  end loop;
               end;

               if F9_N > 0 then
                  Apply_F9 (Model, F9_N);
               end if;
               CC_Tr := Pearson_CC (Model, L.Data.all, L.Mask.all);
            end;
         end;
      end Eval_Leg;

      --  Compute total normalized residual vector + dLOD prior cost
      function Eval_Mean_CC
        (Th : in Real_Vector; Kmax : in Long_Float := 1.0e9) return Long_Float
      is
         Sum_CC : Long_Float := 0.0;
      begin
         for L_Idx in 1 .. N_Legs loop
            declare
               L : Leg_State renames Legs (L_Idx);
               F : constant Real_Vector := L.F_Zero.all + L.J.all * Th;
               M : Real_Vector (1 .. L.N);
               C : Long_Float;
            begin
               Eval_Leg (L, F, M, C, Kmax => Kmax);
               Sum_CC := Sum_CC + C;
            end;
         end loop;
         return Sum_CC / Long_Float (N_Legs);
      end Eval_Mean_CC;

      --  One damped trust-region Gauss-Newton step on shared lpap (Block A)
      procedure Block_A_GN
        (Th     : in out Real_Vector;
         Cur_L  : in Long_Float;
         Step_S : in Long_Float;
         N_Iter : in Positive;
         Kmax   : in Long_Float := 1.0e9)
      is
         --  Total training samples across all legs
         Tot_R : Integer := 0;
         Free_Idx : array (1 .. NP) of Integer;
         NF : Integer := 0;
      begin
         for L_Idx in 1 .. N_Legs loop
            Tot_R := Tot_R + Legs (L_Idx).N_Train;
         end loop;
         for J in 1 .. NP loop
            if Free_M (J) then
               NF := NF + 1;
               Free_Idx (NF) := J;
            end if;
         end loop;
         if NF = 0 or else Tot_R = 0 then
            return;
         end if;

         declare
            R0  : Real_Vector (1 .. Tot_R);
            R1  : Real_Vector (1 .. Tot_R);
            JR_Ptr : Mat_Acc := new Real_Matrix (1 .. Tot_R, 1 .. NF);
            JR  : Real_Matrix renames JR_Ptr.all;
            Mu  : Long_Float := 1.0e-2;

            procedure Fill_Resid
              (Th_In : in Real_Vector; Res : out Real_Vector; Obj : out Long_Float)
            is
               Pos : Integer := 0;
               S_Sq : Long_Float := 0.0;
            begin
               for L_Idx in 1 .. N_Legs loop
                  declare
                     L : Leg_State renames Legs (L_Idx);
                     F : constant Real_Vector := L.F_Zero.all + L.J.all * Th_In;
                     M : Real_Vector (1 .. L.N);
                     C : Long_Float;
                     Cnt, Sy, Sd, My, Md, Vy, Vd : Long_Float := 0.0;
                     W : constant Long_Float :=
                       1.0 / Sqrt (2.0 * Long_Float (L.N_Train * N_Legs));
                  begin
                     Eval_Leg (L, F, M, C, Kmax => Kmax);
                     for I in 1 .. L.N loop
                        if L.Mask (I) > 0.5 then
                           Cnt := Cnt + 1.0;
                           Sy := Sy + M (I);
                           Sd := Sd + L.Data (I);
                        end if;
                     end loop;
                     My := Sy / Cnt; Md := Sd / Cnt;
                     for I in 1 .. L.N loop
                        if L.Mask (I) > 0.5 then
                           Vy := Vy + (M (I) - My) ** 2;
                           Vd := Vd + (L.Data (I) - Md) ** 2;
                        end if;
                     end loop;
                     declare
                        Sy_S : constant Long_Float := Long_Float'Max (Sqrt (Vy / Cnt), 1.0e-12);
                        Sd_S : constant Long_Float := Long_Float'Max (Sqrt (Vd / Cnt), 1.0e-12);
                     begin
                        for I in 1 .. L.N loop
                           if L.Mask (I) > 0.5 then
                              Pos := Pos + 1;
                              Res (Pos) :=
                                ((M (I) - My) / Sy_S - (L.Data (I) - Md) / Sd_S) * W;
                              S_Sq := S_Sq + Res (Pos) * Res (Pos);
                           end if;
                        end loop;
                     end;
                  end;
               end loop;
               declare
                  Dth  : constant Real_Vector := Th_In - Th_Ref;
                  Prior : constant Long_Float :=
                    Cur_L * (Dth * (Phi_G * Dth));
               begin
                  Obj := 0.5 * (S_Sq + Prior);
               end;
            end Fill_Resid;

            Obj0, Obj1 : Long_Float;
            Th_Try : Real_Vector (1 .. NP);
            N_Workers : constant Positive :=
              Integer'Max
                (1, Integer'Min
                   (GEM.Getenv ("NUMBER_OF_PROCESSORS",
                                System.Task_Info.Number_Of_Processors),
                    NF));
         begin
            Fill_Resid (Th, R0, Obj0);
            for It in 1 .. N_Iter loop
               --  Build Jacobian JR (Tot_R x NF) in parallel across N_Workers Ada tasks
               if N_Workers = 1 then
                  for C in 1 .. NF loop
                     declare
                        J_Idx : constant Integer := Free_Idx (C);
                        H     : constant Long_Float := 1.0e-5 * X_Scl (J_Idx);
                        Dummy : Long_Float;
                     begin
                        Th_Try := Th;
                        Th_Try (J_Idx) := Th (J_Idx) + H;
                        Fill_Resid (Th_Try, R1, Dummy);
                        for Row in 1 .. Tot_R loop
                           JR (Row, C) := (R1 (Row) - R0 (Row)) / H;
                        end loop;
                     end;
                  end loop;
               else
                  declare
                     task type JR_Worker is
                        pragma Storage_Size (8_388_608);
                        entry Start (Lo, Hi : in Integer);
                     end JR_Worker;

                     task body JR_Worker is
                        C_Lo, C_Hi : Integer;
                        Th_Loc : Real_Vector (1 .. NP);
                        R_Loc  : Real_Vector (1 .. Tot_R);
                        Dummy  : Long_Float;
                     begin
                        accept Start (Lo, Hi : in Integer) do
                           C_Lo := Lo;
                           C_Hi := Hi;
                        end Start;
                        for C in C_Lo .. C_Hi loop
                           declare
                              J_Idx : constant Integer := Free_Idx (C);
                              H     : constant Long_Float := 1.0e-5 * X_Scl (J_Idx);
                           begin
                              Th_Loc := Th;
                              Th_Loc (J_Idx) := Th (J_Idx) + H;
                              Fill_Resid (Th_Loc, R_Loc, Dummy);
                              for Row in 1 .. Tot_R loop
                                 JR (Row, C) := (R_Loc (Row) - R0 (Row)) / H;
                              end loop;
                           end;
                        end loop;
                     end JR_Worker;

                     Workers : array (1 .. N_Workers) of JR_Worker;
                     Chunk   : constant Integer := (NF + N_Workers - 1) / N_Workers;
                     Lo_C, Hi_C : Integer;
                  begin
                     for W in 1 .. N_Workers loop
                        Lo_C := (W - 1) * Chunk + 1;
                        Hi_C := Integer'Min (W * Chunk, NF);
                        if Lo_C <= NF then
                           Workers (W).Start (Lo_C, Hi_C);
                        else
                           Workers (W).Start (1, 0);
                        end if;
                     end loop;
                  end;
               end if;

               --  Normal equations: (JR^T JR + Prior_H + Mu * diag) d = -(JR^T R0 + Prior_g)
               declare
                  H_Mat : Real_Matrix (1 .. NF, 1 .. NF) := (others => (others => 0.0));
                  Grad  : Real_Vector (1 .. NF) := (others => 0.0);
                  Dth_R : constant Real_Vector := Th - Th_Ref;
                  Pg_Full : constant Real_Vector := Phi_G * Dth_R;
                  P_Scale : constant Long_Float := Cur_L;
                  Step_V  : Real_Vector (1 .. NF);
                  Accepted : Boolean := False;
               begin
                  for C1 in 1 .. NF loop
                     declare
                        S_G : Long_Float := 0.0;
                     begin
                        for Row in 1 .. Tot_R loop
                           S_G := S_G + JR (Row, C1) * R0 (Row);
                        end loop;
                        Grad (C1) := S_G + P_Scale * Pg_Full (Free_Idx (C1));
                        for C2 in C1 .. NF loop
                           declare
                              S_H : Long_Float := 0.0;
                           begin
                              for Row in 1 .. Tot_R loop
                                 S_H := S_H + JR (Row, C1) * JR (Row, C2);
                              end loop;
                              H_Mat (C1, C2) :=
                                S_H + P_Scale * Phi_G (Free_Idx (C1), Free_Idx (C2));
                           end;
                        end loop;
                     end;
                  end loop;
                  for C1 in 1 .. NF loop
                     for C2 in 1 .. C1 - 1 loop
                        H_Mat (C1, C2) := H_Mat (C2, C1);
                     end loop;
                  end loop;

                  --  LM damping trial loop
                  for Try in 1 .. 5 loop
                     declare
                        H_Damp : Real_Matrix := H_Mat;
                     begin
                        for C in 1 .. NF loop
                           H_Damp (C, C) :=
                             H_Mat (C, C) + Mu * Long_Float'Max (H_Mat (C, C), 1.0e-6);
                        end loop;
                        Step_V := Solve_Normal (H_Damp, -Grad, Ridge => 1.0e-12);
                     end;

                     --  Box-clamp to trust region Step_S * X_Scl
                     Th_Try := Th;
                     for C in 1 .. NF loop
                        declare
                           J_Idx : constant Integer := Free_Idx (C);
                           Lim   : constant Long_Float := Step_S * X_Scl (J_Idx);
                           Dv    : Long_Float := Step_V (C);
                        begin
                           if Dv > Lim then Dv := Lim;
                           elsif Dv < -Lim then Dv := -Lim;
                           end if;
                           Th_Try (J_Idx) := Th (J_Idx) + Dv;
                        end;
                     end loop;

                     Fill_Resid (Th_Try, R1, Obj1);
                     if Obj1 < Obj0 - 1.0e-7 then
                        Th := Th_Try;
                        R0 := R1;
                        Obj0 := Obj1;
                        Mu := Long_Float'Max (Mu * 0.4, 1.0e-5);
                        Accepted := True;
                        exit;
                     else
                        Mu := Long_Float'Min (Mu * 4.0, 1.0e4);
                     end if;
                  end loop;

                  exit when not Accepted;
               end;
            end loop;
            Free_Mat (JR_Ptr);
         end;
      end Block_A_GN;

      --  Block B: per-leg continuous damped GN on (Bessel, ltep, IR)
      procedure Block_B_Leg
        (L_Idx  : in Positive;
         Th_Cur : in Real_Vector;
         L_Step : in Long_Float)
      is
         L : Leg_State renames Legs (L_Idx);
         F_In : constant Real_Vector := L.F_Zero.all + L.J.all * Th_Cur;
         --  Pack active continuous parameters into X (up to 16)
         type P_Kind is (K_ImpA, K_ImpB, K_BesK, K_Offs, K_BG, K_LT, K_IR);
         type P_Key is record
            Kind : P_Kind;
            Idx  : Integer := 0;
         end record;
         Keys : array (1 .. 16) of P_Key;
         NK   : Integer := 0;

         procedure Add_Key (K : P_Kind; I : Integer := 0) is
         begin
            NK := NK + 1;
            Keys (NK) := (Kind => K, Idx => I);
         end Add_Key;

         function Get_Vec return Real_Vector is
            V : Real_Vector (1 .. NK);
         begin
            for I in 1 .. NK loop
               case Keys (I).Kind is
                  when K_ImpA => V (I) := L.ImpA;
                  when K_ImpB => V (I) := L.ImpB;
                  when K_BesK => V (I) := L.BesK;
                  when K_Offs => V (I) := L.Offs;
                  when K_BG   => V (I) := L.BG;
                  when K_LT   => V (I) := L.LT (Keys (I).Idx);
                  when K_IR   => V (I) := L.IR;
               end case;
            end loop;
            return V;
         end Get_Vec;

         procedure Set_Vec (V : in Real_Vector) is
         begin
            for I in 1 .. NK loop
               case Keys (I).Kind is
                  when K_ImpA => L.ImpA := V (I);
                  when K_ImpB => L.ImpB := V (I);
                  when K_BesK => L.BesK := V (I);
                  when K_Offs => L.Offs := V (I);
                  when K_BG   => L.BG   := V (I);
                  when K_LT   => L.LT (Keys (I).Idx) := V (I);
                  when K_IR   => L.IR   := V (I);
               end case;
            end loop;
         end Set_Vec;

         M_Buf : Real_Vector (1 .. L.N);
      begin
         if Has_Bessel or else N_Legs = 1 then
            Add_Key (K_ImpA);
            Add_Key (K_ImpB);
            Add_Key (K_BesK);
            Add_Key (K_Offs);
            Add_Key (K_BG);
         end if;
         for J in 1 .. L.NM loop
            if abs L.LT (J) > 0.0 then
               Add_Key (K_LT, J);
            end if;
         end loop;
         Add_Key (K_IR);

         declare
            X0  : Real_Vector := Get_Vec;
            X_S : Real_Vector (1 .. NK);
            R0  : Real_Vector (1 .. L.N_Train);
            R1  : Real_Vector (1 .. L.N_Train);
            JR  : Real_Matrix (1 .. L.N_Train, 1 .. NK);
            C0, C1 : Long_Float;
            Mu  : Long_Float := 1.0e-2;

            procedure Fill_R (X_In : in Real_Vector; Res : out Real_Vector; C_Out : out Long_Float) is
               Pos : Integer := 0;
               Cnt, Sy, Sd, My, Md, Vy, Vd : Long_Float := 0.0;
               W : constant Long_Float := 1.0 / Sqrt (2.0 * Long_Float (L.N_Train));
            begin
               Set_Vec (X_In);
               Eval_Leg (L, F_In, M_Buf, C_Out);
               for I in 1 .. L.N loop
                  if L.Mask (I) > 0.5 then
                     Cnt := Cnt + 1.0; Sy := Sy + M_Buf (I); Sd := Sd + L.Data (I);
                  end if;
               end loop;
               My := Sy / Cnt; Md := Sd / Cnt;
               for I in 1 .. L.N loop
                  if L.Mask (I) > 0.5 then
                     Vy := Vy + (M_Buf (I) - My) ** 2;
                     Vd := Vd + (L.Data (I) - Md) ** 2;
                  end if;
               end loop;
               declare
                  Sy_S : constant Long_Float := Long_Float'Max (Sqrt (Vy / Cnt), 1.0e-12);
                  Sd_S : constant Long_Float := Long_Float'Max (Sqrt (Vd / Cnt), 1.0e-12);
               begin
                  for I in 1 .. L.N loop
                     if L.Mask (I) > 0.5 then
                        Pos := Pos + 1;
                        Res (Pos) := ((M_Buf (I) - My) / Sy_S - (L.Data (I) - Md) / Sd_S) * W;
                     end if;
                  end loop;
               end;
            end Fill_R;

            X_Cur : Real_Vector := X0;
            X_Try : Real_Vector (1 .. NK);
         begin
            for I in 1 .. NK loop
               X_S (I) := Long_Float'Max (abs X0 (I), 0.05);
            end loop;
            Fill_R (X_Cur, R0, C0);

            for It in 1 .. 6 loop
               for C in 1 .. NK loop
                  declare
                     H : constant Long_Float := 1.0e-5 * X_S (C);
                     D_C : Long_Float;
                  begin
                     X_Try := X_Cur;
                     X_Try (C) := X_Cur (C) + H;
                     Fill_R (X_Try, R1, D_C);
                     for Row in 1 .. L.N_Train loop
                        JR (Row, C) := (R1 (Row) - R0 (Row)) / H;
                     end loop;
                  end;
               end loop;

               declare
                  H_Mat : Real_Matrix (1 .. NK, 1 .. NK) := (others => (others => 0.0));
                  Grad  : Real_Vector (1 .. NK) := (others => 0.0);
                  Step_V : Real_Vector (1 .. NK);
                  Accepted : Boolean := False;
               begin
                  for C1_I in 1 .. NK loop
                     for Row in 1 .. L.N_Train loop
                        Grad (C1_I) := Grad (C1_I) + JR (Row, C1_I) * R0 (Row);
                     end loop;
                     Grad (C1_I) := Grad (C1_I) + Leg_Lam * (X_Cur (C1_I) - X0 (C1_I)) / (X_S (C1_I) ** 2);
                     for C2_I in C1_I .. NK loop
                        declare
                           S_H : Long_Float := 0.0;
                        begin
                           for Row in 1 .. L.N_Train loop
                              S_H := S_H + JR (Row, C1_I) * JR (Row, C2_I);
                           end loop;
                           H_Mat (C1_I, C2_I) := S_H;
                        end;
                     end loop;
                     H_Mat (C1_I, C1_I) := H_Mat (C1_I, C1_I) + Leg_Lam / (X_S (C1_I) ** 2);
                  end loop;
                  for C1_I in 1 .. NK loop
                     for C2_I in 1 .. C1_I - 1 loop
                        H_Mat (C1_I, C2_I) := H_Mat (C2_I, C1_I);
                     end loop;
                  end loop;

                  for Try in 1 .. 4 loop
                     declare
                        H_Damp : Real_Matrix := H_Mat;
                     begin
                        for C in 1 .. NK loop
                           H_Damp (C, C) := H_Mat (C, C) + Mu * Long_Float'Max (H_Mat (C, C), 1.0e-6);
                        end loop;
                        Step_V := Solve_Normal (H_Damp, -Grad, Ridge => 1.0e-12);
                     end;
                     X_Try := X_Cur;
                     for C in 1 .. NK loop
                        declare
                           Lim : constant Long_Float := L_Step * X_S (C);
                           Dv  : Long_Float := Step_V (C);
                           Val : Long_Float;
                        begin
                           if Dv > Lim then Dv := Lim;
                           elsif Dv < -Lim then Dv := -Lim;
                           end if;
                           Val := X_Cur (C) + Dv;
                           --  Keep |ltep| >= Min_LT and preserve sign
                           if Keys (C).Kind = K_LT then
                              if X0 (C) > 0.0 and then Val < Min_LT then
                                 Val := Min_LT;
                              elsif X0 (C) < 0.0 and then Val > -Min_LT then
                                 Val := -Min_LT;
                              end if;
                           elsif Keys (C).Kind = K_IR then
                              if Val < -0.95 then Val := -0.95;
                              elsif Val > 0.95 then Val := 0.95;
                              end if;
                           end if;
                           X_Try (C) := Val;
                        end;
                     end loop;

                     Fill_R (X_Try, R1, C1);
                     if C1 > C0 + 1.0e-6 then
                        X_Cur := X_Try;
                        R0 := R1;
                        C0 := C1;
                        Mu := Long_Float'Max (Mu * 0.4, 1.0e-5);
                        Accepted := True;
                        exit;
                     else
                        Mu := Long_Float'Min (Mu * 4.0, 1.0e4);
                     end if;
                  end loop;

                  exit when not Accepted;
               end;
            end loop;

            Set_Vec (X_Cur);
         end;
      end Block_B_Leg;

      --  Block C: greedy 1-swap harmonic screening on Leg L_Idx
      procedure Block_C_Harm
        (L_Idx : in Positive; Th_Cur : in Real_Vector)
      is
         L : Leg_State renames Legs (L_Idx);
         F_In : constant Real_Vector := L.F_Zero.all + L.J.all * Th_Cur;
         M_Buf : Real_Vector (1 .. L.N);
         Base_CC, Cand_CC : Long_Float;
         Best_CC : Long_Float;
         Best_Slot, Best_H : Integer;
      begin
         Eval_Leg (L, F_In, M_Buf, Base_CC);
         for Swap in 1 .. Max_Swaps loop
            Best_CC := Base_CC;
            Best_Slot := 0;
            Best_H := 0;
            for Slot in 1 .. L.NH loop
               declare
                  Old_H : constant Integer := L.Harm (Slot);
               begin
                  for H in 2 .. Max_H loop
                     declare
                        Used : Boolean := False;
                     begin
                        for J in 1 .. L.NH loop
                           if L.Harm (J) = H then
                              Used := True;
                              exit;
                           end if;
                        end loop;
                        if not Used then
                           L.Harm (Slot) := H;
                           Eval_Leg (L, F_In, M_Buf, Cand_CC);
                           if Cand_CC > Best_CC + Min_Gain then
                              Best_CC := Cand_CC;
                              Best_Slot := Slot;
                              Best_H := H;
                           end if;
                        end if;
                     end;
                  end loop;
                  L.Harm (Slot) := Old_H;
               end;
            end loop;
            exit when Best_Slot = 0;
            L.Harm (Best_Slot) := Best_H;
            Base_CC := Best_CC;
         end loop;
      end Block_C_Harm;

      CC0, CC1 : array (1 .. N_Legs) of Long_Float := (others => 0.0);
      Mean0, Mean1, DL0, DL1 : Long_Float;
      Th_Best : Real_Vector (1 .. NP);
      Step_Cur : Long_Float := Init_Step;
      LStep_Cur : Long_Float := Init_LStep;
   begin
      Ada.Text_IO.Put_Line
        ("REFINE: " & N_Legs'Img & " leg(s), lam=" & Fmt (Lam, 3)
         & ", homotopy=" & Boolean'Image (Do_Homotopy)
         & ", rounds=" & Rounds'Img & ", iters=" & Iters'Img);

      --  1. Build dLOD prior Gram matrix Phi_G from dlod3.dat
      declare
         DL_Dat : constant String := GEM.Getenv ("DLOD_DAT", "../dlod3.dat");
         DL_Raw : constant Data_Pairs := Make_Data (DL_Dat);
         Stride : constant Integer := Integer'Max (1, DL_Raw'Length / 2000);
         N_Sub  : constant Integer := (DL_Raw'Length + Stride - 1) / Stride;
         Sub_T  : Data_Pairs (1 .. N_Sub);
         B_Ptr  : Mat_Acc := new Real_Matrix (1 .. N_Sub, 1 .. NP);
         B_Mat  : Real_Matrix renames B_Ptr.all;
         Unit_AP : Long_Periods_Amp_Phase (1 .. NLP) :=
           (others => (Amplitude => 0.0, Phase => 0.0));
         Y_Len  : constant Long_Float := Year_Length (D.B.Year);
      begin
         for I in 1 .. N_Sub loop
            Sub_T (I) := DL_Raw (DL_Raw'First + (I - 1) * Stride);
         end loop;
         for J in 1 .. NLP loop
            Unit_AP (J) := (Amplitude => 1.0, Phase => 0.0);
            declare
               T1 : constant Data_Pairs :=
                 Tide_Sum (Template => Sub_T, Constituents => Unit_AP,
                           Periods => D.A.LP, Ref_Time => 0.0, Scaling => 1.0,
                           Cos_Phase => False, Year_Len => Y_Len);
            begin
               for I in 1 .. N_Sub loop
                  B_Mat (I, 2 * J - 1) := T1 (I).Value;
               end loop;
            end;
            Unit_AP (J) := (Amplitude => 1.0, Phase => Pi / 2.0);
            declare
               T2 : constant Data_Pairs :=
                 Tide_Sum (Template => Sub_T, Constituents => Unit_AP,
                           Periods => D.A.LP, Ref_Time => 0.0, Scaling => 1.0,
                           Cos_Phase => False, Year_Len => Y_Len);
            begin
               for I in 1 .. N_Sub loop
                  B_Mat (I, 2 * J) := T2 (I).Value;
               end loop;
            end;
            Unit_AP (J) := (Amplitude => 0.0, Phase => 0.0);
         end loop;

         --  Mean-center columns and form Phi_G = (1/N_Sub) * B^T B
         for C in 1 .. NP loop
            declare
               S : Long_Float := 0.0;
            begin
               for I in 1 .. N_Sub loop S := S + B_Mat (I, C); end loop;
               S := S / Long_Float (N_Sub);
               for I in 1 .. N_Sub loop B_Mat (I, C) := B_Mat (I, C) - S; end loop;
            end;
         end loop;
         for C1 in 1 .. NP loop
            for C2 in C1 .. NP loop
               declare
                  S : Long_Float := 0.0;
               begin
                  for I in 1 .. N_Sub loop
                     S := S + B_Mat (I, C1) * B_Mat (I, C2);
                  end loop;
                  Phi_G (C1, C2) := S / Long_Float (N_Sub);
                  Phi_G (C2, C1) := Phi_G (C1, C2);
               end;
            end loop;
         end loop;
         Free_Mat (B_Ptr);
      end;

      LPAP_To_Theta (GEM.LTE.LPRef (1 .. NLP), Th_Ref);
      --  Normalize Phi_G so Th_Ref^T * Phi_G * Th_Ref = 1.0 (matches 1 - dLOD_CC scale)
      declare
         Ref_Var : constant Long_Float := Th_Ref * (Phi_G * Th_Ref);
      begin
         if Ref_Var > 1.0e-24 then
            Phi_G := (1.0 / Ref_Var) * Phi_G;
         end if;
      end;
      LPAP_To_Theta (D.B.LPAP, Th_0);
      Th_Best := Th_0;

      for J in 1 .. NLP loop
         declare
            Amp_R : constant Long_Float :=
              Sqrt (Th_Ref (2 * J - 1) ** 2 + Th_Ref (2 * J) ** 2);
            Amp_0 : constant Long_Float :=
              Sqrt (Th_0 (2 * J - 1) ** 2 + Th_0 (2 * J) ** 2);
            Sc : constant Long_Float := Long_Float'Max (Long_Float'Max (Amp_R, Amp_0), 1.0e-3);
         begin
            X_Scl (2 * J - 1) := Sc;
            X_Scl (2 * J)     := Sc;
            if Is_Fixed (D.A.LP (J)) then
               Free_M (2 * J - 1) := False;
               Free_M (2 * J)     := False;
            else
               N_Free := N_Free + 2;
            end if;
         end;
      end loop;

      --  2. Load each leg and compute exact analytic pre-Bessel Jacobian J
      for K in 1 .. N_Legs loop
         if Multi_Leg then
            GEM.Set_Leg (K);
         else
            GEM.Set_Leg (0);
         end if;
         declare
            Idx_File : constant String := GEM.Getenv ("CLIMATE_INDEX", "nino4.dat");
            Raw_Data : constant Data_Pairs := Make_Data (Idx_File);
            N : constant Integer := Raw_Data'Length;
            L : Leg_State renames Legs (K);
            Found : Boolean := False;
            Zero_AP : constant Long_Periods_Amp_Phase (1 .. NLP) :=
              (others => (Amplitude => 0.0, Phase => 0.0));
            Unit_AP : Long_Periods_Amp_Phase (1 .. NLP) := Zero_AP;
            Max_Diff : Long_Float := 0.0;
         begin
            L.N := N;
            L.P := new Shared.Param_S'(D);
            if Multi_Leg then
               Shared.Read_Leg (K, L.P.all, Found, Leg_P);
            end if;
            L.T      := new Real_Vector (1 .. N);
            L.Data   := new Real_Vector (1 .. N);
            L.Mask   := new Real_Vector (1 .. N);
            L.F_Zero := new Real_Vector (1 .. N);
            L.F_Cur  := new Real_Vector (1 .. N);
            L.J      := new Real_Matrix (1 .. N, 1 .. NP);

            for I in 1 .. N loop
               L.T (I)    := Raw_Data (Raw_Data'First + I - 1).Date;
               L.Data (I) := Raw_Data (Raw_Data'First + I - 1).Value;
            end loop;
            if F9_N > 0 then
               Apply_F9_Exact (L.T.all, L.Data.all, F9_N);
            end if;

            L.N_Train := 0;
            for I in 1 .. N loop
               declare
                  Ti : constant Long_Float := L.T (I);
                  In_Win : constant Boolean := (Ti >= Train_Start and then Ti <= Train_End);
                  Use_Row : constant Boolean :=
                    (L.Data (I) /= 0.0) and then (if Exclude_On then not In_Win else In_Win);
               begin
                  if Use_Row then
                     L.Mask (I) := 1.0;
                     L.N_Train := L.N_Train + 1;
                  else
                     L.Mask (I) := 0.0;
                  end if;
               end;
            end loop;
            L.Acc_Ref := L.T (1);
            L.NM      := Integer'Min (NM_Cfg, NLT);
            L.NH      := NH_Cfg;
            L.CW_On   := (L.P.CW1 /= 0.0);
            L.CW_Freq :=
              Year_Length (L.P.B.Year) / (GEM.LTE.Draconic / 2.0)
              - Long_Float'Floor (Year_Length (L.P.B.Year) / (GEM.LTE.Draconic / 2.0));
            L.ImpA := L.P.B.ImpA; L.ImpB := L.P.B.ImpB;
            L.BesK := L.P.B.besK; L.Offs := L.P.B.Offset;
            L.BG   := L.P.B.bg;   L.BG2  := L.P.B.bg2;
            L.IR   := L.P.B.IR;
            for J in 1 .. L.NM loop
               L.LT (J) := L.P.B.LT (J);
            end loop;
            L.Harm := L.P.C;

            --  Exact baseline forcing F_Cur (with real lagC & init) + 104-column
            --  analytic J (with lagC=0, init=0); F_Zero := F_Cur - J * Th_0 so
            --  F_Zero + J * Th = F_Cur + J * (Th - Th_0) is exact locally!
            L.F_Cur.all :=
              Calc_Pre_Bessel (Raw_Data, D.A.LP, D.B.LPAP, L.P.B, L.P.B.init,
                               Strict_IDate, IDate, SPY, Impulse_Only, Imp_Mod12,
                               Use_LagC => True);
            for J in 1 .. NLP loop
               Unit_AP (J) := (Amplitude => 1.0, Phase => 0.0);
               declare
                  C1 : constant Real_Vector :=
                    Calc_Pre_Bessel (Raw_Data, D.A.LP, Unit_AP, L.P.B, 0.0,
                                     Strict_IDate, IDate, SPY, Impulse_Only, Imp_Mod12,
                                     Use_LagC => False);
               begin
                  for I in 1 .. N loop L.J (I, 2 * J - 1) := C1 (I); end loop;
               end;
               Unit_AP (J) := (Amplitude => 1.0, Phase => Pi / 2.0);
               declare
                  C2 : constant Real_Vector :=
                    Calc_Pre_Bessel (Raw_Data, D.A.LP, Unit_AP, L.P.B, 0.0,
                                     Strict_IDate, IDate, SPY, Impulse_Only, Imp_Mod12,
                                     Use_LagC => False);
               begin
                  for I in 1 .. N loop L.J (I, 2 * J) := C2 (I); end loop;
               end;
               Unit_AP (J) := (Amplitude => 0.0, Phase => 0.0);
            end loop;

            L.F_Zero.all := L.F_Cur.all - L.J.all * Th_0;

            --  Verify analytic J via a 1e-4 perturbation spot-check on col 1
            declare
               Th_P : Real_Vector := Th_0;
               AP_P : Long_Periods_Amp_Phase := D.B.LPAP;
               H_Sc : constant Long_Float := 1.0e-4;
               M_Buf : Real_Vector (1 .. N);
            begin
               Th_P (1) := Th_0 (1) + H_Sc;
               Theta_To_LPAP (Th_P, AP_P);
               declare
                  F_Pert : constant Real_Vector :=
                    Calc_Pre_Bessel (Raw_Data, D.A.LP, AP_P, L.P.B, L.P.B.init,
                                     Strict_IDate, IDate, SPY, Impulse_Only, Imp_Mod12,
                                     Use_LagC => True);
                  F_Lin  : constant Real_Vector := L.F_Zero.all + L.J.all * Th_P;
               begin
                  for I in 1 .. N loop
                     Max_Diff := Long_Float'Max (Max_Diff, abs (F_Pert (I) - F_Lin (I)) / H_Sc);
                  end loop;
               end;
               Eval_Leg (L, L.F_Cur.all, M_Buf, CC0 (K));
            end;
            if K = 1 then
               Ada.Text_IO.Put_Line
                 ("  analytic J verified (leg 1 spot-check max diff = "
                  & Fmt (Max_Diff, 9) & ")");
            end if;
         end;
      end loop;
      GEM.Set_Leg (0);

      Mean0 := 0.0;
      for K in 1 .. N_Legs loop Mean0 := Mean0 + CC0 (K); end loop;
      Mean0 := Mean0 / Long_Float (N_Legs);
      DL0 := DLOD_CC (Th_0);
      Ada.Text_IO.Put_Line
        ("  baseline mean=" & Fmt (Mean0) & " dLOD=" & Fmt (DL0)
         & "  [" & Fmt (Long_Float (Ada.Calendar.Clock - T_Start), 2) & " s]");

      --  3. Optional coarse-to-fine FWI homotopy
      if Do_Homotopy then
         declare
            K_Stages : constant array (1 .. 3) of Long_Float := (0.35, 0.90, 2.50);
            S_Stages : constant array (1 .. 3) of Long_Float := (3.0,  2.5,  2.0);
         begin
            for St in K_Stages'Range loop
               Block_A_GN (Th_Best, H_Lam, S_Stages (St), Iters, Kmax => K_Stages (St));
               Ada.Text_IO.Put_Line
                 ("  homotopy stage" & St'Img & " (|k|<=" & Fmt (K_Stages (St), 2)
                  & ", lam=" & Fmt (H_Lam, 2) & "): full-spec mean="
                  & Fmt (Eval_Mean_CC (Th_Best)) & " dLOD=" & Fmt (DLOD_CC (Th_Best)));
            end loop;
         end;
      end if;

      --  4. Alternating Block C (harm) + Block A (lpap) + Block B (per-leg)
      for Rnd in 1 .. Rounds loop
         if Do_Harm then
            for K in 1 .. N_Legs loop
               Block_C_Harm (K, Th_Best);
            end loop;
         end if;

         declare
            Pre_M : constant Long_Float := Eval_Mean_CC (Th_Best);
         begin
            Block_A_GN (Th_Best, Lam, Step_Cur, Iters);
            if Do_Per_Leg then
               for K in 1 .. N_Legs loop
                  Block_B_Leg (K, Th_Best, LStep_Cur);
               end loop;
            end if;
            Mean1 := Eval_Mean_CC (Th_Best);
            if Mean1 > Pre_M then
               Step_Cur  := Long_Float'Min (Step_Cur * 1.5, 2.0);
               LStep_Cur := Long_Float'Min (LStep_Cur * 1.5, 0.25);
            else
               Step_Cur  := Long_Float'Max (Step_Cur * 0.5, 0.02);
               LStep_Cur := Long_Float'Max (LStep_Cur * 0.5, 0.005);
            end if;
            Ada.Text_IO.Put_Line
              ("  round" & Rnd'Img & ": mean " & Fmt (Pre_M) & " -> " & Fmt (Mean1)
               & " dLOD=" & Fmt (DLOD_CC (Th_Best))
               & "  [" & Fmt (Long_Float (Ada.Calendar.Clock - T_Start), 2) & " s]");
         end;
      end loop;

      --  5. Final scores and apply
      Mean1 := 0.0;
      for K in 1 .. N_Legs loop
         declare
            L : Leg_State renames Legs (K);
            F : constant Real_Vector := L.F_Zero.all + L.J.all * Th_Best;
            M : Real_Vector (1 .. L.N);
         begin
            Eval_Leg (L, F, M, CC1 (K));
            Mean1 := Mean1 + CC1 (K);
         end;
      end loop;
      Mean1 := Mean1 / Long_Float (N_Legs);
      DL1 := DLOD_CC (Th_Best);

      Ada.Text_IO.Put ("REFINE per-leg:");
      for K in 1 .. N_Legs loop
         Ada.Text_IO.Put (" " & Fmt (CC0 (K)) & "->" & Fmt (CC1 (K)));
      end loop;
      Ada.Text_IO.New_Line;
      Ada.Text_IO.Put_Line
        ("REFINE summary: mean " & Fmt (Mean0) & " -> " & Fmt (Mean1)
         & "; dLOD " & Fmt (DL0) & " -> " & Fmt (DL1)
         & "  [total " & Fmt (Long_Float (Ada.Calendar.Clock - T_Start), 2) & " s]");

      if Do_Apply and then Mean1 > Mean0 then
         Theta_To_LPAP (Th_Best, D.B.LPAP);
         --  Update leg 1 scalars on primary D
         D.B.ImpA := Legs (1).ImpA; D.B.ImpB := Legs (1).ImpB;
         D.B.besK := Legs (1).BesK; D.B.Offset := Legs (1).Offs;
         D.B.bg   := Legs (1).BG;   D.B.IR := Legs (1).IR;
         for J in 1 .. Legs (1).NM loop
            D.B.LT (J) := Legs (1).LT (J);
         end loop;
         D.C := Legs (1).Harm;

         if Ada.Directories.Exists ("lt.exe.p") then
            Ada.Directories.Copy_File ("lt.exe.p", "lt.exe.p.pre_refine");
         end if;
         Shared.Save (D);

         if Multi_Leg then
            for K in 1 .. N_Legs loop
               declare
                  L : Leg_State renames Legs (K);
                  Img : constant String :=
                    Ada.Strings.Fixed.Trim (Integer'Image (K), Ada.Strings.Both);
                  Fn : constant String := "lt.exe.p." & Img;
               begin
                  if Ada.Directories.Exists (Fn) then
                     Ada.Directories.Copy_File (Fn, Fn & ".pre_refine");
                  end if;
                  Theta_To_LPAP (Th_Best, L.P.B.LPAP);
                  L.P.B.ImpA := L.ImpA; L.P.B.ImpB := L.ImpB;
                  L.P.B.besK := L.BesK; L.P.B.Offset := L.Offs;
                  L.P.B.bg   := L.BG;   L.P.B.IR := L.IR;
                  for J in 1 .. L.NM loop
                     L.P.B.LT (J) := L.LT (J);
                  end loop;
                  L.P.C := L.Harm;
                  Shared.Write_Leg (K, L.P.all, Leg_P);
               end;
            end loop;
         end if;
         Ada.Text_IO.Put_Line ("REFINE applied: updated lt.exe.p (backups *.pre_refine)");
      end if;
   end Run;

end GEM.LTE.Primitives.Refine;
