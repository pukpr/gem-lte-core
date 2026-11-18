with Ada.Numerics.Long_Real_Arrays;
with Ada.Numerics.Long_Elementary_Functions;
with Ada.Long_Float_Text_IO;
with Text_IO;

package body GEM.LTE.Primitives.Inversion is

   use Ada.Numerics.Long_Real_Arrays;
   use Ada.Numerics.Long_Elementary_Functions;

   type Method_Kind is (GN, Viterbi);
   Max_Legs : constant := 9;
   N_Legs : constant Positive :=
     Integer'Max (1, Integer'Min (Max_Legs, Integer'(GEM.Getenv ("LEGS", 0))));
   Method_Asked : constant Method_Kind :=
     (if GEM.Getenv ("INVERSION_METHOD", "GN") = "VITERBI" then Viterbi
      else GN);
   --  VITERBI needs the legs aligned sample by sample; not done (yet)
   Method : constant Method_Kind := (if N_Legs >= 2 then GN else Method_Asked);
   Max_Iter : constant Integer := GEM.Getenv ("INVERSION_ITERS", 4);
   Prior : constant Long_Float := GEM.Getenv ("INVERSION_PRIOR", 0.2);
   Step_SD : constant Long_Float := GEM.Getenv ("INVERSION_STEP_SD", 0.06);
   Min_Gain : constant Long_Float := 1.0E-4;
   --  INVERSION_STOP VALID: keep the base with the best validation score and
   --  stop INVERSION_PATIENCE iterations after it; that base is saved
   Stop_Valid : constant Boolean := GEM.Getenv ("INVERSION_STOP", "TRAIN") = "VALID";
   Patience : constant Integer := GEM.Getenv ("INVERSION_PATIENCE", 2);
   Missing : constant Long_Float := Long_Float'First;  -- score not available
   --  INVERSION_ANCHOR: text file of "period re im" rows (an LPAP solution,
   --  e.g. the dLOD-only one, z = amplitude * exp (i phase)); the lines are
   --  pulled toward it by INVERSION_ANCHOR_W * mean (((z - z_anchor) / Scale)^2),
   --  in the GN normal equations (relative to the base residual sum) and,
   --  converted to correlation units around the seed score, in the training
   --  score that picks and accepts the steps. Lines without a row are free.
   Anchor_File : constant String := GEM.Getenv ("INVERSION_ANCHOR", "");
   Anchor_W : constant Long_Float := GEM.Getenv ("INVERSION_ANCHOR_W", 0.0);
   Anchored : constant Boolean := Anchor_File /= "" and then Anchor_W > 0.0;

   --  space-separated list of numbers
   function Floats (S : String) return Real_Vector is
      V : Real_Vector (1 .. 32);
      N : Natural := 0;
      I : Integer := S'First;
      J : Integer;
   begin
      while I <= S'Last loop
         if S (I) = ' ' then
            I := I + 1;
         else
            J := I;
            while J < S'Last and then S (J + 1) /= ' ' loop
               J := J + 1;
            end loop;
            N := N + 1;
            V (N) := Long_Float'Value (S (I .. J));
            I := J + 1;
         end if;
      end loop;
      return V (1 .. N);
   end Floats;

   Ridges : constant Real_Vector :=
     Floats (GEM.Getenv ("INVERSION_RIDGE", "0.001 0.01 0.1"));
   Steps : constant Real_Vector :=
     Floats (GEM.Getenv ("INVERSION_STEPS", "0.25 0.5 1"));
   NC : constant Natural := Ridges'Length * Steps'Length;

   --  Active parameters: Re/Im of a line's z (always adjacent, Re first),
   --  or one of the two Bessel amplitudes
   type Param_Kind is (Re, Im, Imp_A, Imp_B);
   Max_P : constant := 2 * 256 + 2;
   Kinds : array (1 .. Max_P) of Param_Kind;
   Lines : array (1 .. Max_P) of Natural := (others => 0);
   Eps, Scale : Real_Vector (1 .. Max_P) := (others => 0.0);
   NP : Natural := 0;
   --  anchor value of each Re/Im parameter, and whether it has one
   Z_Anchor : Real_Vector (1 .. Max_P) := (others => 0.0);
   Has_Anchor : array (1 .. Max_P) of Boolean := (others => False);
   N_Anchor : Natural := 0;
   Seed_Score : Long_Float := Missing;

   type Param_Ptr is access Shared.Param_S;
   type Index_List is array (Positive range <>) of Integer;
   type Index_Ptr is access Index_List;
   type Vec_Ptr is access Real_Vector;
   type Mat_Ptr is access Real_Matrix;

   Base : Param_Ptr := null;

   --  Per leg (1 = the leader, or the only task without legs); each leg's
   --  task writes only its own entry. The leader reads them all in Solve,
   --  after the pair's Result entry has waited for every follower's answer.
   type Leg_Store is record
      Tr : Index_Ptr := null;          -- training sample indices
      NT : Natural := 0;
      M0, F0, Target : Vec_Ptr := null;  -- base model, manifold, step target
      J : Mat_Ptr := null;             -- Jacobian (training samples x NP)
      Weight : Long_Float := 1.0;      -- sqrt (LEG_WEIGHT) / base residual rms
      Was_Base : Boolean := True;      -- what this leg's last evaluation was
      Was_Final : Boolean := False;    -- ... and whether it was the final one
      Val, Test : Long_Float := Missing;  -- its validation and test scores
      LW : Long_Float := 1.0;          -- its LEG_WEIGHT
   end record;
   Legs : array (1 .. Max_Legs) of Leg_Store;

   type Phase_Kind is (Want_Base, Probing, Trying);
   Phase : Phase_Kind := Want_Base;
   Probe, Cand, Iter : Natural := 0;
   Base_Score : Long_Float := 0.0;
   Cands : array (1 .. NC) of Param_Ptr := (others => null);
   Cand_Score, Cand_Gain, Cand_Ridge, Cand_Step :
     Real_Vector (1 .. NC) := (others => 0.0);
   Picked : Natural := 0;

   function Is_Base (Leg : Natural := 0) return Boolean is
     (Legs (Integer'Max (1, Leg)).Was_Base);

   --  The best-validation base, re-evaluated at the end, is accepted (saved)
   --  even though its training score is below the last base's
   Final : Boolean := False;
   function Force_Accept (Leg : Natural := 0) return Boolean is
     (Legs (Integer'Max (1, Leg)).Was_Final);

   Best_Val : Long_Float := Missing;
   Best_Val_Iter, Since_Best : Natural := 0;
   Best_Val_D : Param_Ptr := null;

   --  LEG_WEIGHT-weighted mean over the legs that have the score
   function Composite (Test : Boolean) return Long_Float is
      Sum, W : Long_Float := 0.0;
   begin
      for L in 1 .. N_Legs loop
         declare
            X : constant Long_Float := (if Test then Legs (L).Test else Legs (L).Val);
         begin
            if X > -2.0 then
               Sum := Sum + Legs (L).LW * X;
               W := W + Legs (L).LW;
            end if;
         end;
      end loop;
      return (if W > 0.0 then Sum / W else Missing);
   end Composite;

   procedure Put_F (X : Long_Float; Aft : Natural := 4) is
   begin
      Ada.Long_Float_Text_IO.Put (X, Fore => 2, Aft => Aft, Exp => 0);
   end Put_F;

   procedure Put_VT is
      V : constant Long_Float := Composite (False);
      T : constant Long_Float := Composite (True);
   begin
      Text_IO.Put (" valid");
      if V > -2.0 then Put_F (V); else Text_IO.Put ("    n/a"); end if;
      Text_IO.Put (" test");
      if T > -2.0 then Put_F (T); else Text_IO.Put ("    n/a"); end if;
   end Put_VT;

   --  z := z + (DX + i DY) for line L, keeping the sign convention of its
   --  amplitude and its phase within pi of the old one
   procedure Add_Z (D : in out Shared.Param_S; L : Positive; DX, DY : Long_Float) is
      A : constant Long_Float := D.B.LPAP (L).Amplitude;
      P : constant Long_Float := D.B.LPAP (L).Phase;
      X : constant Long_Float := A * Cos (P) + DX;
      Y : constant Long_Float := A * Sin (P) + DY;
      Sg : constant Long_Float := (if A >= 0.0 then 1.0 else -1.0);
      Two_Pi : constant Long_Float := 2.0 * Ada.Numerics.Pi;
      Q : Long_Float;
   begin
      D.B.LPAP (L).Amplitude := Sg * Sqrt (X * X + Y * Y);
      if X /= 0.0 or else Y /= 0.0 then
         Q := Arctan (Sg * Y, Sg * X) - P;
         Q := Q - Two_Pi * Long_Float'Floor (Q / Two_Pi + 0.5);
         D.B.LPAP (L).Phase := P + Q;
      end if;
   end Add_Z;

   procedure Perturb (D : in out Shared.Param_S; K : Positive; E : Long_Float) is
   begin
      case Kinds (K) is
         when Re => Add_Z (D, Lines (K), E, 0.0);
         when Im => Add_Z (D, Lines (K), 0.0, E);
         when Imp_A => D.B.ImpA := D.B.ImpA + E;
         when Imp_B => D.B.ImpB := D.B.ImpB + E;
      end case;
   end Perturb;

   procedure Apply (D : in out Shared.Param_S; V : Real_Vector; Alpha : Long_Float) is
   begin
      for K in 1 .. NP loop
         case Kinds (K) is
            when Re => Add_Z (D, Lines (K), Alpha * V (K), Alpha * V (K + 1));
            when Im => null;
            when Imp_A => D.B.ImpA := D.B.ImpA + Alpha * V (K);
            when Imp_B => D.B.ImpB := D.B.ImpB + Alpha * V (K);
         end case;
      end loop;
   end Apply;

   --  current value of Re/Im parameter K
   function Z_Value (D : Shared.Param_S; K : Positive) return Long_Float is
      A : constant Long_Float := D.B.LPAP (Lines (K)).Amplitude;
      P : constant Long_Float := D.B.LPAP (Lines (K)).Phase;
   begin
      return (if Kinds (K) = Re then A * Cos (P) else A * Sin (P));
   end Z_Value;

   --  mean (((z - z_anchor) / Scale)^2) over the anchored parameters
   function Anchor_P (D : Shared.Param_S) return Long_Float is
      S : Long_Float := 0.0;
   begin
      if N_Anchor = 0 then
         return 0.0;
      end if;
      for K in 1 .. NP loop
         if Has_Anchor (K) then
            S := S + ((Z_Value (D, K) - Z_Anchor (K)) / Scale (K)) ** 2;
         end if;
      end loop;
      return S / Long_Float (N_Anchor);
   end Anchor_P;

   --  the penalty in correlation units: a relative residual change dSSE/SSE0
   --  is about -2 CC0 dCC / (1 - CC0^2) near the seed score CC0
   function Penalty (D : Shared.Param_S) return Long_Float is
      C0 : constant Long_Float :=
        Long_Float'Min (Long_Float'Max (Seed_Score, 0.05), 0.99);
   begin
      return (if N_Anchor = 0 then 0.0
              else Anchor_W * Anchor_P (D) * (1.0 - C0 ** 2) / (2.0 * C0));
   end Penalty;

   procedure Read_Anchor is
      F : Text_IO.File_Type;
      Per, X, Y : Long_Float;
   begin
      Text_IO.Open (F, Text_IO.In_File, Anchor_File);
      while not Text_IO.End_Of_File (F) loop
         Ada.Long_Float_Text_IO.Get (F, Per);
         Ada.Long_Float_Text_IO.Get (F, X);
         Ada.Long_Float_Text_IO.Get (F, Y);
         if not Text_IO.End_Of_File (F) then
            Text_IO.Skip_Line (F);
         end if;
         for K in 1 .. NP loop
            if Kinds (K) in Re .. Im and then
              abs (Base.A.LP (Lines (K)) - Per) <= 1.0E-4 * abs Per
            then
               Z_Anchor (K) := (if Kinds (K) = Re then X else Y);
               if not Has_Anchor (K) then
                  Has_Anchor (K) := True;
                  N_Anchor := N_Anchor + 1;
               end if;
            end if;
         end loop;
      end loop;
      Text_IO.Close (F);
      Text_IO.Put ("INVERSION anchor " & Anchor_File & ":" &
                   Integer'Image (N_Anchor / 2) & " lines, weight");
      Put_F (Anchor_W);
      Text_IO.Put (", start P");
      Put_F (Anchor_P (Base.all), 5);
      Text_IO.New_Line;
      for K in 1 .. NP loop
         if Kinds (K) = Re and then not Has_Anchor (K) then
            Text_IO.Put ("INVERSION anchor: no row for the line of period");
            Put_F (Base.A.LP (Lines (K)), 6);
            Text_IO.New_Line;
         end if;
      end loop;
   end Read_Anchor;

   procedure Setup_Params (D : Shared.Param_S) is
      Amps : Real_Vector (D.B.LPAP'Range);
      NA : Natural := 0;
      Med : Long_Float := 1.0;
   begin
      for L in D.B.LPAP'Range loop
         if D.B.LPAP (L).Amplitude /= 0.0 then
            NA := NA + 1;
            Amps (NA) := abs D.B.LPAP (L).Amplitude;
         end if;
      end loop;
      if NA > 0 then  -- median of the nonzero amplitudes (insertion sort)
         for I in 2 .. NA loop
            declare
               X : constant Long_Float := Amps (I);
               K : Integer := I - 1;
            begin
               while K >= 1 and then Amps (K) > X loop
                  Amps (K + 1) := Amps (K);
                  K := K - 1;
               end loop;
               Amps (K + 1) := X;
            end;
         end loop;
         Med := (if NA mod 2 = 1 then Amps ((NA + 1) / 2)
                 else 0.5 * (Amps (NA / 2) + Amps (NA / 2 + 1)));
      end if;
      for L in D.B.LPAP'Range loop
         if D.B.LPAP (L).Amplitude /= 0.0 then
            for K in Param_Kind range Re .. Im loop
               NP := NP + 1;
               Kinds (NP) := K;
               Lines (NP) := L;
               Scale (NP) := abs D.B.LPAP (L).Amplitude + 0.05 * Med;
               Eps (NP) := 0.01 * Scale (NP);
            end loop;
         end if;
      end loop;
      NP := NP + 1;
      Kinds (NP) := Imp_A;
      Scale (NP) := abs D.B.ImpA + 1.0;
      Eps (NP) := 0.01 * Scale (NP);
      NP := NP + 1;
      Kinds (NP) := Imp_B;
      Scale (NP) := abs D.B.ImpB + 1.0;
      Eps (NP) := 0.01 * Scale (NP);
      Base := new Shared.Param_S'(D);
      Text_IO.Put_Line
        ("INVERSION " & (if Method = GN then "GN" else "VITERBI") & ":" &
         Integer'Image ((NP - 2) / 2) & " lines + ImpA, ImpB;" &
         N_Legs'Img & " leg(s)" &
         (if Method /= Method_Asked then " (VITERBI not available with LEGS)"
          else ""));
      if Anchored then
         Read_Anchor;
      end if;
   end Setup_Params;

   procedure Setup_Leg (L : Positive; Train : Mask) is
      S : Leg_Store renames Legs (L);
   begin
      for I in Train'Range loop
         if Train (I) then
            S.NT := S.NT + 1;
         end if;
      end loop;
      S.Tr := new Index_List (1 .. S.NT);
      S.NT := 0;
      for I in Train'Range loop
         if Train (I) then
            S.NT := S.NT + 1;
            S.Tr (S.NT) := I;
         end if;
      end loop;
      S.M0 := new Real_Vector (1 .. S.NT);
      S.F0 := new Real_Vector (1 .. S.NT);
      S.Target := new Real_Vector (1 .. S.NT);
   end Setup_Leg;

   --  Viterbi path of the manifold correction over all samples (prior only
   --  where not training), the regression held at the base evaluation's
   --  coefficients; the model change at a sample is approximated by the
   --  filter's center weight times the change of its own unlagged value.
   procedure Viterbi_Target
     (Model, Data, Forcing : Data_Pairs; Train : Mask;
      Winding_Period, Center_Weight : Long_Float;
      Pointwise : not null access function
        (I : Integer; F : Long_Float) return Long_Float)
   is
      NT : constant Natural := Legs (1).NT;
      Tr : Index_List renames Legs (1).Tr.all;
      Target : Real_Vector renames Legs (1).Target.all;
      SP : constant Long_Float := Prior * Winding_Period;
      SS : constant Long_Float := Step_SD * Winding_Period;
      DD : constant Long_Float := SS / 6.0;
      H : constant Integer := Integer (Long_Float'Ceiling (3.0 * SP / DD));
      NG : constant Positive := 2 * H + 1;
      Span : constant Integer := Integer (Long_Float'Ceiling (4.0 * SS / DD));
      Lo : constant Integer := Model'First;
      Hi : constant Integer := Model'Last;
      type Back_Array is array (Lo .. Hi, 1 .. NG) of Integer;
      type Back_Ptr is access Back_Array;
      Back : constant Back_Ptr := new Back_Array;
      Cost, New_Cost, Loc : Real_Vector (1 .. NG);
      Path : array (Lo .. Hi) of Integer;
      S2 : Long_Float := 0.0;
      function Grid (G : Integer) return Long_Float is (Long_Float (G - 1 - H) * DD);
      procedure Local (I : Integer) is
         B0, R : Long_Float;
      begin
         for G in 1 .. NG loop
            Loc (G) := (Grid (G) / SP) ** 2;
         end loop;
         if Train (I) then
            B0 := Pointwise (I, Forcing (I).Value);
            R := Data (I).Value - Model (I).Value;
            for G in 1 .. NG loop
               Loc (G) := Loc (G) +
                 (R - Center_Weight *
                    (Pointwise (I, Forcing (I).Value + Grid (G)) - B0)) ** 2 / S2;
            end loop;
         end if;
      end Local;
   begin
      for R in 1 .. NT loop
         S2 := S2 + (Data (Tr (R)).Value - Model (Tr (R)).Value) ** 2;
      end loop;
      S2 := Long_Float'Max (S2 / Long_Float (NT), 1.0E-30);
      Local (Lo);
      Cost := Loc;
      for I in Lo + 1 .. Hi loop
         Local (I);
         for G in 1 .. NG loop
            declare
               Best : Long_Float := Long_Float'Last;
               Arg : Integer := G;
            begin
               for S in -Span .. Span loop
                  if G - S in 1 .. NG then
                     declare
                        V : constant Long_Float :=
                          Cost (G - S) + (Long_Float (S) * DD / SS) ** 2;
                     begin
                        if V < Best then
                           Best := V;
                           Arg := G - S;
                        end if;
                     end;
                  end if;
               end loop;
               New_Cost (G) := Best + Loc (G);
               Back (I, G) := Arg;
            end;
         end loop;
         Cost := New_Cost;
      end loop;
      Path (Hi) := 1;
      for G in 2 .. NG loop
         if Cost (G) < Cost (Path (Hi)) then
            Path (Hi) := G;
         end if;
      end loop;
      for I in reverse Lo + 1 .. Hi loop
         Path (I - 1) := Back (I, Path (I));
      end loop;
      declare
         Sum, Sq : Long_Float := 0.0;
      begin
         for R in 1 .. NT loop
            Target (R) := Grid (Path (Tr (R)));
            Sum := Sum + Target (R);
            Sq := Sq + Target (R) ** 2;
         end loop;
         Text_IO.Put ("INVERSION Viterbi dF sd");
         Put_F (Sqrt (Long_Float'Max (Sq / Long_Float (NT) - (Sum / Long_Float (NT)) ** 2, 0.0)), 5);
         Text_IO.Put (" (fastest winding period");
         Put_F (Winding_Period, 5);
         Text_IO.Put_Line (")");
      end;
   end Viterbi_Target;

   procedure Record_Result
     (D : in Shared.Param_S; Model, Data, Forcing : in Data_Pairs;
      Train : in Mask; Score, Held : in Long_Float; Singular : in Boolean;
      Winding_Period, Center_Weight : in Long_Float;
      Pointwise : not null access function
        (I : Integer; F : Long_Float) return Long_Float;
      Leg : in Natural := 0; Leg_Weight : in Long_Float := 1.0;
      Valid, Test : in Long_Float := Long_Float'First)
   is
      L : constant Positive := Integer'Max (1, Leg);
      Lead : constant Boolean := L = 1;
      S : Leg_Store renames Legs (L);
   begin
      if S.Tr = null then
         Setup_Leg (L, Train);
      end if;
      if Lead and then Base = null then
         Setup_Params (D);
      end if;
      S.Was_Base := Phase = Want_Base;
      S.Was_Final := Final;
      S.Val := Valid;
      S.Test := Test;
      S.LW := Leg_Weight;
      case Phase is
         when Want_Base =>
            declare
               SS : Long_Float := 0.0;
            begin
               for R in 1 .. S.NT loop
                  S.M0 (R) := Model (S.Tr (R)).Value;
                  S.F0 (R) := Forcing (S.Tr (R)).Value;
                  S.Target (R) := Data (S.Tr (R)).Value - S.M0 (R);
                  SS := SS + S.Target (R) ** 2;
               end loop;
               S.Weight := Sqrt (Leg_Weight) /
                 Long_Float'Max (Sqrt (SS / Long_Float (Integer'Max (1, S.NT))), 1.0E-30);
            end;
            if Lead then
               Base.all := D;
               if Seed_Score = Missing then
                  Seed_Score := Score;
               end if;
               Base_Score := Score - Penalty (D);
               if Final then
                  Text_IO.Put ("INVERSION final = iteration" & Best_Val_Iter'Img &
                               " (best validation): train");
               elsif Iter = 0 then
                  Text_IO.Put ("INVERSION seed: train");
               else
                  Text_IO.Put ("INVERSION it" & Iter'Img & " accepted ridge");
                  Put_F (Cand_Ridge (Picked), 3);
                  Text_IO.Put (" step");
                  Put_F (Cand_Step (Picked), 2);
                  Text_IO.Put (": train");
               end if;
               Put_F (Score);
               Text_IO.Put (" held-out");
               Put_F (Held);
               Put_VT;
               if N_Anchor > 0 then
                  Text_IO.Put (" anchor P");
                  Put_F (Anchor_P (D), 5);
                  Text_IO.Put (" penalty");
                  Put_F (Penalty (D), 5);
               end if;
               Text_IO.New_Line;
               if not Final then
                  declare
                     V : constant Long_Float := Composite (False);
                  begin
                     if Best_Val_D = null or else V > Best_Val then
                        Best_Val := V;
                        Best_Val_Iter := Iter;
                        Since_Best := 0;
                        if Best_Val_D = null then
                           Best_Val_D := new Shared.Param_S'(Base.all);
                        else
                           Best_Val_D.all := Base.all;
                        end if;
                     else
                        Since_Best := Since_Best + 1;
                     end if;
                  end;
               end if;
               if Method = Viterbi and then Iter < Max_Iter and then not Final then
                  Viterbi_Target (Model, Data, Forcing, Train, Winding_Period,
                                  Center_Weight, Pointwise);
               end if;
            end if;
         when Probing =>
            if S.J = null then
               S.J := new Real_Matrix (1 .. S.NT, 1 .. NP);
            end if;
            for R in 1 .. S.NT loop
               S.J (R, Probe) :=
                 (if Singular then 0.0
                  elsif Method = GN then (Model (S.Tr (R)).Value - S.M0 (R)) / Eps (Probe)
                  else (Forcing (S.Tr (R)).Value - S.F0 (R)) / Eps (Probe));
            end loop;
         when Trying =>
            if Lead then
               Cand_Score (Cand) :=
                 (if Singular then Long_Float'First else Score - Penalty (D));
               Text_IO.Put ("INVERSION it" & Natural'Image (Iter + 1) & " ridge");
               Put_F (Cand_Ridge (Cand), 3);
               Text_IO.Put (" step");
               Put_F (Cand_Step (Cand), 2);
               Text_IO.Put (": linearized gain");
               Put_F (Cand_Gain (Cand));
               Text_IO.Put (" | train");
               Put_F (Score);
               Text_IO.Put (" held-out");
               Put_F (Held);
               Put_VT;
               if N_Anchor > 0 then
                  Text_IO.Put (" anchor P");
                  Put_F (Anchor_P (D), 5);
               end if;
               Text_IO.New_Line;
            end if;
      end case;
   end Record_Result;

   --  Ridge steps from the rows of every leg (each leg's rows weighted by its
   --  Weight; columns scaled by Scale). The normal equations are accumulated
   --  row by row, so no rows x columns matrix is formed.
   procedure Solve is
      G : Real_Matrix (1 .. NP, 1 .. NP) := (others => (others => 0.0));
      AB : Real_Vector (1 .. NP) := (others => 0.0);
      B2, Tr_G : Long_Float := 0.0;
      Row : Real_Vector (1 .. NP);
      K : Natural := 0;
   begin
      for L in 1 .. N_Legs loop
         declare
            S : Leg_Store renames Legs (L);
            W : constant Long_Float := S.Weight;
         begin
            for R in 1 .. S.NT loop
               for C in 1 .. NP loop
                  Row (C) := W * S.J (R, C) * Scale (C);
               end loop;
               for C1 in 1 .. NP loop
                  for C2 in C1 .. NP loop
                     G (C1, C2) := G (C1, C2) + Row (C1) * Row (C2);
                  end loop;
                  AB (C1) := AB (C1) + Row (C1) * W * S.Target (R);
               end loop;
               B2 := B2 + (W * S.Target (R)) ** 2;
            end loop;
         end;
      end loop;
      for C1 in 1 .. NP loop
         for C2 in 1 .. C1 - 1 loop
            G (C1, C2) := G (C2, C1);
         end loop;
         Tr_G := Tr_G + G (C1, C1);
      end loop;
      --  anchor rows: Anchor_W * B2 * mean ((z + dz - z_anchor) / Scale)^2,
      --  dz / Scale being the scaled unknowns
      if N_Anchor > 0 then
         declare
            WA : constant Long_Float := Anchor_W * B2 / Long_Float (N_Anchor);
         begin
            for K in 1 .. NP loop
               if Has_Anchor (K) then
                  G (K, K) := G (K, K) + WA;
                  AB (K) := AB (K) +
                    WA * (Z_Anchor (K) - Z_Value (Base.all, K)) / Scale (K);
               end if;
            end loop;
         end;
      end if;
      for Lam of Ridges loop
         declare
            V : constant Real_Vector :=
              Solve (G + (Lam * Tr_G / Long_Float (NP)) * Unit_Matrix (NP), AB);
            Dp : Real_Vector (1 .. NP);
            E2 : Long_Float := 0.0;
         begin
            for C in 1 .. NP loop
               Dp (C) := V (C) * Scale (C);
            end loop;
            for L in 1 .. N_Legs loop
               declare
                  S : Leg_Store renames Legs (L);
                  F : Long_Float;
               begin
                  for R in 1 .. S.NT loop
                     F := S.Target (R);
                     for C in 1 .. NP loop
                        F := F - S.J (R, C) * Dp (C);
                     end loop;
                     E2 := E2 + (S.Weight * F) ** 2;
                  end loop;
               end;
            end loop;
            for Alpha of Steps loop
               K := K + 1;
               if Cands (K) = null then
                  Cands (K) := new Shared.Param_S'(Base.all);
               else
                  Cands (K).all := Base.all;
               end if;
               Apply (Cands (K).all, Dp, Alpha);
               Cand_Ridge (K) := Lam;
               Cand_Step (K) := Alpha;
               Cand_Gain (K) := (if B2 > 0.0 then 1.0 - E2 / B2 else 0.0);
            end loop;
         end;
      end loop;
   end Solve;

   procedure Next (D : in out Shared.Param_S; Done : out Boolean) is
      --  Stop; with INVERSION_STOP VALID first go back to the best-validation
      --  base (one more evaluation, accepted and saved) unless it is current
      procedure Finish (Why : String) is
      begin
         Text_IO.Put_Line ("INVERSION: " & Why);
         if Stop_Valid and then Best_Val_D /= null and then Best_Val_Iter /= Iter then
            Final := True;
            Phase := Want_Base;
            D := Best_Val_D.all;
         else
            Done := True;
         end if;
      end Finish;
   begin
      Done := False;
      case Phase is
         when Want_Base =>
            if Final then
               Done := True;
               return;
            elsif Iter >= Max_Iter then
               Finish ("done," & Iter'Img & " iterations");
               return;
            elsif Stop_Valid and then Since_Best >= Patience then
               Finish ("no validation gain for" & Since_Best'Img &
                       " iterations (best at iteration" & Best_Val_Iter'Img & ")");
               return;
            end if;
            Phase := Probing;
            Probe := 1;
            D := Base.all;
            Perturb (D, Probe, Eps (Probe));
         when Probing =>
            if Probe < NP then
               Probe := Probe + 1;
               D := Base.all;
               Perturb (D, Probe, Eps (Probe));
            else
               Solve;
               Phase := Trying;
               Cand := 1;
               D := Cands (Cand).all;
            end if;
         when Trying =>
            if Cand < NC then
               Cand := Cand + 1;
               D := Cands (Cand).all;
            else
               Picked := 1;
               for C in 2 .. NC loop
                  if Cand_Score (C) > Cand_Score (Picked) then
                     Picked := C;
                  end if;
               end loop;
               if Cand_Score (Picked) > Base_Score + Min_Gain then
                  Iter := Iter + 1;
                  Phase := Want_Base;
                  D := Cands (Picked).all;
               else
                  Finish ("no training gain at iteration" &
                          Natural'Image (Iter + 1) & "; stopping");
               end if;
            end if;
      end case;
   end Next;

end GEM.LTE.Primitives.Inversion;
