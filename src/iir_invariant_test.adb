--  IIR_Invariant_Test -- standalone self-test for GEM.LTE.Primitives.IIR,
--  checking two invariants that a production model relies on implicitly
--  whenever the same manifold is evaluated on time series of different
--  lengths (e.g. kS040_W050, which starts in 1950, vs kS040_W050_, the
--  same underlying series extended back before 1900):
--
--  INVARIANT 1 (short/long consistency): if a shorter series and a
--  longer series agree on every real sample in their overlap, then
--  IIR's OWN OUTPUT over that overlap must also agree -- the correlation
--  coefficient computed over a fixed [TRAIN_START, TRAIN_END] window
--  must not depend on how much additional history precedes it. This is
--  only true if the recursion is seeded at the SAME true calendar
--  anchor in both cases. IIR's `Start_Index` search
--  ("first I where Raw(I).Date > Start") silently degrades to
--  Raw'First whenever Start (IDATE) precedes every date in the array --
--  which is exactly what happens for a short series whose own record
--  starts after IDATE. That silently re-seeds the recursion at the
--  short series' own first sample instead of truly at IDATE, breaking
--  the invariant. This test demonstrates the CURRENT behavior directly
--  (see NOTE below) rather than asserting a fix that hasn't been
--  applied to Calc_Forcing's call sites yet.
--
--  INVARIANT 2 (forward/backward round-trip) -- FIXED, this test
--  confirms it: IIR's backward pass ("if Start is inside the series,
--  create pre-history by running backwards") must be the exact
--  mathematical inverse of the forward recursion -- reconstructing an
--  earlier value from a later one must recover exactly what a genuine
--  forward pass from further back would have produced, given the SAME
--  driving signal and the SAME seed value at the shared point. The
--  backward branch used to be a plausible-looking but NOT algebraically
--  exact heuristic (this test originally caught a ~2.5-unit round-trip
--  discrepancy under it); IIR now uses the exact piecewise-closed-form
--  inverse instead (see gem-lte-primitives.adb's IIR for the derivation
--  and its two documented edge cases -- Mem=0 and the genuine
--  non-invertible dead zone of half-width lagC around zero). This test
--  checks the fix directly by comparing a pure long forward pass against
--  a mid-series-seeded pass's backward reconstruction, and gates the
--  overall PASS/FAIL result.
--
--  Invariant 1a is reported but does NOT gate the result -- it is the
--  DEFAULT (STRICT_IDATE absent/False) behavior, kept exactly as-is on
--  purpose for backward compatibility with every existing .resp/.p, so
--  this is a known, permanent characteristic of that mode, not a bug to
--  fix -- see gem-lte-primitives-solution.adb's Strict_IDate flag.
--
--  INVARIANT 3 (STRICT_IDATE's actual mechanism, opt-in): with
--  GEM.LTE.Primitives.Extend_Backward prepending synthetic (date-only)
--  rows before a short series so it reaches back to IDATE before Tide_
--  Sum/Impulse_Amplify/IIR ever run -- the mechanism Calc_Forcing now
--  uses when STRICT_IDATE=true -- short-vs-long consistency must hold
--  AUTOMATICALLY, with the raw config IDATE used as-is, no manual anchor
--  override needed (unlike Invariant 1's default-mode workaround in 1b).
--  Gates the overall result.
--
--  NOT covered by this test (needs the real executable, not this unit
--  test, since it calls IIR/Extend_Backward directly and never goes
--  through GEM.Getenv/resp parsing at all): two bugs found and fixed
--  while wiring STRICT_IDATE into Calc_Forcing, verified instead via
--  end-to-end lt.exe runs --
--    (a) Calc_Forcing must seed the extended computation at Init_Date
--        (the calendar date D.B.init's value actually applies at --
--        defaults to the CURRENT run's own Data_Records'First, which is
--        init's real, already-calibrated meaning for the file a region
--        was fit on), NOT at IDATE directly -- seeding at IDATE with the
--        unchanged init silently reinterprets what init means and
--        produces a materially worse fit (confirmed: CC 0.73 -> 0.57 on
--        kS040_W050 from this alone, with STRICT_IDATE otherwise
--        working correctly).
--    (b) Init_Date's own seeding needs the same "hair below the target"
--        treatment as everywhere else in this file: Start_Index's search
--        is for the first Date STRICTLY GREATER than Start, so passing
--        Init_Date unchanged lands one full sample LATE whenever it
--        exactly equals a real row's date (the common, default case) --
--        confirmed via the SAME symptom as (a) (0.73 -> 0.70, smaller
--        but still real) before being traced to this instead.
--    (c) A NEW resp key (INIT_DATE, to override Init_Date's default when
--        applying one region's manifold to a DIFFERENT, longer file than
--        it was fit on) must be added to gem.adb's Options enum --
--        Read_Response_File only stores keys matching that fixed
--        whitelist, so an unrecognized key is silently dropped, not an
--        error, and Getenv then silently falls through to its default --
--        caught by getting a suspiciously unchanged result before this
--        fix, not a build or runtime failure.
--  With all three fixed: kS040_W050 (F9=1, TRAIN_START=1955,
--  TRAIN_END=2020) and kS040_W050_ evaluated on the SAME already-fitted
--  manifold (no re-optimization) via STRICT_IDATE=true (INIT_DATE=1950.0
--  set explicitly for the long file) reproduce the SAME CC to every
--  displayed digit (0.7301248196099356) and max|diff|=0.0 on both the
--  model and forcing columns -- the deterministic, no-refit path the
--  short/long invariant was always meant to have.
--
--  Usage: run with no arguments. Exits 0 if the GATING checks (1b, 2, 3)
--  pass; exits 1 otherwise, printing every check's max discrepancy
--  either way. This is deliberately independent of -gnata (pragma
--  Assert) so it always runs regardless of build flags.

with Ada.Text_IO; use Ada.Text_IO;
with Ada.Numerics.Long_Elementary_Functions; use Ada.Numerics.Long_Elementary_Functions;
with GEM.LTE.Primitives; use GEM.LTE.Primitives;
with GNAT.OS_Lib;

procedure IIR_Invariant_Test is

   Samples_Per_Year : constant := 12;
   Lag_A : constant Long_Float := 0.99998;   -- mem, close to 1 like real fits (e.g. kS040_W050's ma)
   Lag_C : constant Long_Float := 0.0043;    -- like a real fitted mp
   Any_Failed : Boolean := False;

   --  Synthetic driving signal -- deterministic function of Date only, no
   --  real data involved, standing in for Tide_Sum(...)*Impulse_Delta(...).
   function Synthetic_Raw (Date : Long_Float) return Long_Float is
   begin
      return Sin (2.0 * Ada.Numerics.Pi * Date / 7.3) +
             0.3 * Cos (2.0 * Ada.Numerics.Pi * Date / 3.1) +
             0.1 * Sin (2.0 * Ada.Numerics.Pi * Date / 1.07);
   end Synthetic_Raw;

   --  Mirrors what Tide_Sum does to whatever template it's handed: fully
   --  recompute .Value from .Date alone, for every row -- including any
   --  synthetic rows Extend_Backward just prepended, which otherwise
   --  keep their Value=0.0 placeholder. Calc_Forcing always runs Tide_Sum
   --  immediately after Extend_Backward for exactly this reason.
   function Recompute_From_Dates (D : Data_Pairs) return Data_Pairs is
      Result : Data_Pairs := D;
   begin
      for I in Result'Range loop
         Result (I).Value := Synthetic_Raw (Result (I).Date);
      end loop;
      return Result;
   end Recompute_From_Dates;

   function Make_Series (Date_First, Date_Last : Long_Float) return Data_Pairs is
      N : constant Integer :=
        Integer (Long_Float'Rounding ((Date_Last - Date_First) *
                 Long_Float (Samples_Per_Year))) + 1;
      Result : Data_Pairs (1 .. N);
      D : Long_Float;
   begin
      for I in Result'Range loop
         D := Date_First + Long_Float (I - 1) / Long_Float (Samples_Per_Year);
         Result (I) := (Date => D, Value => Synthetic_Raw (D));
      end loop;
      return Result;
   end Make_Series;

   function Max_Abs_Diff_Over_Overlap
     (A, B : in Data_Pairs; From_Date : Long_Float) return Long_Float
   is
      Worst : Long_Float := 0.0;
      J : Integer := B'First;
      Diff : Long_Float;
   begin
      for I in A'Range loop
         if A (I).Date >= From_Date - 1.0e-9 then
            --  find matching row in B by nearest date (both on the same
            --  monthly grid, so this is an exact match, not an
            --  interpolation)
            while J < B'Last and then B (J).Date < A (I).Date - 1.0e-9 loop
               J := J + 1;
            end loop;
            if abs (B (J).Date - A (I).Date) < 1.0e-6 then
               Diff := abs (A (I).Value - B (J).Value);
               if Diff > Worst then
                  Worst := Diff;
               end if;
            end if;
         end if;
      end loop;
      return Worst;
   end Max_Abs_Diff_Over_Overlap;

   procedure Report (Name : String; Worst : Long_Float; Tolerance : Long_Float;
                      Gates_Result : Boolean := True) is
   begin
      if Worst <= Tolerance then
         Put_Line ("PASS  " & Name & "  max_diff=" & Worst'Img &
                    "  (tolerance " & Tolerance'Img & ")");
      else
         Put_Line ((if Gates_Result then "FAIL  " else "FAIL (known issue, non-gating)  ") &
                    Name & "  max_diff=" & Worst'Img &
                    "  EXCEEDS tolerance " & Tolerance'Img);
         if Gates_Result then
            Any_Failed := True;
         end if;
      end if;
   end Report;

   Long_First : constant Long_Float := 1856.0;
   Long_Last  : constant Long_Float := 2023.0;
   Short_First : constant Long_Float := 1950.0;
   IDate : constant Long_Float := 1880.0;   -- mirrors production's typical IDATE

begin
   Put_Line ("=== IIR_Invariant_Test ===");

   --------------------------------------------------------------------
   --  INVARIANT 1: short series vs long series, same overlap, same
   --  nominal IDATE -- must produce identical IIR output over the
   --  overlap. NOTE: this is EXPECTED TO FAIL under the current
   --  Start_Index clamping behavior whenever IDate < Short_First; that
   --  is precisely the bug this test exists to make visible and
   --  regression-proof, not a fix already applied to the call sites in
   --  gem-lte-primitives-solution.adb.
   --------------------------------------------------------------------
   declare
      Long_Series  : constant Data_Pairs := Make_Series (Long_First, Long_Last);
      Short_Series : constant Data_Pairs := Make_Series (Short_First, Long_Last);
      Long_Out  : constant Data_Pairs :=
        IIR (Long_Series, lagA => Lag_A, lagC => Lag_C, iA => 0.0, Start => IDate);
      Short_Out_As_Fit : constant Data_Pairs :=
        IIR (Short_Series, lagA => Lag_A, lagC => Lag_C, iA => 0.0, Start => IDate);
      Worst : constant Long_Float :=
        Max_Abs_Diff_Over_Overlap (Long_Out, Short_Out_As_Fit, Short_First);
   begin
      Put_Line ("[Invariant 1a, known issue, non-gating] naive: same nominal IDate (" &
                 IDate'Img & "), IDate precedes the short series' own start (" &
                 Short_First'Img & ")");
      Report ("short-vs-long IIR output over the overlap (naive IDate)",
              Worst, 1.0e-9, Gates_Result => False);

      --  The corrected anchor: NOTE this does NOT mean changing the
      --  short series' own Start -- Start_Index := "first I where
      --  Raw(I).Date > Start" already lands on Raw'First for ANY Start
      --  below the short series' own first date, so re-trying the short
      --  series alone with a different (still-too-early) Start changes
      --  nothing (Invariant 1a already IS that case). The actual fix is
      --  the other direction: evaluate the LONG series with Start
      --  forced to the training window's own first date too, instead of
      --  the true IDate -- a manual stand-in for what Calc_Forcing now
      --  does AUTOMATICALLY (and more correctly -- using the raw config
      --  IDATE, no manual override) when STRICT_IDATE=true; see
      --  Invariant 3 below for that actual mechanism.
      declare
         Anchor : constant Long_Float := Short_First - 1.0e-3;
         Long_Out_Anchored : constant Data_Pairs :=
           IIR (Long_Series, lagA => Lag_A, lagC => Lag_C, iA => 0.0, Start => Anchor);
         Worst2 : constant Long_Float :=
           Max_Abs_Diff_Over_Overlap (Long_Out_Anchored, Short_Out_As_Fit, Short_First);
      begin
         Put_Line ("[Invariant 1b] corrected: LONG series' own seed anchored " &
                    "at the training window's first date instead of the raw IDate config value");
         Report ("short-vs-long IIR output over the overlap (anchored seed)",
                 Worst2, 1.0e-9);
      end;
   end;

   --------------------------------------------------------------------
   --  INVARIANT 2: forward/backward round-trip. A pure long forward
   --  pass (seeded at the very first sample, so no backward
   --  reconstruction is exercised at all) must match a second pass that
   --  is seeded MID-SERIES (using the first pass's own value there) and
   --  reconstructs everything before that point via IIR's backward
   --  branch.
   --------------------------------------------------------------------
   declare
      Full_Series : constant Data_Pairs := Make_Series (Long_First, Long_Last);
      Forward_Only : constant Data_Pairs :=
        IIR (Full_Series, lagA => Lag_A, lagC => Lag_C, iA => 0.0,
             Start => Long_First - 1.0e-3);  -- seeds at the very first row
      Mid_Date : constant Long_Float := 1950.0;
      Seed_At_Mid : Long_Float := 0.0;
   begin
      for I in Forward_Only'Range loop
         if abs (Forward_Only (I).Date - Mid_Date) < 1.0e-6 then
            Seed_At_Mid := Forward_Only (I).Value;
         end if;
      end loop;

      declare
         --  Start_Index searches for the first Date STRICTLY GREATER
         --  than Start -- passing Mid_Date itself would seed one row
         --  LATE (at Mid_Date + 1 month), off by one from where
         --  Seed_At_Mid was actually read. Mirror the same
         --  hair-below-the-target convention used everywhere else in
         --  this test so Start_Index lands exactly on Mid_Date.
         Seeded_Mid : constant Data_Pairs :=
           IIR (Full_Series, lagA => Lag_A, lagC => Lag_C, iA => Seed_At_Mid,
                Start => Mid_Date - 1.0e-3);
         Worst : Long_Float := 0.0;
         Diff : Long_Float;
      begin
         for I in Full_Series'Range loop
            if Full_Series (I).Date < Mid_Date - 1.0e-9 then
               Diff := abs (Forward_Only (I).Value - Seeded_Mid (I).Value);
               if Diff > Worst then
                  Worst := Diff;
               end if;
            end if;
         end loop;
         Put_Line ("[Invariant 2] backward reconstruction (from " &
                    Mid_Date'Img & ") vs the true earlier forward-computed values");
         Report ("backward IIR reconstruction == forward integration", Worst, 1.0e-9);
      end;
   end;

   --------------------------------------------------------------------
   --  INVARIANT 3: STRICT_IDATE's actual mechanism (Extend_Backward),
   --  used automatically now by Calc_Forcing -- short-vs-long consistency
   --  with the RAW config IDate, no manual anchor override, unlike
   --  Invariant 1b's stand-in. Recompute_From_Dates mirrors Calc_Forcing's
   --  own Tide_Sum step, which always runs immediately after
   --  Extend_Backward and overwrites every row's Value from its Date --
   --  omitting it here would leave the newly-synthesized rows at their
   --  Value=0.0 placeholder and fail this test for a reason that has
   --  nothing to do with Extend_Backward or IIR themselves (caught during
   --  development: an earlier version of this test did exactly that and
   --  showed a spurious ~20-unit mismatch).
   --------------------------------------------------------------------
   declare
      Long_Series  : constant Data_Pairs := Make_Series (Long_First, Long_Last);
      Short_Series : constant Data_Pairs := Make_Series (Short_First, Long_Last);
      Extended_Short : constant Data_Pairs :=
        Recompute_From_Dates
          (Extend_Backward (Short_Series, IDate, Long_Float (Samples_Per_Year)));
      Long_Out : constant Data_Pairs :=
        IIR (Long_Series, lagA => Lag_A, lagC => Lag_C, iA => 0.0, Start => IDate);
      Short_Out_Extended : constant Data_Pairs :=
        IIR (Extended_Short, lagA => Lag_A, lagC => Lag_C, iA => 0.0, Start => IDate);
      Worst : constant Long_Float :=
        Max_Abs_Diff_Over_Overlap (Long_Out, Short_Out_Extended, Short_First);
   begin
      Put_Line ("[Invariant 3] Extend_Backward + IIR, raw IDate (" & IDate'Img &
                 "), no manual override");
      Report ("short(extended)-vs-long IIR output over the overlap", Worst, 1.0e-9);
   end;

   Put_Line ("===========================");
   if Any_Failed then
      Put_Line ("RESULT: at least one GATING invariant FAILED -- see above.");
      GNAT.OS_Lib.OS_Exit (1);
   else
      Put_Line ("RESULT: all gating invariants passed " &
                 "(1a's known, non-gating issue is unrelated -- see above).");
      GNAT.OS_Lib.OS_Exit (0);
   end if;
end IIR_Invariant_Test;
