package GEM.LTE.Primitives.Shared is

   Max_Harmonics : constant := 30;

   type Param_A (NLP, NLT : Integer) is record
      k0 : Long_Float;
      level : Long_Float;
      LP : Long_Periods (1 .. NLP);  -- Fixed
      LTAP : Modulations_Amp_Phase (1 .. NLT); -- Calculated
   end record;

   type Param_B (NLP, NLT : Integer)
   is record                    -- Random walk values
      Offset : Long_Float;   -- Integrated trend on forcing
      bg : Long_Float;   -- Background of impulse
      ImpA : Long_Float;   -- Amplitude of sin impulse
      ImpB : Long_Float;   -- Phase of sin impulse
      ImpC : Long_Float;   -- Power of sin impulse
      DelA : Long_Float;   -- Amp of delta impulse
      DelB : Long_Float;   -- Phase of delta impulse
      Asym : Long_Float;   -- Asymmetry of semi-annual impulse

      Ann1 : Long_Float;   -- Annual sine coefficient
      Ann2 : Long_Float;   -- Annual cosine coefficient
      Sem1 : Long_Float;   -- Semi-annual sine coefficient
      Sem2 : Long_Float;   -- Semi-annual cosine coefficient
      IR : Long_Float;   -- Impulse pass-through
      Year : Long_Float;   -- Year correction (in days)

      mA : Long_Float;   -- 1st order response
      mP : Long_Float;   -- 2nd order response
      shiftT : Long_Float;   -- Starting time correction
      init : Long_Float;   -- Initial value

      LPAP : Long_Periods_Amp_Phase (1 .. NLP);
      LT : Modulations (1 .. NLT);
   end record;
   
   --  ==========================================================================
   --  MEMORY LAYOUT DOCUMENTATION
   --  ==========================================================================
   --  Param_B is designed for contiguous memory layout to support array overlay
   --  in the random descent optimization algorithm.
   --
   --  STRUCTURE (for NLP=29, NLT=11):
   --    - 18 scalar Long_Float fields (Offset through init)
   --    - LPAP array: NLP Amp_Phase records = NLP * 2 Long_Floats (58 floats)
   --    - LT array: NLT Long_Floats (11 floats)
   --    TOTAL: 18 + 58 + 11 = 87 Long_Float values
   --
   --  SIZE CALCULATION:
   --    With GNAT compiler, unconstrained arrays have dope information (bounds)
   --    Record size = discriminants + scalar fields + array data
   --    For overlay: Size_Shared = Record'Size / Long_Float'Size - 1
   --      The -1 accounts for discriminant storage (NLP, NLT)
   --
   --  REPRESENTATION:
   --    No explicit representation clause is added because GNAT naturally
   --    packs arrays contiguously with no gaps. Adding explicit component
   --    clauses would require compile-time knowledge of NLP and NLT, which
   --    are runtime discriminants.
   --
   --  VERIFICATION:
   --    The Size_Shared calculation in gem-lte-primitives-solution.adb
   --    verifies at runtime that the layout is as expected.
   --  ==========================================================================

   -- Record structure used for optimization sharing
   -- This could use the benefit of Ada record representation clauses,
   -- but as it is with GNAT compiler, there is no dope information
   -- and so is essentially a linear array of long_floats

   type Param_S (NLP, NLT : Integer) is record
      A : Param_A (NLP, NLT);
      B : Param_B (NLP, NLT);
      C : Ns (1 .. Max_Harmonics) := (others => 0);
   end record;

   -- don't think this is necessary but watch for cases where the
   -- compiler may optimize away reads and treat the data as unused
   -- pragma Volatile(Param_S);

   --
   -- This is the approach for delivering the input parameters to the
   -- processing threads. The call to Put will release the Get per thread
   --
   procedure Put (P : in Param_S);

   function Get (N_Tides, N_Modulations : in Integer) return Param_S;

   --
   -- from a file -- name of file executable
   --
   procedure Save (P : in Param_S);

   procedure Load (P : in out Param_S);

   procedure Dump (D : in Param_S);

   --  Persist the winding amplitudes/phases (tidal + harmonic) computed by
   --  the regression at the end of a run, as a standalone JSON file
   --  alongside the existing .p/.par output -- named "<exe>.windings.json".
   --  Purely additive: mirrors the "---- LTE ----" stdout block, does not
   --  read back in anywhere or affect Save/Load/Dump or the optimization
   --  loop. Intended to let downstream tooling collect amplitude/phase
   --  data across many runs/indices without scraping stdout.
   --
   --  Ctx (optional) adds a "secular_context" block giving Trend/Accel the
   --  dates they depend on. The MLR (Regression_Factors) fits, and LTE
   --  evaluates, the secular part as
   --     Trend * t + Accel * (t - Accel_Ref)**2      (t = decimal year)
   --  with Accel_Ref = the first date of the array the MLR was handed.
   --  Under EXCLUDE that is the record start: TRAIN_START/TRAIN_END then
   --  bound the excluded TEST interval, not the training span.
   type Secular_Context is record
      Known          : Boolean    := False;
      Accel_Ref      : Long_Float := 0.0;
      Fit_Start      : Long_Float := 0.0;  -- first/last date regressed on
      Fit_End        : Long_Float := 0.0;
      Interval_Start : Long_Float := 0.0;  -- TRAIN_START/TRAIN_END dates
      Interval_End   : Long_Float := 0.0;
      Record_Start   : Long_Float := 0.0;
      Record_End     : Long_Float := 0.0;
      Exclude        : Boolean    := False;
      Enclosing      : Boolean    := False;
      Coverage       : Long_Float := 1.0;
      Accel_Enabled  : Boolean    := True;  -- ACCEL env/resp setting
      Aero_Enabled   : Boolean    := False; -- AERO file loaded
      Aero_Coef      : Long_Float := 0.0;   -- fitted aerosol coefficient
      Aero_Lag       : Long_Float := 0.0;   -- AERO_LAG, years
      Alpha          : Long_Float := 0.0;   -- ALPHA gain fraction (0 = off)
      M_Min          : Long_Float := 0.0;   -- manifold minimum used with it
   end record;
   No_Secular_Context : constant Secular_Context := (others => <>);

   procedure Save_Windings
     (Trend, Accel, K0, Level, IR : in Long_Float;
      M                           : in Modulations;
      MAP                         : in Modulations_Amp_Phase;
      NM, NH                      : in Integer;
      B                           : in Param_B;
      Ctx                         : in Secular_Context := No_Secular_Context);

end GEM.LTE.Primitives.Shared;
