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
      ImpA   : Long_Float;   -- Amplitude of sin impulse
      ImpB   : Long_Float;   -- Phase of sin impulse
      DelA   : Long_Float;   -- Amplitude of delta impulse
      DelB   : Long_Float;   -- Phase of delta impulse
      Asym   : Long_Float;   -- Asymmetry of semi-annual impulse
      Ann1   : Long_Float;   -- Annual cosine amplitude
      Ann2   : Long_Float;   -- Annual cosine phase
      Sem1   : Long_Float;   -- Semi-annual cosine amplitude
      Sem2   : Long_Float;   -- Semi-annual cosine phase
      mA     : Long_Float;   -- 1st order IIR feedback
      mP     : Long_Float;   -- 2nd order IIR feedback
      shiftT : Long_Float;   -- Starting time correction
      init   : Long_Float;   -- IIR initial value

      LPAP   : Long_Periods_Amp_Phase (1 .. NLP);
      LT     : Modulations (1 .. NLT);
   end record;
   
   --  ==========================================================================
   --  MEMORY LAYOUT DOCUMENTATION
   --  ==========================================================================
   --  Param_B is designed for contiguous memory layout to support array overlay
   --  in the random descent optimization algorithm.
   --
   --  STRUCTURE (for NLP=29, NLT=11):
   --    - 14 scalar Long_Float fields (Offset through init)
   --    - LPAP array: NLP Amp_Phase records = NLP * 2 Long_Floats (58 floats)
   --    - LT array: NLT Long_Floats (11 floats)
   --    TOTAL: 14 + NLP*2 + NLT = 8 + 58 + 11 = 77 Long_Float values
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

end GEM.LTE.Primitives.Shared;
