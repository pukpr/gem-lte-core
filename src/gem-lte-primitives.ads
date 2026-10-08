package GEM.LTE.Primitives is

   type Pair is  -- A pair of values, such as in monthly temperature time-series
   record
      Date, Value : Long_Float;
   end record;

   type Data_Pairs is array (Integer range <>) of Pair;
   Empty_Data : constant Data_Pairs (1 .. 0) := (1 .. 0 => (0.0, 0.0));

   -- corresponding to a climate index such as ENSO
   function Make_Data (Name : in String) return Data_Pairs;

   -- Save results
   procedure Save
     (Model, Data, Forcing : in Data_Pairs;
      File_Name : in String := "lte_results.csv";
      IR : in Long_Float := 0.0);

   --
   -- Main algorithms
   --

   -- Infinite Impulse Response -- integrator
   function IIR
     (Raw : in Data_Pairs; lagA,  lagC : in Long_Float; -- lagB,
      iA : in Long_Float := 0.0; -- , iB, iC
      Start : in Long_Float := Long_Float'First) -- ; mA, mB : in Long_Float := 0.0)
      return Data_Pairs;

   --  Prepends synthetic (date-only, Value=0.0) rows before D'First, from
   --  D's own first date back to (at least) Target_First_Date, spaced at
   --  1/Sampling_Per_Year. A pure function of dates only -- for use
   --  ahead of a chain of equally pure-date functions (Tide_Sum, then an
   --  impulse gate, then IIR) so that chain can genuinely integrate from
   --  an early anchor date (e.g. IDATE) even when the real record does
   --  not reach back that far, rather than relying on IIR's own
   --  Start_Index search silently clamping to D'First (see IIR's own
   --  comment). Returns D unchanged (same bounds, no copy needed by the
   --  caller to detect this) if D already starts at or before
   --  Target_First_Date.
   function Extend_Backward
     (D : in Data_Pairs; Target_First_Date : in Long_Float;
      Sampling_Per_Year : in Long_Float) return Data_Pairs;

   -- Finite Impulse Response -- smoother
   function FIR
     (Raw : in Data_Pairs; Behind, Current, Ahead : in Long_Float)
      return Data_Pairs;

   -- Impulse modulation
   generic
      with function Impulse (Time : in Long_Float) return Long_Float;
   function Amplify
     (Raw : in Data_Pairs; Offset, Ramp, Start : in Long_Float)
      return Data_Pairs;

   -- LTE models = Superposition of tides + Laplace's Tidal Eqn modulation

   function Tide_Sum
     (Template : in Data_Pairs; Constituents : in Long_Periods_Amp_Phase;
      Periods : in Long_Periods; Ref_Time : in Long_Float := 0.0;
      Scaling : in Long_Float := 1.0; Cos_Phase : in Boolean := True;
      Year_Len : in Long_Float := Year_Length; Integ : in Long_Float := 0.0;
      Lod : in Long_Float := 0.0; Rho : in Long_Float := 1.0) return Data_Pairs;

   --  Accel_Ref (default Long_Float'First, a sentinel meaning "use
   --  Forcing'First's own date", today's exact behavior): the calendar
   --  date the Accel*(Date-Ref)**2 term is anchored at. Regression_
   --  Factors fits Accel assuming the reference is ITS OWN (possibly
   --  TRAIN_START-restricted) first date -- when LTE is then called to
   --  build Model over a DIFFERENT, larger array (e.g. a longer
   --  companion series under STRICT_IDATE/INIT_DATE), Forcing'First's
   --  date no longer matches that training reference, and the fitted
   --  Accel gets applied around the wrong zero-point -- confirmed
   --  directly: a 94-year reference mismatch (1856 vs the true 1950
   --  training reference) produced a ~2.3-unit Model discrepancy on
   --  kS040_W050 vs kS040_W050_ despite identical Forcing and nearly
   --  identical regression coefficients. Pass the SAME date
   --  Regression_Factors actually used (Forcing(First).Date at the
   --  call site, using that call's own TRAIN_START-resolved First) to
   --  fix this whenever Model is built over an array other than the one
   --  Accel was fit on.
   --  Quantity the ALPHA window is tested on at index I: the manifold
   --  value, or (ALPHA_SLOPE > 0) its signed slope; see the body.
   function Alpha_Signal (F : Data_Pairs; I : Integer) return Long_Float;

   --  X minus its seasonal pattern (monthly climatology minus its mean);
   --  used for the ANNUAL_DITHER anomaly comparison.
   function Remove_Climatology (X : Data_Pairs) return Data_Pairs;

   function LTE
     (Forcing : in Data_Pairs; Wave_Numbers : in Modulations;
      Amp_Phase : in Modulations_Amp_Phase;
      Offset, K0, Trend, Accel : in Long_Float := 0.0;
      NonLin : in Long_Float := 1.0;
      Annual : in Annual_Harmonics := (0.0, 0.0, 0.0, 0.0);
      Third : in Long_Float := 0.0;
      Accel_Ref : in Long_Float := Long_Float'First;
      Aero : in Long_Float := 0.0;
      --  ALPHA piecewise gain: each winding term is doubled while the
      --  manifold F is within an Alpha fraction of its minimum negative
      --  excursion (F - M_Min <= Alpha*|M_Min|), unchanged otherwise.
      --  Alpha = 0.0 (default) = no effect.
      --  Regression_Factors must be given the SAME Alpha and M_Min.
      Alpha : in Long_Float := 0.0;
      M_Min : in Long_Float := 0.0)
      return Data_Pairs;

   --  AERO (file path, default unset = off): volcanic aerosol series, two
   --  columns "date value" (e.g. stratospheric aerosol optical depth, or
   --  a sparse eruption profile, per the CSALT model). When loaded,
   --  Regression_Factors fits one extra coefficient (Aero) on
   --  Aero_Value(t), and LTE adds Aero * Aero_Value(t) to the model.
   --  AERO_LAG (years, default 0.0) shifts the series later in time:
   --  Aero_Value(t) = series(t - AERO_LAG), linearly interpolated, and 0.0
   --  outside the file's date range.
   function Aero_On return Boolean;
   function Aero_Value (Date : in Long_Float) return Long_Float;
   function Aero_Lag return Long_Float;

   -- Query to determine if a Tidal Constituent value should not be changed
   -- This uses a float comparison and is really only used for tidal periods
   function Is_Fixed (Value : in Long_Float) return Boolean;

   procedure Regression_Factors
     (Data_Records : in Data_Pairs;  -- Time series
   --First, Last,  -- Training Interval

      NM : in Positive; -- # modulations
      Forcing : in Data_Pairs;  -- Value @ Time
   -- Factors_Matrix : in out Matrix;

      DBLT : in Periods; DALTAP : out Amp_Phases; DALEVEL : out Long_Float;
      DAK0 : out Long_Float; Secular_Trend : in out Long_Float;
      Accel : out Long_Float;
      Accel_Ref : out Long_Float; -- date the Accel column is centred on
      Aero : out Long_Float;      -- aerosol coefficient (0 if AERO unset)
      Singular : out Boolean;
      Annual : out Annual_Harmonics;
      Third : in Long_Float := 0.0;
      IR : in Long_Float := 0.0;
      --  NONLIN power p: with p /= 1.0 each winding's two basis columns
      --  are sign(sin)|sin|**p and sign(cos)|cos|**p, the same waveform
      --  LTE evaluates (see LTE), so the fitted amplitudes are least-
      --  squares optimal for it. p = 1.0 is the ordinary sin/cos basis.
      NonLin : in Long_Float := 1.0;
      Alpha : in Long_Float := 0.0;   -- see LTE
      M_Min : in Long_Float := 0.0);

   --
   -- Utility procedures
   --
   
   function Winding_Agreement
     (Model, Data, Manifold : in Data_Pairs) return Long_Float;

   -- Correlation coefficient
   function CC (X, Y : in Data_Pairs) return Long_Float;

   -- Complexity-Invariant Distance
   function CID (X, Y : in Data_Pairs) return Long_Float;

   -- Zero-crossing metric
   function Xing (X, Y : in Data_Pairs) return Long_Float;

   -- Derivative CC
   function DER_CC (Model, Data : in Data_Pairs) return Long_Float;

   -- Dynamic Time Warp, Sakoe_Chiba_Optimized
   function DTW_Distance
     (X, Y : in Data_Pairs; Window_Size : Positive) return Long_Float;

   -- Earth-Mover Distance metric/Wasserstein
   function EMD
     (X, Y : in Data_Pairs; Derivative : in Boolean := False)
      return Long_Float;

   -- RMS
   function RMS
     (X, Y : in Data_Pairs; Ref, Offset : in Long_Float) return Long_Float;
   -- similar to RMS
   function Scaled_Error_Metric (X, Y : in Data_Pairs) return Long_Float;

   -- Correlation coefficient on spectrum
   function FT_CC (Model, Data, Forcing : in Data_Pairs) return Long_Float;

   -- Hoyer Peakedness
   function Hoyer_Spectral_Peak (Model, Data, Forcing : in Data_Pairs) return Long_Float;

   -- Minimum Entropy
   function Min_Entropy_Power_Spectrum
     (X, Y : in Data_Pairs) return Long_Float;

   function Is_Minimum_Entropy return Boolean;

   -- Dumps to stdIO all the data up to time corresponding to run_time
   procedure Dump
     (Model, Data : in Data_Pairs; Run_Time : Long_Float := 200.0);

   -- Halts the running threads
   procedure Stop;
   function Halted return Boolean;
   procedure Continue;

   -- 3 point median
   function Median (Raw : in Data_Pairs) return Data_Pairs;

   -- rectangular window of width = 2*Lobe_Width+1
   function Window
     (Raw : in Data_Pairs; Lobe_Width : in Positive) return Data_Pairs;

   NL : constant Boolean := True;
   procedure Put
     (Value : in Long_Float; Text : in String := "";
      New_Line : in Boolean := False);

   -- 9 point centered filter
   function Filter9Point (Raw : in Data_Pairs) return Data_Pairs;

end GEM.LTE.Primitives;
