--  ============================================================================
--  GEM.LTE.Primitives.Solution - Parallel Optimization Solver
--  ============================================================================
--
--  PURPOSE:
--    Implements a multi-threaded parameter optimization engine for climate
--    modeling using GEM-LTE (Laplace's Tidal Equation). Spawns worker threads
--    that search parameter space in parallel, with a monitor for inter-thread
--    coordination and best-result tracking.
--
--  ARCHITECTURE:
--    - Thread-based parallelism: Maps worker tasks to available CPU cores
--    - Protected Monitor: Thread-safe coordination of best metric tracking
--    - Random Descent: Stochastic optimization algorithm for parameter search
--    - Dipole Model: Core climate modeling procedure (ENSO, IOD, etc.)
--
--  OPTIMIZATION STRATEGY:
--    Each thread performs random descent search in parameter space, testing
--    different tidal forcing configurations and LTE coefficients. The Monitor
--    tracks the globally best metric (typically correlation coefficient) and
--    allows threads to coordinate without explicit locking. When a thread
--    finds a better solution, it updates the Monitor, which triggers other
--    threads to potentially adjust their search strategy.
--
--  KEY ALGORITHMS:
--    1. CompareRef: Validates model against reference tidal data (dlod)
--    2. Thread Task: Worker thread mapped to specific CPU core
--    3. Monitor Protected Object: Thread-safe best-metric tracking
--    4. Dipole_Model: Main climate modeling computation with parameter fitting
--
--  CONFIGURATION:
--    Extensive environment variable configuration system (via GEM.Getenv)
--    controls all aspects: metrics (CC/RMS/DTW/EMD), training intervals,
--    optimization thresholds, LTE parameters, filtering, etc.
--
--  ============================================================================

with Text_IO;
with Ada.Numerics.Long_Elementary_Functions;
with Ada.Long_Float_Text_IO;
with GEM.Random_Descent;
with GNAT.Ctrl_C;
with System.Task_Info;
with System.Multiprocessors;
with GEM.LTE.Primitives.Shared;
with GEM.LTE.Primitives.Param_B_Overlay;
with Ada.Exceptions;
with GNAT.Traceback.Symbolic;
with Ada.IO_Exceptions;
with GNAT.OS_Lib;

package body GEM.LTE.Primitives.Solution is

   Is_Split : constant Boolean := GEM.Getenv ("SPLIT_TRAINING", False);
   Split_Low : constant Boolean := GEM.Getenv ("SPLIT_LOW", True);
   Alternate : constant Boolean := GEM.Getenv ("ALTERNATE", False);
   Monotonic_Increase : constant Boolean :=
     GEM.Getenv ("SECULAR", True); -- FALSE
   Trigger : Long_Float := GEM.Getenv ("TRIGGER", 0.999_999);
   Ratio : constant Long_Float := GEM.Getenv ("RATIO", 0.0);
   NLoops : constant Integer := GEM.Getenv ("NLOOPS", 20); --100
   Test_Only : constant Boolean := GEM.Getenv ("TEST_ONLY", False);

   --  Compare model output against reference tidal data (dLOD - day length)
   --  Used for validation rather than optimization - generates CSV files
   --  showing how well the tidal reconstruction matches reference data.
   function CompareRef
     (LP : in Long_Periods; LPRef, AP : in Long_Periods_Amp_Phase)
      return Long_Float
   is
      dLOD_dat : String := GEM.Getenv ("DLOD_DAT", "../dlod3.dat");
      D : Data_Pairs := Make_Data (dLOD_dat);
      Ref, R2, M2 : Data_Pairs := D;
      Metric : Long_Float := -1.0;
      Ref_Time : Long_Float := 0.0;
      File : Text_IO.File_Type;
   begin
      R2 :=
        Tide_Sum
          (Template => Ref, Constituents => LPRef, Periods => LP,
           Ref_Time => Ref_Time, Scaling => 1.0, Cos_Phase => False);

      M2 :=
        Tide_Sum
          (Template => Ref, Constituents => AP, Periods => LP,
           Ref_Time => Ref_Time, Scaling => 1.0, Cos_Phase => False);
      Metric := Long_Float'Max (CC (R2, M2), Metric);
      Text_IO.Create (File, Text_IO.Out_File, "dlod_compare.csv");
      for I in R2'Range loop
         Text_IO.Put_Line
           (File,
            R2 (I).Date'Img & ", " & R2 (I).Value'Img & ", " &
            M2 (I).Value'Img);
      end loop;
      Text_IO.Close (File);
      Text_IO.Create (File, Text_IO.Out_File, "dlod_ref.dat");
      for I in R2'Range loop
         Text_IO.Put_Line
           (File, R2 (I).Date'Img & ASCII.HT & M2 (I).Value'Img);
      end loop;
      Text_IO.Close (File);
      return Metric;
   exception
      when Constraint_Error =>
         Text_IO.Put_Line ("No reference compare");
         return 0.0;
   end CompareRef;

   --  Worker task type: Maps to specific CPU core for true parallelism
   --  Each thread performs independent parameter search with periodic
   --  coordination via the Monitor. Uses large stack (1GB) for deep
   --  recursive calculations and large data arrays.
   task type Thread
     (CPU : System.Multiprocessors.CPU; N_Tides, N_Modulations : Integer) is
      pragma CPU (CPU);
      pragma Storage_Size (1_000_000_000);
      entry Start;
   end Thread;
   type Thread_Access is access Thread;

   --  Spawns N worker threads and maps them round-robin to available CPU cores
   --  Installs Ctrl+C handler for clean shutdown. Each thread gets assigned
   --  a CPU and starts in deferred mode (waits for Start entry call).
   procedure Start
     (N_Tides, N_Modulations : in Integer; Number_of_Threads : Positive := 1)
   is
      TA : array (1 .. Number_of_Threads) of Thread_Access;
      use System.Multiprocessors;
      Num : constant Positive := System.Task_Info.Number_Of_Processors;
      CPU : Positive;
   begin
      GNAT.Ctrl_C.Install_Handler (Handler => GEM.LTE.Primitives.Stop'Access);
      for I in 1 .. Number_of_Threads loop
         CPU := I mod Num + 1;
         TA (I) :=
           new Thread
             (CPU => CPU_Range (CPU), N_Tides => N_Tides,
              N_Modulations => N_Modulations);
         TA (I).Start;
         delay 0.1; -- let them gradually start up
      end loop;
   end Start;

   Worst_Case : constant Long_Float := GEM.Getenv ("STARTING_METRIC", 0.001);

   --  =========================================================================
   --  Monitor: Protected object for inter-thread coordination
   --
   --  PURPOSE:
   --    Provides thread-safe tracking of the globally best metric found by
   --    any thread. Threads report their results via Check(), and can query
   --    current best via Status(). Uses entry barrier to efficiently notify
   --    waiting threads only when state changes (avoiding busy-wait).
   --
   --  COORDINATION STRATEGY:
   --    - Threads compete to find best metric (typically correlation coeff)
   --    - Only metric improvements trigger Status entry release
   --    - Percentage calculation allows threads to gauge their relative quality
   --    - Optional trigger threshold (TRIGGER env var) can halt all threads
   --      when target metric reached
   --  =========================================================================
   protected Monitor is
      procedure Check
        (Metric : in Long_Float;
         OOB : in Long_Float; -- Out-of-band
         Client : in Integer;
         Count : in Long_Integer;
         Best : out Boolean;  -- accessing thread deemed best
         BestClient : out Integer;
         Percentage : out Integer); -- % of best metric if not best
      procedure Client
        (ID : out Integer); -- registering a thread ID, called once
      entry Status
        (Metric : out Long_Float;   -- used by a monitoring thread, i.e. main
         OOB : out Long_Float;
         Client : out Integer;      -- returns client thread w/ best metric
         Cycle : out Long_Integer); -- and the cycle count it is on
      procedure Stop;
      procedure Reset;
   private
      --  Current best metric and associated metadata
      Best_Metric : Long_Float :=
        Worst_Case; -- so doesn't cause overflow for %
      Best_OOB : Long_Float := Worst_Case;
      Client_Index : Integer := 1;
      Waiting : Boolean := True; -- triggers status only if value changes
      Best_Client : Integer := -1;
      Best_Count : Long_Integer := 0;
   end Monitor;

   protected body Monitor is
      --  Check if this thread's metric is the new best
      --  Implements trigger-based early termination if metric exceeds threshold
      procedure Check
        (Metric : in Long_Float;
         OOB : in Long_Float; -- Out-of-band
         Client : in Integer;
         Count : in Long_Integer; Best : out Boolean; BestClient : out Integer;
         Percentage : out Integer)
      is
         M : Long_Float := Metric + Ratio * OOB;
      begin
         if M >= Best_Metric then
            Waiting := not (M > Best_Metric);
            Best_Metric := M;
            Best_OOB := OOB;
            Best := True;
            Best_Client := Client;
            Best_Count := Count;
            if Trigger > 0.0 and Best_OOB > Trigger then
               Primitives.Stop;
            end if;
         else
            Best := False;
         end if;
         if Best_Metric > Worst_Case then
            Percentage := Integer (100.0 * M / Best_Metric);
         else
            Percentage := 0;
         end if;
         BestClient := Best_Client;
      exception
         when Constraint_Error =>
            Best := False;
            BestClient := Best_Client;
            Percentage := 0;
      end Check;

      --  Assign unique ID to registering thread
      procedure Client (ID : out Integer) is
      begin
         ID := Client_Index;
         Client_Index := Client_Index + 1;
      end Client;

      --  Blocking call - only releases when best metric improves
      --  Used by monitoring thread (main) to track progress without polling
      entry Status
        (Metric : out Long_Float;
         OOB : out Long_Float; -- Out-of-band
         Client : out Integer;
         Cycle : out Long_Integer)
        when not Waiting is
      begin
         Waiting := True;
         Metric := Best_Metric;
         OOB := Best_OOB;
         Client := Best_Client;
         Cycle := Best_Count;
      end Status;

      --  Force Status entry to release (for clean shutdown)
      procedure Stop is
      begin
         Waiting :=
           False; -- necesary to allow a clean exit when program halted
      end Stop;

      --  Reset best metric tracking (used in alternating exclude mode)
      procedure Reset is
      begin
         Best_Metric := 0.0;
         Best_OOB := 0.0;
      end Reset;

   end Monitor;

   --  Query current best metric (blocking call - waits for improvement)
   --  Returns formatted string with thread ID, iteration count, metrics
   function Status return String is
      Metric, OOB : Long_Float;
      Client : Integer;
      Cycle : Long_Integer;
      S1, S2 : String (1 .. 10);
   begin
      Monitor.Status (Metric, OOB, Client, Cycle);
      Ada.Long_Float_Text_IO.Put (S1, Metric, Aft => 5, Exp => 0);
      Ada.Long_Float_Text_IO.Put (S2, OOB, Aft => 5, Exp => 0);
      return "Status: " & Client'Img & Cycle'Img & S1 & S2;
   end Status;

   -----------------------------------
   --  Thread body: Continuously runs Dipole_Model optimization
   --  Optional ALTERNATE mode switches between include/exclude training
   --  data on each iteration (for validation testing)
   -----------------------------------
   task body Thread is
      ID : Integer;
      Name : String :=
        GEM.Getenv (Name => "CLIMATE_INDEX", Default => "nino4.dat");
      Split : Boolean :=
        GEM.Getenv (Name => "SPLIT_TRAINING", Default => False);
      Exclude : Boolean := GEM.Getenv ("EXCLUDE", False);
   begin
      Monitor.Client (ID);
      Text_IO.Put_Line (Name & " for Thread #" & ID'Img & Thread.CPU'Img);
      accept Start;
      loop
         Dipole_Model
           (N_Tides => Thread.N_Tides, N_Modulations => Thread.N_Modulations,
            ID => ID, File_Name => Name, Split_Training => Split,
            Exclude => Exclude);
         exit when not Alternate;
         Exclude := not Exclude;
         Monitor.Reset;
         delay 1.0 + 1.0 * Duration (ID);
         Continue;

      end loop;
   exception
      when Ada.IO_Exceptions.Name_Error =>
         Text_IO.Put_Line (Name & " index not found?");
         GNAT.OS_Lib.OS_Exit (0);
   end Thread;

   --  =========================================================================
   --  Dipole_Model: Core climate modeling and parameter optimization procedure
   --
   --  PURPOSE:
   --    Fits a climate dipole model (e.g., ENSO, IOD, AMO) to observed data
   --    using tidal forcing + Laplace's Tidal Equation. Performs stochastic
   --    search through parameter space to maximize correlation coefficient
   --    (or minimize RMS, or other configurable metrics).
   --
   --  ALGORITHM:
   --    1. Load climate index data (e.g., NINO3.4 SST anomalies)
   --    2. Initialize parameters from shared state or defaults
   --    3. Main optimization loop:
   --       a. Generate tidal forcing from long-period constituents
   --       b. Apply LTE modulation (wave equation solution)
   --       c. Perform multivariate regression for modulation harmonics
   --       d. Calculate metric (CC, RMS, DTW, etc.) vs. observed data
   --       e. If improved, update shared state and notify Monitor
   --       f. Use random descent to adjust parameters and iterate
   --    4. Periodically save best parameters to disk
   --
   --  CONFIGURATION:
   --    Highly configurable via 50+ environment variables controlling:
   --    - Metric type (CC, RMS, DTW, EMD, etc.)
   --    - Training interval (dates, split, exclude)
   --    - Optimization parameters (spread, threshold, max loops)
   --    - LTE physics (filters, nonlinearity, decay, symmetry)
   --    - Advanced features (Mathieu mode, full-wave rectification, etc.)
   --
   --  CONCURRENCY:
   --    Multiple threads run this procedure concurrently, each with its own
   --    random search trajectory. Monitor coordinates best-result tracking.
   --  =========================================================================
   procedure Dipole_Model
     (N_Tides, N_Modulations : in Integer; ID : in Integer := 0;
      File_Name : in String := "nino4.dat";
      Split_Training : in Boolean := False; Exclude : in Boolean := False)
   is
      package LEF renames Ada.Numerics.Long_Elementary_Functions;

      -- Helper function to parse NH: supports "6" (count) or "1 1 1 1 1 1" (list)
      function Parse_NH (NH_Str : String) return Ns is
         Harms_Temp : constant Ns := S_to_I (NH_Str);
      begin
         -- If NH is a single integer > 1 (count), convert to array of 1's
         if Harms_Temp'Length = 1 and then Harms_Temp (1) > 1 then
            return (1 .. Harms_Temp (1) => 1);
         else
            -- Old format: space-separated list of integers
            return Harms_Temp;
         end if;
      end Parse_NH;

      Data_Records : Data_Pairs := Make_Data (File_Name);
      DR : Data_Pairs := Data_Records;

      function Impulse (Time : Long_Float) return Long_Float;

      function Impulse_Amplify is new Amplify (Impulse => Impulse);

      Impulses : Data_Pairs := Data_Records;
      Forcing : Data_Pairs := Data_Records;
      --  TODO: Can remove - F_Model was used for experimental blending of
      --  forcing calculation methods, never fully implemented. Line 837
      --  shows the intended 1% new / 99% old mixing logic that was abandoned.
      --F_Model  : Data_Pairs := Data_Records;
      Model : Data_Pairs := Data_Records;
      KeepModel : Data_Pairs := Data_Records;

      Best : Boolean := False;
      Percentage : Integer;
      Best_Client : Integer;

      D : Shared.Param_S (N_Tides, N_Modulations) :=
        GEM.LTE.Primitives.Shared.Get (N_Tides, N_Modulations);
      DKeep : Shared.Param_S (N_Tides, N_Modulations) := D;
      D0 : constant Shared.Param_S := D;  -- reference

      Maximum_Loops : constant Long_Integer :=
        GEM.Getenv ("MAXLOOPS", 100_000);
      Threshold : constant Integer := GEM.Getenv ("THRESHOLD", 99);
      Spread_Min : constant Long_Float :=
        GEM.Getenv ("SPREAD_MIN", 0.000_000_1);
      Spread_Max : constant Long_Float := GEM.Getenv ("SPREAD_MAX", 0.1);
      Spread_Cycle : constant Long_Float :=
        GEM.Getenv ("SPREAD_CYCLE", 1_000.0);
      Catchup : constant Boolean :=
        GEM.Getenv ("THRESHOLD_ACTION", "RESTART") = "CATCHUP";
      RMS_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "RMS";
      ZC_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "ZC";
      FT_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "FT";
      DTW_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "DTW";
      Winding_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "W";
      CID_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "CID";
      CTW_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "CTW";
      DER_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "DER";
      EMD_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "EMD";
      Hoy_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "HOYER";
      DTW_CC : constant Boolean := GEM.Getenv ("METRIC", "CC") = "DC";
      EMD_CC : constant Boolean := GEM.Getenv ("METRIC", "CC") = "EC";
      EMD_DER : constant Boolean := GEM.Getenv ("METRIC", "CC") = "ED";
      SEM_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "SEM";
      DTW_Window : constant Integer :=
        GEM.Getenv ("DTW", 1); -- 100 is about ~5% of the time-series length
      Sampling_Per_Year : constant Long_Float := GEM.Getenv ("SAMPLING", 12.0);
      Filter : constant Long_Float := GEM.Getenv ("FILTER", 0.333_333_33);
      MLR_On : constant Boolean := GEM.Getenv ("MLR", False); -- wrong name
      Forcing_Only : constant Boolean := GEM.Getenv ("FORCING", False);
      Filter9Pt : constant Integer := GEM.Getenv ("F9", 0);
      Climate_Trend : constant Boolean := GEM.Getenv ("TREND", False);
      NonLin : constant Long_Float := GEM.Getenv ("NONLIN", 1.0);
      Decay : constant Long_Float := GEM.Getenv ("DECAY", 1.0);
      Lock_Freq : constant Boolean := GEM.Getenv ("LOCKF", False);
      Local_Max : constant Boolean := GEM.Getenv ("LOCAL", False);
      Lock_Tidal : constant Boolean := GEM.Getenv ("LOCKT", False);
      Lock_T_Amp : constant Boolean := GEM.Getenv ("LOCKA", False);
      --  UNUSED CONSTANT (Compiler Warning):
      --  ImpMonth loaded from environment but never referenced in code.
      --
      --  RATIONALE: Was intended for hardcoded monthly impulse timing
      --  (see commented code lines 497-499 below). Testing showed the
      --  DPos calculation from DelB parameter is more flexible and general.
      --
      Vary_Initial : constant Boolean := GEM.Getenv ("VI", False);
      Initial_Conditions_Date : constant Long_Float :=
        GEM.Getenv ("IDATE", 0.0);
      RMS_Data : Long_Float := 0.0;

      function Metric (X, Y, Z : in Data_Pairs) return Long_Float is
      begin
         if RMS_Metric then
            return RMS (X, Y, RMS_Data, 0.0);
            --  TODO: Can remove - Experimental hybrid metric combining RMS and CC
            --  was tested but never adopted. Simple RMS alone was sufficient.
            -- return (RMS(X,Y,RMS_Data, 0.0) + CC(X,Y))/2.0;
         elsif CID_Metric then
            return CC (X, Y) * CID (X, Y);
         elsif ZC_Metric then
            return Xing (X, Y);
         elsif DER_Metric then
            return DER_CC (X, Y);
         elsif FT_Metric then
            return FT_CC (X, Y, Z);
         elsif Winding_Metric then
            return (Winding_Agreement (X, Y, Z) + CC(X,Y))/2.0;
         elsif DTW_Metric then
            return DTW_Distance (X, Y, DTW_Window);
         elsif EMD_Metric then
            return EMD (X, Y);
         elsif Hoy_Metric then
            return -- Ada.Numerics.Long_Elementary_Functions.Sqrt ()
                   (Hoyer_Spectral_Peak (X, Y, Z) + CC (X, Y)) * 0.5;
         elsif DTW_CC then
            return
              Ada.Numerics.Long_Elementary_Functions.Sqrt
                (DTW_Distance (X, Y, DTW_Window) * CC (X, Y));
         elsif CTW_Metric then
            return
              Ada.Numerics.Long_Elementary_Functions.Sqrt
                (Long_Float'Max (DTW_Distance (X, Y, DTW_Window), 0.0) *
                 CID (X, Y));
         elsif EMD_CC then
            return (EMD (X, Y) + CC (X, Y)) * 0.5;
         elsif EMD_DER then
            return (EMD (X, Y, Derivative => True) + CC (X, Y)) * 0.5;
         elsif SEM_Metric then
            return Scaled_Error_Metric (X, Y);  -- Y should be data
         elsif MLR_On then
            return CC (X, Y) * Min_Entropy_Power_Spectrum (Z, Y);
         elsif Is_Minimum_Entropy then
            return Min_Entropy_Power_Spectrum (X, Y); -- or X = Z
         else
            return CC (X, Y);
         end if;
      end Metric;

      function Impulse_Delta (Time : Long_Float) return Long_Float is
         Value : Long_Float;
         -- Impulses will occur on a month for monthly data
         Trunc : Integer :=
           Integer (((Time - Long_Float'Floor (Time)) * Sampling_Per_Year));
         DPos : Integer := Integer ((abs (D.B.DelB) * Sampling_Per_Year));
      begin
         if Trunc = DPos then
            Value := D.B.DelA;
         elsif Trunc = (DPos + Integer (Sampling_Per_Year) / 2) mod 12 then
            Value := D.B.Asym;
         else
            Value := 0.0;
         end if;
         return Value;
      end Impulse_Delta;

      function Impulse (Time : Long_Float) return Long_Float is
         Value : Long_Float;
      begin
         Value := Impulse_Delta (Time);
         return Value;
      end Impulse;

--      function Annual_Add
--        (Model : in Data_Pairs; Polarity : Long_Float := 1.0) return Data_Pairs
--      is
--         M : Data_Pairs := Model;
--         Pi : Long_Float := Ada.Numerics.Pi;
--         Annual_Factor, Semi_Factor : Long_Float := 1.0;
--         use Ada.Numerics.Long_Elementary_Functions;
--      begin
--         for I in Model'Range loop
--            M (I).Value :=
--              M (I).Value +
--              Polarity * Annual_Factor * D.B.Ann1 * Cos (2.0 * Pi * M (I).Date + D.B.Ann2) +
--              Polarity * Annual_Factor * Semi_Factor * D.B.Sem1 *
--                Cos (4.0 * Pi * M (I).Date + D.B.Sem2);
--         end loop;
--         return M;
--      end Annual_Add;

      function Bessel
        (Model : in Data_Pairs; 
         eS, eC, k, k2, eS2, eC2 : in Long_Float) return Data_Pairs
      is
         M : Data_Pairs := Model;
         Pi : Long_Float := Ada.Numerics.Pi;
         use Ada.Numerics.Long_Elementary_Functions;
      begin
         for I in Model'Range loop
            M (I).Value := M (I).Value + eS*Sin(2.0 * Pi * k * M (I).Value) + 
                                         eC*Cos(2.0 * Pi * k * M (I).Value) +
                                     eS2*eS*Sin(2.0 * Pi * k2* M (I).Value) + 
                                     eC2*eC*Cos(2.0 * Pi * k2* M (I).Value);
                               -- Secondary*eS*Sin(2.0 * Pi * k2* M (I).Value) + 
                               -- Secondary*eC*Cos(2.0 * Pi * k2* M (I).Value);
         end loop;
         return M;
      end Bessel;

      
      procedure Put_CC
        (Val1, Val2 : in Long_Float; Counter : in Long_Integer;
         Thread : in Integer)
      is
      begin
         Text_IO.Put (GEM.Getenv ("METRIC", "CC"));
         Ada.Long_Float_Text_IO.Put (Val1, Fore => 4, Aft => 10, Exp => 0);
         Ada.Long_Float_Text_IO.Put (Val2, Fore => 4, Aft => 10, Exp => 0);
         Text_IO.Put_Line ("  " & Thread'Img & Counter'Img);
      end Put_CC;

      der : Long_Float;
      CorrCoeff, Old_CC, Prior_Best_CC : Long_Float := 0.0;
      CorrCoeffP : Long_Float;
      CorrCoeffTest : Long_Float := 0.0;
      Progress_Cycle, Spread : Long_Float;
      Counter : Long_Integer := -1;

      function Find_Index (Time : in Long_Float) return Integer is
         Index : Integer;
      begin
         for I in Data_Records'Range loop
            exit when Data_Records (I).Date > Time; -- Finding indices at time
            Index := I;
         end loop;
         return Index;
      end Find_Index;

      TS : Long_Float :=
        GEM.Getenv ("TRAIN_START", Data_Records (Data_Records'First).Date);
      TE : Long_Float :=
        GEM.Getenv ("TRAIN_END", Data_Records (Data_Records'Last).Date);

      First : Integer := Find_Index (TS);
      Last : Integer := Find_Index (TE);
      Mid : Integer := (First + Last) / 2;

      function Exclude_Metric return Long_Float is
         X : Data_Pairs :=
           Model (Model'First .. First) & Model (Last .. Model'Last);
         Y : Data_Pairs :=
           Data_Records (Data_Records'First .. First) &
           Data_Records (Last .. Data_Records'Last);
         Z : Data_Pairs :=
           Forcing (Forcing'First .. First) & Forcing (Last .. Forcing'Last);
      begin
         return Metric (X, Y, Z);
      end Exclude_Metric;

      function Excluded (D : Data_Pairs) return Data_Pairs is
      begin
         if Exclude then
            return D (D'First .. First) & D (Last .. D'Last);
         else
            return D (First .. Last);
         end if;
      end Excluded;

      Max_Harmonics : Positive := GEM.Getenv ("MAXH", 1_000);

      ------------------------------------------------------------------------
   -- Array overlay for random descent optimization (Type-Safe Version)
   -- The Param_B record is overlaid with an array to allow Walker.Markov
   -- to randomly modify parameters. This is documented and verified via
   -- the Param_B_Overlay package for maximum safety.
      ------------------------------------------------------------------------
      NM : constant Integer := GEM.Getenv ("NM", N_Modulations);
      
      -- Full overlay size for verification (entire D.B structure)
      Full_Size : constant Positive :=
        GEM.LTE.Primitives.Param_B_Overlay.Overlay_Size (D.B.NLP, D.B.NLT);
      
      -- Reduced size for Walker search based on NM (only search used LT entries)
      -- Size = 18 scalars + (NLP * 2) LPAP + NM LT entries
      Size_Shared : constant Positive :=
        GEM.Getenv ("DSIZE", 18 + (D.B.NLP * 2) + NM);

      package Walker is new GEM.Random_Descent
        (Fixed => Is_Fixed, Set_Range => Size_Shared,
         Harmonic_Range => Max_Harmonics);

      Set, Keep, Set0 : Walker.LF_Array (1 .. Size_Shared);
      
      --  Suppress overlay address clause warning (intentional unsafe operation)
      pragma Warnings (Off, "overlay changes scalar storage order");
      for Set'Address use D.B.Offset'Address;
      pragma Warnings (On, "overlay changes scalar storage order");
      
      --  SAFETY NOTE: This address clause creates an array view of Param_B.
      --  The layout is verified by Param_B_Overlay.Verify_Layout below.
      --  Field positions are documented via named constants in the overlay package.
      --  Walker only searches 1..Size_Shared (NM-limited), not the full array.
      ------------------------------------------------------------------------
      Harms : Ns := Parse_NH (GEM.Getenv ("NH", ""));
      Harms_Keep : Ns := Harms;
      NH : Integer := Harms'Length;

      Secular_Trend : Long_Float := 0.0;
      --  TODO: Can remove - Single exclamation mark marker used during
      --  development for code bookmarking, no functional significance.
      --!
      Singular : Boolean;
      M : Modulations (1 .. NM + NH);
      MAP : Modulations_Amp_Phase (1 .. NM + NH);
      Accel : Long_Float;
      Annual_Cycle : Annual_Harmonics;
      Keep_Initial_Value, Init_Value0 : Long_Float := 0.0;

      function Calc_Forcing return Data_Pairs is
         F : Data_Pairs := Forcing;
      begin
         der := 1.0 - D.B.mA; -- keeps the integrator stable

         Impulses :=
           Impulse_Amplify
             (Raw =>
                Tide_Sum
                  (Template => Data_Records, Constituents => D.B.LPAP,
                   Periods => D.A.LP, Ref_Time => 0.0, Scaling => 0.0,
                   Year_Len => Year_Length, Integ => D.B.ShiftT),
              Offset => 0.0, Ramp => 0.0,
              Start => Data_Records (Data_Records'First).Date);

            F :=
              IIR
                (Raw => Impulses, lagA => der, -- lagB => D.B.mA,
                 lagC => D.B.mP, iA => D.B.init, -- iB => 0.0, iC => 0.0, 
                 Start => Initial_Conditions_Date);

         return F;
      end Calc_Forcing;

   begin
      if Data_Records'Length = 0 then
         Text_IO.Put_Line (File_Name & " empty or not found");
         GNAT.OS_Lib.OS_Exit (0);
      end if;

      for I in First .. Last loop
         RMS_Data :=
           RMS_Data + Data_Records (I).Value * Data_Records (I).Value;
      end loop;
      RMS_Data := Ada.Numerics.Long_Elementary_Functions.Sqrt (RMS_Data);
      Old_CC := 0.0;
      
      -- Verify overlay layout before optimization begins (use full size)
      GEM.LTE.Primitives.Param_B_Overlay.Verify_Layout (D.B, Full_Size);
      
      -- Optional debug: print field-by-field mapping (set OVERLAY_DEBUG=1)
      declare
         --  Create view of Walker.LF_Array as Param_B_Overlay.LF_Array for debug
         Debug_Set : GEM.LTE.Primitives.Param_B_Overlay.LF_Array (1 .. Size_Shared);
         for Debug_Set'Address use Set'Address;
      begin
         GEM.LTE.Primitives.Param_B_Overlay.Debug_Print_Field_Mapping (D.B, Debug_Set);
      end;
      
      Walker.Reset;
      if Filter9Pt > 0 then
         for F in 1 .. Filter9Pt loop
            Data_Records := Filter9Point (Data_Records);
         end loop;
      end if;
      Text_IO.Put_Line ("Catchup mode enabled:" & Boolean'Image (Catchup));
      Keep := Set;
      Set0 := Set;
      Init_Value0 := D.B.init;

      for I in 1 .. Harms'Length loop
         exit when D.C (I) = 0;
         Harms (I) := D.C (I);
      end loop;

      loop
         Counter := Counter + 1;
         delay 0.0; -- context switching point if multi-processing not avilable

         -- Tidal constituents summed, amplified by impulse, and LTE modulated
         Forcing := Calc_Forcing;

         M (1 .. NM) := D.B.LT (1 .. NM);
         if Lock_Freq then
            M (NM) := 1.0 / (Decay * D.B.mP);
         end if;
         MAP (1 .. NM) := D.A.LTAP (1 .. NM);
         for I in 1 .. NH loop
            M (NM + I) := Long_Float (Harms (I)) * M (NM);
         end loop;
         
         -- Check for duplicate periods that could cause singular regression matrix
         declare
            Has_Duplicates : Boolean := False;
         begin
            -- Check if harmonic periods match any base periods
            for I in 1 .. NH loop
               for J in 1 .. NM loop
                  if abs (M (NM + I) - M (J)) < 1.0E-10 then
                     Text_IO.Put_Line ("Error: Harmonic period matches base period - cannot perform regression");
                     Has_Duplicates := True;
                     exit;
                  end if;
               end loop;
               exit when Has_Duplicates;
            end loop;
            
            -- Check for duplicate harmonic periods
            if not Has_Duplicates and NH > 1 then
               for I in 1 .. NH - 1 loop
                  for J in I + 1 .. NH loop
                     if abs (M (NM + I) - M (NM + J)) < 1.0E-10 then
                        Text_IO.Put_Line ("Error: Duplicate harmonic periods detected - cannot perform regression");
                        Has_Duplicates := True;
                        exit;
                     end if;
                  end loop;
                  exit when Has_Duplicates;
               end loop;
            end if;
            
            if Has_Duplicates then
               Singular := True;
               return;
            end if;
         end;
         
         if MLR_On or not (Forcing_Only or Is_Minimum_Entropy) then
            if Climate_Trend then
               Secular_Trend := 1.0;
            else
               Secular_Trend := 0.0;
            end if;
            DR := Data_Records;
--            DR := Annual_Add (DR, -1.0);
            if NM = 1 then
               Forcing := Bessel(Forcing, D.B.ImpA, D.B.ImpB, M(NM), 0.0, 0.0, 0.0);
            else
               if Lock_Freq then
                  Forcing := Bessel(Forcing, D.B.ImpA, D.B.ImpB, M(NM-1), M(NM), D.B.ImpC, D.B.IR);
               else
                  Forcing := Bessel(Forcing, D.B.ImpA, D.B.ImpB, M(NM), M(NM-1), D.B.ImpC, D.B.IR);
               end if;
            end if;
            Regression_Factors
              (Data_Records => Excluded (DR), -- Time series
               Forcing => Excluded (Forcing),  -- Value @ Time
               NM => NM + NH, -- # modulations
               DBLT => M, --D.B.LT,
               DALTAP => MAP, --D.A.LTAP,
               DALEVEL => D.A.level,
               DAK0 => D.A.k0, 
               Secular_Trend => Secular_Trend, 
               Accel => Accel,
               Singular => Singular, 
               Annual => Annual_Cycle,
               Third => 0.0 
            );
         else
            Singular := False;
         end if;

         CorrCoeff := 0.0;
            if Forcing_Only or (Is_Minimum_Entropy and not MLR_On) then
               Model := Forcing;
            else
               Model :=
                 LTE
                   (Forcing => Forcing,
                    Wave_Numbers =>
                      M (1 .. NM+NH), --D.B.Lt(1..NM),
                    Amp_Phase => MAP (1 .. NM+NH), --D.A.LTAP,
                    Offset => D.A.level, 
                    K0 => D.A.k0, 
                    Trend => Secular_Trend,
                    Accel => Accel, 
                    NonLin => NonLin,
                    Annual => Annual_Cycle,
                    Third => 0.0);  
               if Monotonic_Increase then
                  Secular_Trend := abs Secular_Trend;
                  Accel := abs Accel;
               end if;

--               Model := Annual_Add (Model);
            end if;


            if Filter9Pt > 0 then
               for F in 1 .. Filter9Pt loop
                  Model := Filter9Point (Model);
               end loop;
            else
               -- extra filtering, 2 equal-weighted 3-point box windows creating triangle
               Model := Median (Model);
               Model :=
                 FIR
                   (FIR (Model, Filter, 1.0 - 2.0 * Filter, Filter), Filter,
                    1.0 - 2.0 * Filter, Filter);
            end if;

            -- pragma Debug ( Dump(Model, Data_Records, Run_Time) );

            if Split_Training then
               if Split_Low then
                  CorrCoeff :=
                    Metric
                      (Model (First .. Mid), Data_Records (First .. Mid),
                       Forcing (First .. Mid));
                  CorrCoeffTest :=
                    Metric
                      (Model (Mid .. Last), Data_Records (Mid .. Last),
                       Forcing (Mid .. Last));
               else
                  CorrCoeffTest :=
                    Metric
                      (Model (First .. Mid), Data_Records (First .. Mid),
                       Forcing (First .. Mid));
                  CorrCoeff :=
                    Metric
                      (Model (Mid .. Last), Data_Records (Mid .. Last),
                       Forcing (Mid .. Last));
               end if;
            else
               if Exclude then
                  CorrCoeffP :=
                    Exclude_Metric; -- ( Model(First..Last), Data_Records(First..Last), Forcing(First..Last));
               else
                  CorrCoeffP :=
                    Metric
                      (Model (First .. Last), Data_Records (First .. Last),
                       Forcing (First .. Last));
               end if;
            end if;

            CorrCoeff := CorrCoeffP;

         if not Split_Training then
            if Exclude then  -- calculate OOB
               CorrCoeffP :=
                 Metric
                   (Model (First .. Last), Data_Records (First .. Last),
                    Forcing (First .. Last));
            else
               CorrCoeffP :=
                 Exclude_Metric; -- ( Model(First..Last), Data_Records(First..Last), Forcing(First..Last));
            end if;
         end if;

         -- Register the results with a monitor
         if Split_Training then
            Monitor.Check
              (CorrCoeffTest, CorrCoeff, ID, Counter, Best, Best_Client,
               Percentage);
         else
            Monitor.Check
              (CorrCoeff, CorrCoeffP, ID, Counter, Best, Best_Client,
               Percentage);
         end if;

         if ID = Best_Client then
            Counter := 1; -- no use penalizing thread in the lead
         end if;

         if CorrCoeff > Old_CC then
            Old_CC := CorrCoeff;
            Keep := Set;
            KeepModel := Model;
            DKeep := D;
            if Catchup then  -- save it for other threads to reset from
               GEM.LTE.Primitives.Shared.Put (D);
            end if;
            Harms_Keep := Harms;
         elsif Local_Max and CorrCoeff > Prior_Best_CC then
            -- Don't revert to Keep values
            Prior_Best_CC := CorrCoeff;
         else
            -- Go back to starting point (Keep) if local max not retained
            if Vary_Initial then
               D.B.init := Keep_Initial_Value;
            else
               Set := Keep;
               Harms := Harms_Keep;
            end if;
         end if;

         if Singular or
           ((not Best) and Counter > Maximum_Loops and Percentage < Threshold)
         then
            Text_IO.Put_Line ("Resetting" & ID'Img & Percentage'Img & "%");
            if Catchup then
               D := GEM.LTE.Primitives.Shared.Get (N_Tides, N_Modulations);
            else -- Restart
               D := D0; -- load back reference model parameters
            end if;
            Counter := 1;
            Old_CC := 0.0;
         end if;

         exit when Halted;

         Progress_Cycle := Long_Float (Counter);
         -- This slowly oscillates to change the size of the step to hopefully
         -- help it escape local minima, every N cycles
         if Counter = 0 then
            Spread := 0.0;
         else
            Spread :=
              Spread_Min +
              Spread_Max * (1.0 - LEF.Cos (Progress_Cycle / Spread_Cycle));
         end if;
         if Test_Only then
            exit;
         elsif Vary_Initial then
            Walker.Markov (D.B.init, Keep_Initial_Value, Spread, Init_Value0);
         else
            if Lock_Tidal then
               declare
                  DBLAP : constant Amp_Phases := D.B.LPAP;
               begin
                  Walker.Markov (Set, Keep, Spread, Set0);
                  D.B.LPAP := DBLAP;
               end;
            elsif Lock_T_Amp then
               declare
                  DBLAP : constant Amp_Phases := D.B.LPAP;
               begin
                  Walker.Markov (Set, Keep, Spread, Set0);
                  for I in DBLAP'Range loop
                     D.B.LPAP (I).Amplitude := DBLAP (I).Amplitude;
                  end loop;
               end;
            else
               Walker.Markov (Set, Keep, Spread, Set0);
            end if;
            Walker.Random_Harmonic (Harms, Harms_Keep);
         end if;

      end loop;
      Monitor.Stop;
      if Test_Only or Best_Client = ID then

         -- Text_IO.Put_Line("### " & File_Name);
         -- Walker.Dump(Keep); -- Print results of last best evaluation,
         GEM.LTE.Primitives.Shared.Dump (DKeep);

         Text_IO.Put_Line ("---- LTE ----");
         Put (Secular_Trend, " :trend:", NL);
         Put (Accel, " :accel:", NL);
         Put (D.A.k0, " :K0:", NL);
         Put (D.A.level, " :level:", NL);
         for I in 1 .. NM loop
            Put (M (I), ", ");  --D.B.LT
            Put (MAP (I).Amplitude, ", "); --D.A.LTAP
            Put (MAP (I).Phase, Integer (M (I) / M (NM))'Img, NL);
         end loop;
         for I in NM + 1 .. NM + NH loop
            Put (M (I), ", ");  --D.B.LT
            Put (MAP (I).Amplitude, ", "); --D.A.LTAP
            Put (MAP (I).Phase, Integer (M (I) / M (NM))'Img, NL);
            if I - NM <= DKeep.C'Last then
               DKeep.C (I - NM) := Integer (M (I) / M (NM));
            end if;
         end loop;
         Text_IO.Put_Line ("```");

         GEM.LTE.Primitives.Shared.Save (DKeep);
         Save (KeepModel, Data_Records, Forcing);    -- saves to file
         if Split_Training then
            Put_CC (CorrCoeff, CorrCoeffTest, Counter, ID);
         else
            CorrCoeffP := Exclude_Metric;
            CorrCoeff :=
              Metric
                (Model (First .. Last), Data_Records (First .. Last),
                 Forcing (First .. Last));
            Put_CC (CorrCoeff, CorrCoeffP, Counter, ID);
         end if;

         CorrCoeff := CompareRef (D.A.LP, LPRef, D.B.LPAP);
         Put (CorrCoeff, ":dLOD:   ");
         Put (Year_Length, ":Yr: " & File_Name);

      else
         null; -- Text_IO.Put_Line("Exited " & ID'Img);
      end if;

   exception
      when E : others =>
         Text_IO.Put_Line
           ("Solution err: " & Ada.Exceptions.Exception_Information (E));
         -- The following may need a debug-specifi compiler switch to activate
         Text_IO.Put_Line (GNAT.Traceback.Symbolic.Symbolic_Traceback (E));

   end Dipole_Model;

   function Check_Every_N_Loops return Integer is
   begin
      if Is_Split then
         return 1;
      else
         return NLoops;
      end if;
   end Check_Every_N_Loops;

   procedure Set_Trigger (Level : Integer) is
   begin
      Trigger := Long_Float (Level) / 10.0;
   end Set_Trigger;

end GEM.LTE.Primitives.Solution;
