with Ada.Strings.Fixed;
with GEM.Zonal;
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
with GEM.dLOD;
with GNAT.Ctrl_C;
with System.Task_Info;
with System.Multiprocessors;
with GEM.LTE.Primitives.Shared;
with GEM.LTE.Primitives.Param_B_Overlay;
with Ada.Exceptions;
with GNAT.Traceback.Symbolic;
with Ada.IO_Exceptions;
with GNAT.OS_Lib;
with Ada.Directories;

package body GEM.LTE.Primitives.Solution is

   Is_Split : constant Boolean := GEM.Getenv ("SPLIT_TRAINING", False);
   Split_Low : constant Boolean := GEM.Getenv ("SPLIT_LOW", True);
   Alternate : constant Boolean := GEM.Getenv ("ALTERNATE", False);
   --  SECULAR (formerly Monotonic_Increase) removed 2026-09-29: it took
   --  abs() of Secular_Trend/Accel AFTER LTE had already built the model
   --  from the signed MLR values, so it never constrained the model --
   --  it only flipped the sign of what was printed and saved to
   --  lt.exe.windings.json whenever the MLR found a negative trend or
   --  accel (seen on kS020_W010: saved accel +5.6e-6, model used -5.6e-6).
   Trigger : Long_Float := GEM.Getenv ("TRIGGER", 0.999_999);
   Ratio : constant Long_Float := GEM.Getenv ("RATIO", 0.0);
   NLoops : constant Integer := GEM.Getenv ("NLOOPS", 20); --100
   Test_Only : constant Boolean := GEM.Getenv ("TEST_ONLY", False);

   --  LEGS = 2: the worker tasks are divided into two leg groups. Tasks come
   --  in pairs, a leader on leg 1 and a follower on leg 2. Each leg takes its
   --  own settings (index file, SAMPLING, TRAIN_START, TRAIN_END, EXCLUDE,
   --  ... anything a response file can hold) from "<exe>.resp.1" and
   --  "<exe>.resp.2" ahead of the common response file. The leader walks
   --  the parameters; every candidate is evaluated by both tasks of the pair
   --  at the same time, each on its own leg, and the pair's Exchange (below)
   --  combines the two scores into the one composite the search accepts or
   --  rejects on. One parameter file is saved; results and windings are
   --  saved per leg (lte_results.1.csv, lte_results.2.csv, ...).
   --  Default 0: no legs, everything as before.
   --  LEGS = N (2 .. 9) generalizes the pair to a group of N tasks: one
   --  leader on leg 1 and a follower on each of legs 2 .. N, with settings
   --  from "<exe>.resp.1" .. "<exe>.resp.N". The composite is the
   --  LEG_WEIGHT-weighted mean over all N legs.
   Max_Group : constant := 9;
   Legs_N : constant Integer :=
     Integer'Min (Max_Group, Integer'(GEM.Getenv ("LEGS", 0)));
   Paired : constant Boolean := Legs_N >= 2;
   Group : constant Positive := Integer'Max (1, Legs_N);
   subtype Follower_Leg is Integer range 2 .. Max_Group;
   Active_Threads : Natural := 0;  -- set by Start

   --  SELECT_TOP = N (default 0 = off, the old rule: the search with the best
   --  training score saves). With N > 0 the choice is made once, when every
   --  search has stopped: of the (at most) N best training scores that lie
   --  within SELECT_TIE of the best one -- a statistical tie on training --
   --  the one with the highest held-out score is saved. Works with or
   --  without legs (with legs the scores are the pair's composites).
   Select_Top : constant Integer := GEM.Getenv ("SELECT_TOP", 0);
   Select_Tie : constant Long_Float := GEM.Getenv ("SELECT_TIE", 0.01);
   Max_Clients : constant := 512;
   type Score_List is array (1 .. Max_Clients) of Long_Float;

   function Assign_Leg (L : Natural) return Boolean is
   begin
      GEM.Set_Leg (L);
      return True;
   end Assign_Leg;

   --  Compare model output against reference tidal data (dLOD - day length)
   --  Used for validation rather than optimization - generates CSV files
   --  showing how well the tidal reconstruction matches reference data.
   function CompareRef
     (LP : in Long_Periods; AP : in Long_Periods_Amp_Phase;
      Year_Correction : in Long_Float; Rho : in Long_Float := 1.0;
      ZTau : in Long_Float := 1.0; ZDist : in Long_Float := 1.0;
      ZTilt : in Long_Float := 0.0; SsaA : in Long_Float := 0.0;
      SsaP : in Long_Float := 0.0)
      return Long_Float
   is
      dLOD_dat : String := GEM.Getenv ("DLOD_DAT", "../dlod3.dat");
      D : Data_Pairs := Make_Data (dLOD_dat);
      Ref, R2, M2 : Data_Pairs := D;
      Metric : Long_Float := -1.0;
      Ref_Time : Long_Float := 0.0;
      LPRef : constant Long_Periods_Amp_Phase :=
        GEM.dLOD (dLOD_dat, Year_Correction);
      File : Text_IO.File_Type;
   begin
      R2 :=
        Tide_Sum
          (Template => Ref, Constituents => LPRef, Periods => LP,
           Ref_Time => Ref_Time, Scaling => 1.0, Cos_Phase => False,
           Year_Len => Year_Length (Year_Correction));

      M2 :=
        Tide_Sum
          (Template => Ref, Constituents => AP, Periods => LP,
           Ref_Time => Ref_Time, Scaling => 1.0, Cos_Phase => False,
           Year_Len => Year_Length (Year_Correction), Rho => Rho,
           ZTau => ZTau, ZDist => ZDist, ZTilt => ZTilt, SsaA => SsaA,
           SsaP => SsaP);
      --  TIDES=ZONAL: Tide_Sum returns the generator for both calls, so
      --  the gate compares the generator directly with the observed dLOD.
      if GEM.Zonal.On then
         Metric := CC (D, M2);
      else
         Metric := Long_Float'Max (CC (R2, M2), Metric);
      end if;
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
     (CPU : System.Multiprocessors.CPU; N_Tides, N_Modulations : Integer;
      Leg, Pair, Index : Natural) is
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
      --  Paired legs need an even number of tasks (at least one pair)
      NT : constant Positive :=
        (if Paired
         then Integer'Max
                (Group, Number_of_Threads - Number_of_Threads mod Group)
         else Number_of_Threads);
      TA : array (1 .. NT) of Thread_Access;
      use System.Multiprocessors;
      Num : constant Positive := System.Task_Info.Number_Of_Processors;
      CPU : Positive;
   begin
      GNAT.Ctrl_C.Install_Handler (Handler => GEM.LTE.Primitives.Stop'Access);
      Active_Threads := NT;
      for I in 1 .. NT loop
         CPU := I mod Num + 1;
         TA (I) :=
           new Thread
             (CPU => CPU_Range (CPU), N_Tides => N_Tides,
              N_Modulations => N_Modulations,
              Leg => (if Paired then (I - 1) mod Group + 1 else 0),
              Pair => (if Paired then (I - 1) / Group + 1 else 0),
              Index => I);
         TA (I).Start;
         delay 0.1; -- let them gradually start up
      end loop;
   end Start;

   Worst_Case : constant Long_Float := GEM.Getenv ("STARTING_METRIC", 0.001);

   --  Same value enso_opt.adb used to size the thread pool via Start's own
   --  Number_of_Threads argument -- re-read here (rather than threading a
   --  new parameter through Start/Thread/Dipole_Model) since GEM.Getenv's
   --  resp/env state is fixed for the process's lifetime, so this is
   --  guaranteed to match. Used only by the VALIDATE lockbox below, to
   --  know how many final reports to wait for.
   Validate_Thread_Count : constant Positive :=
     GEM.Getenv ("NUMBER_OF_PROCESSORS", System.Task_Info.Number_Of_Processors);

   --  Final-candidate lockbox (see Monitor.Report_Final/Winner below):
   --  every thread's own search, accept/reject and Catchup resets are
   --  completely unaffected by this -- it only changes which thread's
   --  finished candidate actually gets saved at the very end, judged by
   --  performance on the held-out [Last, D'Last] segment rather than by
   --  whichever thread had the best TRAINING-window fit. Deliberately a
   --  one-time, end-of-run decision, not a live search signal (see the
   --  TNSTAAFL discussion this was designed against). Default False:
   --  identical to today's Best_Client-based save. Package-level (rather
   --  than local to Dipole_Model) so the live-progress Status function
   --  below can also see it, to decide whether to show a validate column.
   Validate : constant Boolean := GEM.Getenv ("VALIDATE", False);

   --  Access types for the VALIDATE lockbox (see Monitor.Report_Final
   --  below): Param_S is discriminated and Data_Pairs/Modulations/
   --  Modulations_Amp_Phase are unconstrained, none of which a protected
   --  object can hold directly -- same reason GEM.LTE.Primitives.Shared's
   --  own Server uses an access type for Param_S.
   type Final_Param_P is access all GEM.LTE.Primitives.Shared.Param_S;
   type Final_Model_P is access all Data_Pairs;
   type Final_M_P is access all Modulations;
   type Final_MAP_P is access all Modulations_Amp_Phase;

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
   --  Serializes the final save block. In TEST_ONLY every thread saves, and
   --  two threads writing lt.exe.p / lte_results.csv / dlod_ref.dat at once
   --  (while the other reads dlod_ref.dat back) aborted mid-save with
   --  END_ERROR, leaving a 0-byte lte_results.csv (26 of 89 re-scores,
   --  2026-10-01). Release is a no-op unless the caller is the owner, so the
   --  exception handler can call it unconditionally.
   protected Save_Mutex is
      entry Seize (Who : in Integer);
      procedure Release (Who : in Integer);
   private
      Owner : Integer := -1;
   end Save_Mutex;

   protected body Save_Mutex is
      entry Seize (Who : in Integer) when Owner = -1 is
      begin
         Owner := Who;
      end Seize;
      procedure Release (Who : in Integer) is
      begin
         if Owner = Who then
            Owner := -1;
         end if;
      end Release;
   end Save_Mutex;

   protected Monitor is
      procedure Check
        (Metric : in Long_Float;
         OOB : in Long_Float; -- Out-of-band
         Client : in Integer;
         Count : in Long_Integer;
         Validate_Live : in Long_Float; -- informational only, see below
         Best : out Boolean;  -- accessing thread deemed best
         BestClient : out Integer;
         Percentage : out Integer); -- % of best metric if not best
      procedure Client
        (ID : out Integer); -- registering a thread ID, called once
      entry Status
        (Metric : out Long_Float;   -- used by a monitoring thread, i.e. main
         OOB : out Long_Float;
         Client : out Integer;      -- returns client thread w/ best metric
         Cycle : out Long_Integer;  -- and the cycle count it is on
         Validate_Live : out Long_Float); -- see Check's own comment
      procedure Stop;
      --  TEST_ONLY: each worker calls this once when it has finished (scored,
      --  saved, printed). When all have, the program is halted and the
      --  main thread's Status call released, so a test run ends by itself
      --  instead of waiting to be killed.
      procedure Test_Done;
      procedure Reset;
      --  The task currently holding the best score (for a follower, which
      --  never calls Check itself, to learn whether its pair leads).
      procedure Leader (Client : out Integer);
      --  SELECT_TOP: every search reports its kept candidate's training and
      --  held-out scores once, when it stops; Selected then tells each
      --  task (a follower asks for its leader) whether its candidate is
      --  the one to save. It blocks until all searches have reported.
      procedure Final_Report
        (Client : in Integer; Train, Held_Out : in Long_Float);
      entry Selected (Client : in Integer; Yes : out Boolean);

      --  VALIDATE lockbox (see Dipole_Model's own use of it): each thread
      --  calls this exactly once, when it stops searching, handing over
      --  its own best-found candidate and that candidate's score on the
      --  held-out validation segment. Validation NEVER feeds Check/
      --  Best_Metric above -- it only ever decides which already-finished
      --  candidate gets saved, never steers the search itself (see the
      --  TNSTAAFL discussion this was designed against: blending it into
      --  the live optimization target would just make it a second,
      --  differently-weighted training signal, silently burning the one
      --  thing a held-out check is for). Am_I_Last tells the calling
      --  thread whether it is the one responsible for performing the
      --  actual file save, once every thread has reported.
      procedure Report_Final
        (D : in GEM.LTE.Primitives.Shared.Param_S;
         Model : in Data_Pairs;
         Forcing : in Data_Pairs;
         M : in Modulations;
         MAP : in Modulations_Amp_Phase;
         Trend, Accel, Aero : in Long_Float;
         Validate_Score : in Long_Float;
         Am_I_Last : out Boolean);
      procedure Winner
        (D : out GEM.LTE.Primitives.Shared.Param_S;
         Model : out Data_Pairs;
         Forcing : out Data_Pairs;
         M : out Modulations;
         MAP : out Modulations_Amp_Phase;
         Trend, Accel, Aero : out Long_Float);

      --  Non-blocking peek at the SAME three numbers the live Status
      --  display shows (Best_Metric/Best_OOB/Best_Validate_Live) --
      --  unlike the Status entry above, this returns immediately rather
      --  than waiting for the next improvement. Used by the VALIDATE=TRUE
      --  exit report so what gets printed at exit is exactly a copy of
      --  the last live display, for the leader, rather than a value
      --  freshly (and inconsistently) recomputed from Model/KeepModel.
      procedure Current
        (Metric : out Long_Float; OOB : out Long_Float;
         Validate_Live : out Long_Float);
   private
      --  Current best metric and associated metadata
      Best_Metric : Long_Float :=
        Worst_Case; -- so doesn't cause overflow for %
      Best_OOB : Long_Float := Worst_Case;
      Client_Index : Integer := 1;
      Waiting : Boolean := True; -- triggers status only if value changes
      Tests_Done : Natural := 0;
      Fin_Train : Score_List := (others => Long_Float'First);
      Fin_Held : Score_List := (others => Long_Float'First);
      Reports : Natural := 0;
      Chosen : Integer := 0;
      Best_Client : Integer := -1;
      Best_Count : Long_Integer := 0;

      --  Live, informational-only running MAXIMUM of every thread's own
      --  validate score (each thread's own KeepModel, scored per
      --  Validate_Metric) seen so far -- purely for the progress display
      --  (Status below), so a run can be watched to see the best
      --  generalization found by ANY thread yet, not just whichever
      --  thread currently happens to be leading on training. Tracked
      --  unconditionally in Check below (independent of Best_Metric's
      --  own update), never feeding Check's own accept/reject decision.
      Best_Validate_Live : Long_Float := 0.0;

      --  VALIDATE lockbox state -- deliberately separate from Best_Metric/
      --  Best_Client above; the two never mix.
      Threads_Reported : Natural := 0;
      Best_Validate_Score : Long_Float := Long_Float'First;
      Best_Validate_D : Final_Param_P := null;
      Best_Validate_Model : Final_Model_P := null;
      --  Forcing was missing from this lockbox entirely (found 2026-09-29
      --  while cross-checking experiments/Feb2026/lte_forward.py's
      --  reimplementation against the live Ada pipeline): Winner swaps
      --  D/Model/M/MAP with whichever THREAD actually had the best
      --  Validate_Score, which need not be the reporting thread itself --
      --  but the reporting thread's own local KeepForcing was never part
      --  of that exchange, so Forcing := KeepForcing (added earlier the
      --  same day, see Dipole_Model) still only resynced within a single
      --  thread's own accept history, not across the cross-thread swap
      --  Winner performs. Confirmed directly: a fresh, single-threaded,
      --  zero-perturbation replay of a saved lt.exe.p reproduced the
      --  live Ada Calc_Forcing/Bessel output exactly (to 1e-13), yet
      --  still disagreed with that same run's own saved lte_results.csv
      --  Forcing column by up to ~80% relative on cells where the
      --  reporting thread and the true Validate winner differed.
      Best_Validate_Forcing : Final_Model_P := null;
      Best_Validate_M : Final_M_P := null;
      Best_Validate_MAP : Final_MAP_P := null;
      Best_Validate_Trend : Long_Float := 0.0;
      Best_Validate_Accel : Long_Float := 0.0;
      Best_Validate_Aero : Long_Float := 0.0;
   end Monitor;

   protected body Monitor is
      --  Check if this thread's metric is the new best
      --  Implements trigger-based early termination if metric exceeds threshold
      procedure Check
        (Metric : in Long_Float;
         OOB : in Long_Float; -- Out-of-band
         Client : in Integer;
         Count : in Long_Integer;
         Validate_Live : in Long_Float; Best : out Boolean;
         BestClient : out Integer;
         Percentage : out Integer)
      is
         M : Long_Float := Metric + Ratio * OOB;
      begin
         --  BUG FIX: this used to live inside the "if M >= Best_Metric"
         --  branch below, so the displayed "validate" number was really
         --  "whatever the current TRAINING leader's own validate score
         --  happens to be" -- since that leader's Validate_Live and OOB
         --  are both reported from the SAME triggering call, they tended
         --  to move in lockstep, making validate look like it was just
         --  mirroring test. Tracking it here, unconditionally, makes it
         --  a genuine running maximum across EVERY thread's own report --
         --  including threads that never become the training leader --
         --  matching what Report_Final/Winner actually select at the
         --  end, rather than an incidental byproduct of train's own
         --  record-setting. Does not change Waiting/the display's own
         --  refresh cadence (still gated by Best_Metric improving), only
         --  the VALUE shown when it does refresh.
         if Validate_Live > Best_Validate_Live then
            Best_Validate_Live := Validate_Live;
         end if;
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
         Cycle : out Long_Integer;
         Validate_Live : out Long_Float)
        when not Waiting is
      begin
         Waiting := True;
         Metric := Best_Metric;
         OOB := Best_OOB;
         Client := Best_Client;
         Cycle := Best_Count;
         Validate_Live := Best_Validate_Live;
      end Status;

      --  Force Status entry to release (for clean shutdown)
      procedure Stop is
      begin
         Waiting :=
           False; -- necesary to allow a clean exit when program halted
      end Stop;

      procedure Test_Done is
      begin
         Tests_Done := Tests_Done + 1;
         if Tests_Done >=
           (if Active_Threads > 0 then Active_Threads
            else Validate_Thread_Count)
         then
            GEM.LTE.Primitives.Stop;
            Waiting := False;
         end if;
      end Test_Done;

      procedure Leader (Client : out Integer) is
      begin
         Client := Best_Client;
      end Leader;

      procedure Final_Report
        (Client : in Integer; Train, Held_Out : in Long_Float) is
      begin
         if Client in Score_List'Range then
            Fin_Train (Client) := Train;
            Fin_Held (Client) := Held_Out;
         end if;
         Reports := Reports + 1;
      end Final_Report;

      entry Selected (Client : in Integer; Yes : out Boolean)
        when Reports >=
          (if Paired then Active_Threads / Group else Active_Threads)
      is
      begin
         if Chosen = 0 then
            declare
               Top : Long_Float := Long_Float'First;
               Used : array (Score_List'Range) of Boolean := (others => False);
               Pick : Integer;
            begin
               for I in Score_List'Range loop
                  Top := Long_Float'Max (Top, Fin_Train (I));
               end loop;
               Text_IO.Put_Line
                 ("SELECT_TOP" & Integer'Max (1, Select_Top)'Img &
                  ": best training score" &
                  Top'Img & ", tie width" & Select_Tie'Img);
               for N in 1 .. Integer'Max (1, Select_Top) loop
                  Pick := 0;
                  for I in Score_List'Range loop
                     if not Used (I) and then Fin_Train (I) > Long_Float'First
                       and then Fin_Train (I) >= Top - Select_Tie
                       and then (Pick = 0
                                 or else Fin_Train (I) > Fin_Train (Pick))
                     then
                        Pick := I;
                     end if;
                  end loop;
                  exit when Pick = 0;
                  Used (Pick) := True;
                  Text_IO.Put_Line
                    ("   candidate" & Pick'Img & ": training" &
                     Fin_Train (Pick)'Img & ", held out" &
                     Fin_Held (Pick)'Img);
                  if Chosen = 0 or else Fin_Held (Pick) > Fin_Held (Chosen)
                  then
                     Chosen := Pick;
                  end if;
               end loop;
               if Chosen = 0 then
                  Chosen := -1;  -- nothing to save
               end if;
               Text_IO.Put_Line ("   selected" & Chosen'Img);
            end;
         end if;
         Yes := Client = Chosen;
      end Selected;

      --  Reset best metric tracking (used in alternating exclude mode)
      procedure Reset is
      begin
         Best_Metric := 0.0;
         Best_OOB := 0.0;
         Best_Validate_Live := 0.0;
         --  Also reset the VALIDATE lockbox: without this, a second
         --  Alternate round would inherit the previous round's already-
         --  satisfied Threads_Reported count (so Am_I_Last would fire on
         --  the very first report) and its stale winner.
         Threads_Reported := 0;
         Best_Validate_Score := Long_Float'First;
         Best_Validate_D := null;
         Best_Validate_Model := null;
         Best_Validate_Forcing := null;
         Best_Validate_M := null;
         Best_Validate_MAP := null;
         Best_Validate_Trend := 0.0;
         Best_Validate_Accel := 0.0;
         Best_Validate_Aero := 0.0;
      end Reset;

      procedure Report_Final
        (D : in GEM.LTE.Primitives.Shared.Param_S;
         Model : in Data_Pairs;
         Forcing : in Data_Pairs;
         M : in Modulations;
         MAP : in Modulations_Amp_Phase;
         Trend, Accel, Aero : in Long_Float;
         Validate_Score : in Long_Float;
         Am_I_Last : out Boolean)
      is
      begin
         --  Always reallocate rather than reuse-in-place: Model/M/MAP are
         --  unconstrained, and different reports are only guaranteed to
         --  match in length within a single run, not necessarily safe to
         --  assign into a previous allocation's storage.
         if Best_Validate_D = null or else Validate_Score > Best_Validate_Score
         then
            Best_Validate_Score := Validate_Score;
            Best_Validate_D := new GEM.LTE.Primitives.Shared.Param_S'(D);
            Best_Validate_Model := new Data_Pairs'(Model);
            Best_Validate_Forcing := new Data_Pairs'(Forcing);
            Best_Validate_M := new Modulations'(M);
            Best_Validate_MAP := new Modulations_Amp_Phase'(MAP);
            Best_Validate_Trend := Trend;
            Best_Validate_Accel := Accel;
            Best_Validate_Aero := Aero;
         end if;
         Threads_Reported := Threads_Reported + 1;
         Am_I_Last := Threads_Reported >= Validate_Thread_Count;
      end Report_Final;

      procedure Winner
        (D : out GEM.LTE.Primitives.Shared.Param_S;
         Model : out Data_Pairs;
         Forcing : out Data_Pairs;
         M : out Modulations;
         MAP : out Modulations_Amp_Phase;
         Trend, Accel, Aero : out Long_Float)
      is
      begin
         D := Best_Validate_D.all;
         Model := Best_Validate_Model.all;
         Forcing := Best_Validate_Forcing.all;
         M := Best_Validate_M.all;
         MAP := Best_Validate_MAP.all;
         Trend := Best_Validate_Trend;
         Accel := Best_Validate_Accel;
         Aero := Best_Validate_Aero;
      end Winner;

      procedure Current
        (Metric : out Long_Float; OOB : out Long_Float;
         Validate_Live : out Long_Float)
      is
      begin
         Metric := Best_Metric;
         OOB := Best_OOB;
         Validate_Live := Best_Validate_Live;
      end Current;

   end Monitor;

   --  Query current best metric (blocking call - waits for improvement)
   --  Returns formatted string with thread ID, iteration count, metrics
   function Status return String is
      Metric, OOB, Validate_Live : Long_Float;
      Client : Integer;
      Cycle : Long_Integer;
      S1, S2, S3 : String (1 .. 10);
   begin
      Monitor.Status (Metric, OOB, Client, Cycle, Validate_Live);
      Ada.Long_Float_Text_IO.Put (S1, Metric, Aft => 5, Exp => 0);
      Ada.Long_Float_Text_IO.Put (S2, OOB, Aft => 5, Exp => 0);
      --  Validate score sits between the training (S1) and held-out test
      --  (S2) readouts, and appears at all only when VALIDATE is on --
      --  its mere presence in the live display is the tell that VALIDATE
      --  is active (see Monitor.Check's own comment: purely informational,
      --  never fed back into any thread's own accept/reject).
      if Validate then
         Ada.Long_Float_Text_IO.Put (S3, Validate_Live, Aft => 5, Exp => 0);
         return "Status: " & Client'Img & Cycle'Img & S1 & " V:" & S3 & S2;
      else
         return "Status: " & Client'Img & Cycle'Img & S1 & S2;
      end if;
   end Status;

   --  Exchange: one per pair of leg tasks (LEGS = 2). The leader posts each
   --  candidate (Propose) before evaluating it on its own leg; the follower
   --  takes it (Fetch), evaluates it on the other leg and posts its score
   --  (Answer). Result hands the leader the composite of the two scores,
   --  weighted by each leg's LEG_WEIGHT, and Verdict hands the follower the
   --  same number, so both tasks make the same accept/reject decision and
   --  keep the same candidate.
   type Flags is array (Follower_Leg) of Boolean;
   type Leg_Scores is array (Follower_Leg) of Long_Float;

   protected type Exchange is
      procedure Propose
        (D : in Shared.Param_S; H : in Ns; Threshold : in Long_Float);
      entry Fetch (Follower_Leg)
        (D : out Shared.Param_S; H : out Ns; Threshold : out Long_Float;
         Quit : out Boolean);
      procedure Answer (L : in Follower_Leg; Score, OOB, Weight : in Long_Float);
      entry Result
        (Score, OOB, Weight : in Long_Float;
         Composite, Composite_OOB : out Long_Float);
      entry Verdict (Follower_Leg) (Composite : out Long_Float);
      --  LEG_COMMON: the leader's regressed amplitudes and phases for this
      --  candidate, for the followers to apply to the common windings.
      procedure Post_Amps (A : in Modulations_Amp_Phase);
      entry Get_Amps (Follower_Leg) (A : out Modulations_Amp_Phase);
      procedure Close;          -- leader has finished: release the followers
      procedure Follower_Gone (L : in Follower_Leg);  -- that follower failed
   private
      Cand : Final_Param_P := null;
      HB : Ns (1 .. 100) := (others => 0);
      Thr : Long_Float := 0.0;
      Have, Answered, Decided, Gone : Flags := (others => False);
      Amps_Ready : Flags := (others => False);
      AB : Modulations_Amp_Phase (1 .. 400) := (others => (0.0, 0.0));
      Closed : Boolean := False;
      Pending : Natural := 0;  -- followers that have not yet answered
      SF, PF, WF : Leg_Scores := (others => 0.0);
      Comp : Long_Float := 0.0;
   end Exchange;

   protected body Exchange is
      procedure Propose
        (D : in Shared.Param_S; H : in Ns; Threshold : in Long_Float) is
      begin
         if Cand = null then
            Cand := new Shared.Param_S'(D);
         else
            Cand.all := D;
         end if;
         HB (1 .. H'Length) := H;
         Thr := Threshold;
         Pending := 0;
         for L in 2 .. Group loop
            if not Gone (L) then
               Have (L) := True;
               Answered (L) := False;
               Pending := Pending + 1;
            end if;
         end loop;
      end Propose;

      entry Fetch (for L in Follower_Leg)
        (D : out Shared.Param_S; H : out Ns; Threshold : out Long_Float;
         Quit : out Boolean)
        when Have (L) or Closed is
      begin
         if Have (L) then
            D := Cand.all;
            H := HB (1 .. H'Length);
            Threshold := Thr;
            Quit := False;
            Have (L) := False;
         else
            Threshold := 0.0;
            Quit := True;
         end if;
      end Fetch;

      procedure Answer (L : in Follower_Leg; Score, OOB, Weight : in Long_Float)
      is
      begin
         SF (L) := Score;
         PF (L) := OOB;
         WF (L) := Weight;
         if not Answered (L) then
            Answered (L) := True;
            Pending := Pending - 1;
         end if;
      end Answer;

      entry Result
        (Score, OOB, Weight : in Long_Float;
         Composite, Composite_OOB : out Long_Float)
        when Pending = 0 is
         W : Long_Float := Weight;
         S : Long_Float := Weight * Score;
         P : Long_Float := Weight * OOB;
      begin
         for L in 2 .. Group loop
            if Answered (L) and then not Gone (L) then
               W := W + WF (L);
               S := S + WF (L) * SF (L);
               P := P + WF (L) * PF (L);
            end if;
         end loop;
         if W > 0.0 then
            Composite := S / W;
            Composite_OOB := P / W;
         else
            Composite := Score;
            Composite_OOB := OOB;
         end if;
         Comp := Composite;
         for L in 2 .. Group loop
            Decided (L) := True;
         end loop;
      end Result;

      entry Verdict (for L in Follower_Leg) (Composite : out Long_Float)
        when Decided (L) or Closed is
      begin
         Composite := Comp;
         Decided (L) := False;
      end Verdict;

      procedure Post_Amps (A : in Modulations_Amp_Phase) is
      begin
         AB (1 .. A'Length) := A;
         Amps_Ready := (others => True);
      end Post_Amps;

      entry Get_Amps (for L in Follower_Leg) (A : out Modulations_Amp_Phase)
        when Amps_Ready (L) or Closed is
      begin
         A := AB (1 .. A'Length);
         Amps_Ready (L) := False;
      end Get_Amps;

      procedure Close is
      begin
         Closed := True;
      end Close;

      procedure Follower_Gone (L : in Follower_Leg) is
      begin
         Gone (L) := True;
         Have (L) := False;
         if not Answered (L) and then Pending > 0 then
            Pending := Pending - 1;
         end if;
      end Follower_Gone;
   end Exchange;

   Max_Pairs : constant := 256;
   Exchanges : array (1 .. Max_Pairs) of Exchange;

   -----------------------------------
   --  Thread body: Continuously runs Dipole_Model optimization
   --  Optional ALTERNATE mode switches between include/exclude training
   --  data on each iteration (for validation testing)
   -----------------------------------
   task body Thread is
      --  Must come first: the settings read below are the leg's own.
      Leg_Assigned : constant Boolean := Assign_Leg (Thread.Leg);
      pragma Unreferenced (Leg_Assigned);
      ID : Integer;
      Name : String :=
        GEM.Getenv (Name => "CLIMATE_INDEX", Default => "nino4.dat");
      Split : Boolean :=
        GEM.Getenv (Name => "SPLIT_TRAINING", Default => False);
      Exclude : Boolean := GEM.Getenv ("EXCLUDE", False);
   begin
      if Paired then
         ID := Thread.Index;  -- group p: leader first, then its followers
         Text_IO.Put_Line
           (Name & " for Thread #" & ID'Img & Thread.CPU'Img & " leg" &
            Thread.Leg'Img & " pair" & Thread.Pair'Img &
            (if Exclude then " excluding" else " training on") & " " &
            GEM.Getenv ("TRAIN_START", "start") & " .. " &
            GEM.Getenv ("TRAIN_END", "end") & " sampling " &
            GEM.Getenv ("SAMPLING", "12.0"));
      else
         Monitor.Client (ID);
         Text_IO.Put_Line (Name & " for Thread #" & ID'Img & Thread.CPU'Img);
      end if;
      accept Start;
      loop
         Dipole_Model
           (N_Tides => Thread.N_Tides, N_Modulations => Thread.N_Modulations,
            ID => ID, File_Name => Name, Split_Training => Split,
            Exclude => Exclude, Leg => Thread.Leg, Pair => Thread.Pair);
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
      Split_Training : in Boolean := False; Exclude : in Boolean := False;
      Leg : in Natural := 0; Pair : in Natural := 0)
   is
      --  LEGS: role of this task in its pair (both False without legs)
      Leader : constant Boolean := Paired and then Leg = 1;
      Follower : constant Boolean := Paired and then Leg >= 2;
      Leader_ID : constant Integer := (if Follower then ID - (Leg - 1) else ID);
      Leg_Weight : constant Long_Float := GEM.Getenv ("LEG_WEIGHT", 1.0);
      --  LEG_IR (set in a leg's own response file, meant for follower legs):
      --  this leg uses the given one-year lag coefficient instead of the
      --  shared, searched IR. The lag term is Model(t) - IR * Model(t - 1 yr),
      --  so LEG_IR 1.0 makes the leg the one-year difference of its response.
      --  LEG_COMMON = H (default 0 = off), for legs that are the same index:
      --  a follower does not regress the harmonics above H itself; it applies
      --  the amplitudes and phases the leader regressed on its own leg
      --  (times this leg's LEG_SCALE), and regresses only the remaining
      --  terms on what is left of its data. Its score then says how well
      --  the leader's fast response holds on this leg.
      Common_Split : constant Integer := GEM.Getenv ("LEG_COMMON", 0);
      Leg_Scale : constant Long_Float := GEM.Getenv ("LEG_SCALE", 1.0);
      --  LEG_SCORE FAST (a follower with LEG_COMMON): the leg is scored on
      --  the leader's fast part alone -- the metric between that part of
      --  the model and the data with the leg's own remaining model removed,
      --  over the leg's training interval -- instead of on the whole model.
      --  It measures whether the fast response calibrated on the leader's
      --  leg is present, relative to the manifold, on this leg.
      Fast_Score : constant Boolean :=
        Follower and then Common_Split > 0
        and then GEM.Getenv ("LEG_SCORE", "FULL") = "FAST";
      Has_Leg_IR : constant Boolean := GEM.Getenv ("LEG_IR", "") /= "";
      Leg_IR : constant Long_Float :=
        (if Has_Leg_IR then GEM.Getenv ("LEG_IR", 0.0) else 0.0);
      Results_File : constant String := "lte_results" & GEM.Leg_Tag & ".csv";
      Follower_Quit : Boolean := False;
      Pair_Threshold : Long_Float := 0.0;
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
      --  Caches the SAME accepted iteration's Forcing (post-Bessel) that
      --  KeepModel/DKeep were derived from -- see the Am_I_Last resync
      --  below (2026-09-29 fix).
      KeepForcing : Data_Pairs := Data_Records;

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
      --  METRIC=FB: band-limited spectral correlation (see FB_CC), cycles/year
      FB_Metric : constant Boolean := GEM.Getenv ("METRIC", "CC") = "FB";
      FB_Lo : constant Long_Float := GEM.Getenv ("FB_LO", 0.5);
      FB_Hi : constant Long_Float := GEM.Getenv ("FB_HI", 1.0);
      FB_Lo2 : constant Long_Float := GEM.Getenv ("FB_LO2", 0.0);
      FB_Hi2 : constant Long_Float := GEM.Getenv ("FB_HI2", 0.0);
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
      --  SAMPLING /= 12 (2026-10-10): data on the grid year + i / SAMPLING.
      --  SPY is the number of samples in the year: the impulse slot wraps at
      --  SPY and the lag term IR reaches back SPY samples (one year). The
      --  integrator's leak mA and friction mP are per MONTH in the .p file;
      --  Step_Scale converts them to per sample, so a monthly fit gives the
      --  same manifold at any sampling. With SAMPLING = 12 nothing changes.
      SPY : constant Integer := Integer (Sampling_Per_Year);
      Monthly : constant Boolean := SPY = 12;
      Step_Scale : constant Long_Float := 12.0 / Sampling_Per_Year;
      Filter : constant Long_Float := GEM.Getenv ("FILTER", 0.333_333_33);
      MLR_On : constant Boolean := GEM.Getenv ("MLR", False); -- wrong name
      Forcing_Only : constant Boolean := GEM.Getenv ("FORCING", False);
      Filter9Pt : constant Integer := GEM.Getenv ("F9", 0);
      Climate_Trend : constant Boolean := GEM.Getenv ("TREND", False);
      NonLin : constant Long_Float := GEM.Getenv ("NONLIN", 1.0);
      Decay : constant Long_Float := GEM.Getenv ("DECAY", 1.0);
      Lock_Freq : constant Boolean := GEM.Getenv ("LOCKF", False);
      --  When set, the DR pre-emphasis below is skipped (regression target
      --  stays raw Data_Records) and Regression_Factors instead folds the
      --  lag-12 delay differential directly into its basis columns -- see
      --  the note at Regression_Factors' Uncompensated_Mode in
      --  gem-lte-primitives.adb for why this is the jointly-optimal
      --  version of what DR only approximates. Same flag also gates the
      --  Save-time modulation-power-spectrum correction.
      Uncompensated : constant Boolean := GEM.Getenv ("UNCOMPENSATED", False);
      Local_Max : constant Boolean := GEM.Getenv ("LOCAL", False);
      Lock_Tidal : constant Boolean := GEM.Getenv ("LOCKT", False);
      Lock_T_Amp : constant Boolean := GEM.Getenv ("LOCKA", False);
      --  ZONAL_ADJ: which of the ZONAL adjustables the search may move, a
      --  list of names from "ztau zdist ztilt ssa" (default none: they stay
      --  at their .p values, so TIDES=ZONAL behaves as before).
      Zonal_Adj : constant String := GEM.Getenv ("ZONAL_ADJ", "");
      Free_ZTau : constant Boolean := Ada.Strings.Fixed.Index (Zonal_Adj, "ztau") > 0;
      Free_ZDist : constant Boolean := Ada.Strings.Fixed.Index (Zonal_Adj, "zdist") > 0;
      Free_ZTilt : constant Boolean := Ada.Strings.Fixed.Index (Zonal_Adj, "ztilt") > 0;
      Free_Ssa : constant Boolean := Ada.Strings.Fixed.Index (Zonal_Adj, "ssa") > 0;
      Sv_ZTau, Sv_ZDist, Sv_ZTilt, Sv_SsaA, Sv_SsaP : Long_Float := 0.0;
      --  IMPULSE_MOD12: see Impulse_Delta (impulse month taken modulo the year).
      --  default on under TIDES=ZONAL, so the semiannual (delA) impulse can fire
      Impulse_Mod12 : constant Boolean := GEM.Getenv ("IMPULSE_MOD12", GEM.Zonal.On);
      --  IMPULSE_EXACT (default FALSE, 2026-10-10): the tide that an impulse
      --  picks up is evaluated at the impulse's own date, year + frac (|delB|)
      --  (and half a year later for the asym impulse), instead of at the date
      --  of the sample the impulse falls on. The impulse still lands on the
      --  nearest sample, but its size no longer depends on how the data are
      --  sampled, so one .p gives the same manifold for monthly, 5-day or
      --  daily data. Delta impulses with the constituent table only.
      Impulse_Exact : constant Boolean := GEM.Getenv ("IMPULSE_EXACT", False);
      --  Constrains ltep (the "winding"/LT modulation period) entries to
      --  snap to the nearest value in a caller-supplied canonical set
      --  (backbone harmonics/subharmonics) after each Markov step that
      --  touches one, rather than letting them drift freely -- unlike
      --  LOCKT/LOCKA above, which freeze the LPAP tidal-constituent
      --  array, this is specific to ltep and leaves everything else
      --  (including which candidate is picked, and every scalar/LPAP
      --  parameter) fully free. See Canon_Candidates below and its
      --  population just after this procedure's own `begin`.
      Lock_Winding : constant Boolean := GEM.Getenv ("LOCKW", False);
      Canon_Backbone : constant Long_Float := GEM.Getenv ("CANON_BACKBONE", 0.207);
      Canon_Harmonics : constant Ns :=
        S_to_I (GEM.Getenv ("CANON_HARMONICS",
                "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 18 19 23 27"));
      Canon_Subharmonic : constant Integer := GEM.Getenv ("CANON_SUBHARMONIC", 11);
      Impulse_Only  : constant Boolean := GEM.Getenv ("IMPULSE", True);
      Has_Friction : constant Boolean := GEM.Getenv ("FRICTION", False);
      Jerk : constant Long_Float := GEM.Getenv ("JERK", 0.0);
      --  Read independently from GEM.LTE.Primitives' own Ridge_Lambda (that
      --  one lives in the parent package BODY, not visible here) — used
      --  only to report the effective value that Regression_Coefficients
      --  actually applied for this run; see the "---- LTE ----" printout.
      Ridge_Lambda : constant Long_Float := GEM.Getenv ("RIDGE", 0.0);
            
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
      --  The calendar date D.B.init's value actually applies at.
      --  Defaults to THIS run's own Data_Records'First -- for the file a
      --  region was originally fit on, that IS init's true meaning (it's
      --  exactly where IIR's Start_Index already lands today, IDATE
      --  clamping or not), so the default reproduces today's output
      --  bit-for-bit with no config needed. Override only when applying
      --  one region's already-fitted manifold to a DIFFERENT (e.g.
      --  longer) file than it was fit on -- set it to that original
      --  file's own first date so `init` keeps its real meaning instead
      --  of being silently reinterpreted as "value at THIS file's first
      --  date", which would be wrong whenever the two files differ.
      --  Read with a sentinel (no real calendar date is negative) purely
      --  to detect whether it was explicitly set at all, for
      --  Strict_IDate's own default just below -- IDATE and INIT_DATE
      --  answer different questions (how far back the physics should
      --  integrate, vs. what date an already-fitted init corresponds to)
      --  and neither can substitute for the other in general: 35 of this
      --  project's own region configs -- every "_"-suffixed full-record
      --  companion (kN_baltic, kN020_E050_, kS040_W050_, kRedSea,
      --  kPersianGulf, ...) -- already have IDATE falling INSIDE their
      --  own record rather than before it, so their own init is already
      --  calibrated to mean "value near IDATE", not "value at the
      --  record's own start"; unconditionally defaulting Init_Date to
      --  Data_Records'First for everyone would silently break those.
      Init_Date_Raw : constant Long_Float := GEM.Getenv ("INIT_DATE", -1.0);
      Init_Date_Explicit : constant Boolean := Init_Date_Raw >= 0.0;
      Init_Date : constant Long_Float :=
        (if Init_Date_Explicit then Init_Date_Raw
         else Data_Records (Data_Records'First).Date);
      --  Default False: every existing .resp (none of which set this key)
      --  keeps today's behavior byte-for-byte -- IIR's own Start_Index
      --  search silently clamps to the record's own first sample whenever
      --  IDATE precedes it, which is what every currently-fitted `init`
      --  value's meaning implicitly assumes. Opting a NEW fit into True
      --  makes Calc_Forcing genuinely integrate from IDATE via
      --  Extend_Backward first -- correct and now safe to do (IIR's
      --  backward pass is a verified exact inverse, see
      --  iir_invariant_test.adb). Does NOT require re-fitting `init`:
      --  its real, already-calibrated meaning is "value at Init_Date"
      --  above, not "value at IDATE" -- Calc_Forcing seeds the extended
      --  computation there, unchanged, and lets the exact backward pass
      --  derive everything from IDATE to Init_Date on its own,
      --  deterministically, from the SAME init every existing fit
      --  already has. Also implied by Init_Date_Explicit: setting
      --  INIT_DATE only ever makes sense together with the extended
      --  computation, so it turns Strict_IDate on by itself -- one flag,
      --  not two, for the cross-file case.
      Strict_IDate : constant Boolean :=
        GEM.Getenv ("STRICT_IDATE", False) or else Init_Date_Explicit;
      RMS_Data : Long_Float := 0.0;

      --  Quadrature-flip regularizer (2026-09-29, per user request): the
      --  lowest-value winding (~0.01-0.02, BACKBONE/11 subharmonic family)
      --  has its amplitude AND phase decided purely by OLS regression --
      --  its period rivals the record length, so sin/cos are nearly
      --  collinear and the regression can land on either quadrature with
      --  no algebraic way to choose (documented at length in HANDOFF.md,
      --  "subharmonic quadrature metastability"). MANIFOLD points at a
      --  reference manifold file (a mean/median manifold across the
      --  collective, written by a sweep script, NOT read/written by this
      --  program) sampled on the SAME date grid as this cell's own
      --  Data_Records -- every Metric call's own Z argument is always
      --  some slice of Forcing (the manifold, lte_results.csv column 4)
      --  taken from that same array's own bounds, so Manifold_Ref can be
      --  indexed by Z'First..Z'Last directly with no separate date-
      --  matching step. Off (no behavior change) whenever MANIFOLD is
      --  unset or the file doesn't exist -- deliberately checked with
      --  Ada.Directories.Exists BEFORE ever calling Make_Data, since
      --  Make_Data's own failure path is GNAT.OS_Lib.OS_Exit(0) (silent
      --  process termination, not a catchable exception) for a missing
      --  file, which would be far too destructive for what's meant to be
      --  an optional regularizer.
      Manifold_Ref_File : constant String := GEM.Getenv ("MANIFOLD", "");
      Manifold_Regularize : constant Boolean :=
        Manifold_Ref_File'Length > 0
        and then Ada.Directories.Exists (Manifold_Ref_File);
      Manifold_Ref : constant Data_Pairs :=
        (if Manifold_Regularize then Make_Data (Manifold_Ref_File)
         else Empty_Data);
      --  TIDES=HYBRID: ridge weight on the table's correction lines (see GEM.Zonal)
      Hybrid_Lambda : constant Long_Float := GEM.Getenv ("HYBRID_LAMBDA", 1.0);
      --  LOD_TERM: mixes a little LOD (time integral of the dLOD-rate forcing)
      --  into the climate forcing, amplitude-corrected per line (see GEM.Zonal)
      Lod_Term : constant Long_Float := GEM.Getenv ("LOD_TERM", 0.0);

      function Metric (X, Y, Z : in Data_Pairs) return Long_Float is
         Raw_Score : Long_Float;
      begin
         if RMS_Metric then
            Raw_Score := RMS (X, Y, RMS_Data, 0.0);
            --  TODO: Can remove - Experimental hybrid metric combining RMS and CC
            --  was tested but never adopted. Simple RMS alone was sufficient.
            -- Raw_Score := (RMS(X,Y,RMS_Data, 0.0) + CC(X,Y))/2.0;
         elsif CID_Metric then
            Raw_Score := CC (X, Y) * CID (X, Y);
         elsif ZC_Metric then
            Raw_Score := Xing (X, Y);
         elsif DER_Metric then
            Raw_Score := DER_CC (X, Y);
         elsif FT_Metric then
            Raw_Score := FT_CC (X, Y, Z);
         elsif FB_Metric then
            Raw_Score := FB_CC (X, Y, FB_Lo, FB_Hi, FB_Lo2, FB_Hi2);
         elsif Winding_Metric then
            Raw_Score := (Winding_Agreement (X, Y, Z) + CC(X,Y))/2.0;
         elsif DTW_Metric then
            Raw_Score := DTW_Distance (X, Y, DTW_Window);
         elsif EMD_Metric then
            Raw_Score := EMD (X, Y);
         elsif Hoy_Metric then
            Raw_Score := -- Ada.Numerics.Long_Elementary_Functions.Sqrt ()
                   (Hoyer_Spectral_Peak (X, Y, Z) + CC (X, Y)) * 0.5;
         elsif DTW_CC then
            Raw_Score :=
              Ada.Numerics.Long_Elementary_Functions.Sqrt
                (DTW_Distance (X, Y, DTW_Window) * CC (X, Y));
         elsif CTW_Metric then
            Raw_Score :=
              Ada.Numerics.Long_Elementary_Functions.Sqrt
                (Long_Float'Max (DTW_Distance (X, Y, DTW_Window), 0.0) *
                 CID (X, Y));
         elsif EMD_CC then
            Raw_Score := (EMD (X, Y) + CC (X, Y)) * 0.5;
         elsif EMD_DER then
            Raw_Score := (EMD (X, Y, Derivative => True) + CC (X, Y)) * 0.5;
         elsif SEM_Metric then
            Raw_Score := Scaled_Error_Metric (X, Y);  -- Y should be data
         elsif MLR_On then
            Raw_Score := CC (X, Y) * Min_Entropy_Power_Spectrum (Z, Y);
         elsif Is_Minimum_Entropy then
            Raw_Score := Min_Entropy_Power_Spectrum (X, Y); -- or X = Z
         else
            Raw_Score := CC (X, Y);
         end if;

         --  Multiplicative penalty: Manifold_CC~1.0 (this run's own
         --  manifold slice agrees with the collective reference) leaves
         --  Raw_Score untouched; Manifold_CC~0.0 (orthogonal, i.e. the
         --  classic quadrature-swap signature) or negative (anti-
         --  correlated) crushes it toward/through zero, so a thread's own
         --  accept/reject (Old_CC, Monitor.Check, all downstream of this
         --  single choke point every Metric caller already goes through)
         --  naturally steers away from the wrong quadrature without
         --  needing a second, separate search objective.
         if Manifold_Regularize then
            begin
               declare
                  Manifold_CC : constant Long_Float :=
                    CC (Z, Manifold_Ref (Z'First .. Z'Last));
               begin
                  Raw_Score := Raw_Score * Long_Float'Max (0.0, Manifold_CC);
               end;
            exception
               when Constraint_Error =>
                  --  Manifold_Ref doesn't cover this particular Z slice
                  --  (e.g. a shorter/misaligned reference file) -- no
                  --  reference available for this range, so no penalty.
                  null;
            end;
         end if;

         --  TIDES=HYBRID: penalize the correction lines by their share of
         --  the generator's variance, so a line only grows if it pays.
         if GEM.Zonal.Hybrid and then Hybrid_Lambda > 0.0 then
            declare
               S : Long_Float := 0.0;
            begin
               for I in D.B.LPAP'Range loop
                  S := S + 0.5 * D.B.LPAP (I).Amplitude**2;
               end loop;
               Raw_Score := Raw_Score * Long_Float'Max
                 (0.0, 1.0 - Hybrid_Lambda * S / GEM.Zonal.Forcing_Variance);
            end;
         end if;

         return Raw_Score;
      end Metric;

      function Impulse_Delta (Time : Long_Float) return Long_Float is
         Value : Long_Float;
         -- Impulses will occur on a month for monthly data
         Trunc : Integer :=
           Integer (((Time - Long_Float'Floor (Time)) * Sampling_Per_Year));
         --  IMPULSE_MOD12=TRUE: reduce the impulse month modulo the year. Without it,
         --  |DelB| >= ~0.96 gives DPos >= 12, which no month index (0..11) can equal,
         --  so the DelA impulse never fires and only the half-year companion (Asym,
         --  at (DPos + 6) mod 12) does. Default FALSE keeps existing fits as fitted.
         DPos : Integer :=
           (if Impulse_Mod12
            then Integer ((abs (D.B.DelB) * Sampling_Per_Year)) mod Integer (Sampling_Per_Year)
            else Integer ((abs (D.B.DelB) * Sampling_Per_Year)));
      begin
         if Trunc = DPos then
            Value := D.B.DelA;
         elsif Trunc = (DPos + Integer (Sampling_Per_Year) / 2) mod SPY then
            Value := D.B.Asym;
         else
            Value := 0.0;
         end if;
         return Value;
      end Impulse_Delta;

      --  Annual_Impulse: a genuine once-a-year DELTA (not a smooth Ann1/
      --  Ann2 sinusoid) with its own independent amplitude, D.B.ImpC --
      --  previously a dead pass-through (threaded into Bessel's k2
      --  parameter, which the active Bessel body never actually reads).
      --  Reuses Impulse_Delta's own DPos month-of-year slot (from DelB)
      --  rather than adding a second free phase parameter -- this tests
      --  specifically whether a sharp annual kick coincident with the
      --  existing tidal-gating impulse explains real variance the smooth
      --  seasonal terms don't, not a independently-phased new cycle.
      --  Motivated directly by kN020_E050: Ann1/Ann2/Sem1/Sem2 together
      --  account for only ~10% of the real data's peak-to-peak
      --  excursion, so the seasonal-cycle SHAPE the smooth harmonics
      --  can't reach is a real, unexplained gap worth testing an
      --  impulsive (not sinusoidal) model against. Inert by construction
      --  when ImpC = 0.0 (the default for every existing fit).
      --
      --  DPos is taken mod Sampling_Per_Year -- Impulse_Delta's own DPos
      --  (reused here) is NOT range-limited to a valid Trunc slot
      --  (0..Sampling_Per_Year-1) and a fitted |DelB|*Sampling_Per_Year
      --  routinely lands above that (e.g. kN020_E050's own DelB=1.126
      --  gives a raw DPos of 13, matching no Trunc value at all, which
      --  would make Annual_Impulse permanently a no-op for that region
      --  regardless of ImpC). Impulse_Delta's own second branch already
      --  wraps its analogous position with `mod 12` for exactly this
      --  reason; this mirrors that fix for the primary slot instead of
      --  inheriting the same latent bug.
      function Annual_Impulse (Time : Long_Float) return Long_Float is
         Trunc : Integer :=
           Integer (((Time - Long_Float'Floor (Time)) * Sampling_Per_Year));
         DPos : Integer :=
           Integer ((abs (D.B.DelB) * Sampling_Per_Year)) mod
             Integer (Sampling_Per_Year);
      begin
         if Trunc = DPos then
            return D.B.ImpC;
         else
            return 0.0;
         end if;
      end Annual_Impulse;

      function Impulse_Delta_Smear (Time : Long_Float) return Long_Float is
         Value : Long_Float;
         -- Impulses will occur on a month for monthly data
         Trunc : Integer :=
           Integer (((Time - Long_Float'Floor (Time)) * Sampling_Per_Year));
         --  IMPULSE_MOD12=TRUE: reduce the impulse month modulo the year. Without it,
         --  |DelB| >= ~0.96 gives DPos >= 12, which no month index (0..11) can equal,
         --  so the DelA impulse never fires and only the half-year companion (Asym,
         --  at (DPos + 6) mod 12) does. Default FALSE keeps existing fits as fitted.
         DPos : Integer :=
           (if Impulse_Mod12
            then Integer ((abs (D.B.DelB) * Sampling_Per_Year)) mod Integer (Sampling_Per_Year)
            else Integer ((abs (D.B.DelB) * Sampling_Per_Year)));
         DPos2 : Integer := (DPos + Integer (Sampling_Per_Year) / 2) mod SPY;
      begin
         if Trunc = DPos2-1 then
            if D.B.DelA < 0.0 then
               Value := -D.B.DelA;
            else
               Value := 0.0;
            end if;
         elsif Trunc = Dpos2 then
            Value := D.B.Asym;
         elsif Trunc = DPos2+1 then
            if D.B.DelA > 0.0 then
               Value := D.B.DelA;
            else
               Value := 0.0;
            end if;
         else
            Value := 0.0;
         end if;
         return Value;
      end Impulse_Delta_Smear;

      function Impulse (Time : Long_Float) return Long_Float is
         Value : Long_Float;
      begin
         if Impulse_Only then
            Value := Impulse_Delta (Time);
         else
            Value := Impulse_Delta_Smear (Time);
         end if;
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
                                     eS2*eS*Sin(4.0 * Pi * k * M (I).Value) + 
                                     eC2*eC*Cos(4.0 * Pi * k * M (I).Value);
                                     --eS2*eS*Sin(2.0 * Pi * k2* M (I).Value) + 
                                     --eC2*eC*Cos(2.0 * Pi * k2* M (I).Value);
         end loop;
         return M;
      end Bessel;

      --  LAYER (default 0 = off): two-stage manifold. The LOWEST winding
      --  M(1) is taken out of the MLR and instead builds a new manifold
      --  from the tidal one F (after its usual first-order Bessel
      --  correction by impA/impB, which keep their meaning),
      --     LAYER=1:  F2 = eS*sin(2 pi M(1) F) + eC*cos(2 pi M(1) F)
      --     LAYER=2:  F2 = F + the same term
      --  and the remaining windings M(2..) are regressed on F2. eS/eC are
      --  the offs/bg search parameters (otherwise only the Bessel stage's
      --  second-order factors, dropped in this mode), so the lowest
      --  winding's amplitude and phase (and M(1) itself) are found
      --  iteratively by the search, not by the regression. Needs NM >= 2.
      --  LAYER_MIX=TRUE (default FALSE) with LAYER=1 mixes back a weight W
      --  of the Bessel-corrected tidal manifold, F2 = W*Fb + (layer term),
      --  with W = impC (otherwise the once-a-year model impulse, which
      --  contributed nothing measurable); W = 0 reproduces the pure layer.
      --  Opt-in so that existing fits, which carry small nonzero impC
      --  values as impulses, keep their meaning. LAYER=2 is W = 1 and
      --  keeps impC as the annual impulse.
      Layer_Mix : constant Boolean := GEM.Getenv ("LAYER_MIX", False);
      --  LAYER_HP (months, default 0 = off): the mixed-back Fb is first
      --  HIGH-PASSED, Fb minus its centred LAYER_HP-month running mean
      --  (shortened at the ends), so W restores the fast tidal content
      --  the layer suppresses at its crest without tilting the slow
      --  (~120-yr) part that sets the layer's quantized levels.
      Layer_HP : constant Integer := GEM.Getenv ("LAYER_HP", 0);

      --  IR_SIDE (default 0 = symmetric): one-sided 12-month delay
      --  differential. +1 applies the IR term only where the value a year
      --  earlier was WARM, -1 only where it was COLD; "warm" means above
      --  the series' own centred IR_BASE-month running mean (default 121,
      --  so the trend does not decide the sign). Applied consistently to
      --  the data side (DR = D + IR*D(t-12)) and the model side
      --  (Model - IR*Model(t-12)), each judged on its own series.
      --  Compensated mode only (UNCOMPENSATED folds IR into the regression
      --  basis symmetrically, so IR_SIDE is ignored there).
      IR_Side : constant Integer := GEM.Getenv ("IR_SIDE", 0);

      --  ANNUAL_DITHER a (default 0 = off), ANNUAL_DITHER_PHASE phi
      --  (radians), ANNUAL_DITHER_STAGE (1 = after the Bessel/layer
      --  transform, on the manifold the windings see; 0 = before it, on
      --  the tidal manifold): add a*cos(2 pi t + phi) to that manifold. The
      --  model output and the regression columns are then de-seasonalized
      --  (Remove_Climatology), so only non-annual nonlinear products of
      --  the dummy annual signal reach the fit.
      Dither_A : constant Long_Float := GEM.Getenv ("ANNUAL_DITHER", 0.0);
      Dither_Phase : constant Long_Float :=
        GEM.Getenv ("ANNUAL_DITHER_PHASE", 0.0);
      Dither_Stage : constant Integer := GEM.Getenv ("ANNUAL_DITHER_STAGE", 1);
      function Add_Dither (F : Data_Pairs) return Data_Pairs is
         R : Data_Pairs := F;
         use Ada.Numerics.Long_Elementary_Functions;
      begin
         for I in F'Range loop
            R (I).Value := F (I).Value + Dither_A *
              Cos (2.0 * Ada.Numerics.Pi * F (I).Date + Dither_Phase);
         end loop;
         return R;
      end Add_Dither;
      IR_Base : constant Integer := GEM.Getenv ("IR_BASE", 121);

      --  True where X(I) is on the side IR_SIDE selects (always True for 0)
      type Side_Mask is array (Integer range <>) of Boolean;
      function IR_Mask (X : Data_Pairs) return Side_Mask is
         R : Side_Mask (X'Range) := (others => True);
         Cum : array (X'First - 1 .. X'Last) of Long_Float;
         H : constant Integer := Integer'Max (1, IR_Base / 2);
         A, B : Integer;
         Mean : Long_Float;
      begin
         if IR_Side = 0 or else Uncompensated then
            return R;
         end if;
         Cum (X'First - 1) := 0.0;
         for I in X'Range loop
            Cum (I) := Cum (I - 1) + X (I).Value;
         end loop;
         for I in X'Range loop
            A := Integer'Max (X'First, I - H);
            B := Integer'Min (X'Last, I + H);
            Mean := (Cum (B) - Cum (A - 1)) / Long_Float (B - A + 1);
            R (I) := (if IR_Side > 0 then X (I).Value > Mean
                      else X (I).Value < Mean);
         end loop;
         return R;
      end IR_Mask;
      Layer : constant Integer := GEM.Getenv ("LAYER", 0);
      --  LAYER_K (default 0 = off): wavenumber of the impA/impB Bessel
      --  correction in LAYER mode. By default that is M(NM) (M(NM-1)
      --  with LOCKF), i.e. it moves when NM changes; reducing NM to a
      --  single regressed winding would then re-key the correction to
      --  that winding and change the layered manifold. LAYER_K pins it
      --  (e.g. to the old backbone).
      --  With NM >= 3 (and LAYER_K unset) the Bessel wavenumber is instead
      --  the SECOND winding M(2), which then also leaves the MLR: M(1)
      --  builds the layer, M(2) keys the Bessel correction, M(3..) (and
      --  integer harmonics of M(NM)) are regressed -- all three found
      --  independently by the search.
      Layer_K : constant Long_Float := GEM.Getenv ("LAYER_K", 0.0);

      function Layered
        (Model : in Data_Pairs; eS, eC, k : in Long_Float;
         W_Mix : in Long_Float := 0.0) return Data_Pairs
      is
         M : Data_Pairs := Model;
         Pi : constant Long_Float := Ada.Numerics.Pi;
         use Ada.Numerics.Long_Elementary_Functions;
         W : Long_Float;
         --  Mixed-back part: Fb itself, or Fb minus its running mean.
         Mix : array (Model'Range) of Long_Float;
      begin
         if Layer = 1 and then W_Mix /= 0.0 and then Layer_HP > 1 then
            declare
               Cum : array (Model'First - 1 .. Model'Last) of Long_Float;
               H : constant Integer := Layer_HP / 2;
               A, B : Integer;
            begin
               Cum (Model'First - 1) := 0.0;
               for I in Model'Range loop
                  Cum (I) := Cum (I - 1) + Model (I).Value;
               end loop;
               for I in Model'Range loop
                  A := Integer'Max (Model'First, I - H);
                  B := Integer'Min (Model'Last, I + H);
                  Mix (I) := Model (I).Value
                    - (Cum (B) - Cum (A - 1)) / Long_Float (B - A + 1);
               end loop;
            end;
         else
            for I in Model'Range loop
               Mix (I) := Model (I).Value;
            end loop;
         end if;
         for I in Model'Range loop
            W := eS * Sin (2.0 * Pi * k * M (I).Value)
               + eC * Cos (2.0 * Pi * k * M (I).Value);
            M (I).Value :=
              (if Layer = 2 then M (I).Value + W else W_Mix * Mix (I) + W);
         end loop;
         return M;
      end Layered;

      function Frictional (Model : in Data_Pairs;
                           Offset, Saturate, Eps : in Long_Float) return Data_Pairs
      is
         M : Data_Pairs := Model;
         use Ada.Numerics.Long_Elementary_Functions;
         k : Long_Float := 0.0001 + abs Saturate;
         Y : Long_Float;
      begin
         for I in Model'Range loop
            --Y := abs (M (I).Value - abs Offset);
            Y := sqrt( (abs(M (I).Value - abs Offset))**2.25 + abs Eps );
            --Y := sqrt (abs ( (M (I).Value - abs Offset) )) + abs Eps;
            --Y := abs (M (I).Value - abs Offset) + abs Eps;
            M (I).Value := 1.0/k*(1.0 - exp(-k*Y));
         end loop;
         return M;
      end Frictional;

      
      --  Val_Validate sits between Val1 (train) and Val2 (test), and prints
      --  at all only when VALIDATE is on (its presence is itself the tell)
      --  -- mirrors Status's own live "V:" column, but this is the final,
      --  once-only report of the score that actually decided the save.
      procedure Put_CC
        (Val1, Val2 : in Long_Float; Counter : in Long_Integer;
         Thread : in Integer; Val_Validate : in Long_Float := 0.0)
      is
      begin
         Text_IO.Put (GEM.Getenv ("METRIC", "CC"));
         Ada.Long_Float_Text_IO.Put (Val1, Fore => 4, Aft => 10, Exp => 0);
         if Validate then
            Text_IO.Put (" V:");
            Ada.Long_Float_Text_IO.Put
              (Val_Validate, Fore => 4, Aft => 10, Exp => 0);
         end if;
         Ada.Long_Float_Text_IO.Put (Val2, Fore => 4, Aft => 10, Exp => 0);
         Text_IO.Put_Line ("  " & Thread'Img & Counter'Img);
      end Put_CC;

      der : Long_Float;
      CorrCoeff, Old_CC, Prior_Best_CC : Long_Float := 0.0;
      --  True once this thread's search has accepted at least one
      --  candidate (CorrCoeff > Old_CC, below) -- i.e. KeepModel/DKeep
      --  have been updated away from their own declaration-time
      --  initializers (a plain copy of the real Data_Records/D0, kept
      --  only to give them the right bounds/shape before anything real
      --  exists). If this never happens across the whole run, see the
      --  deadlock guard right after the main loop.
      Ever_Accepted : Boolean := False;
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

      --  COVERAGE (0.0 .. 1.0, default 1.0 = current unregularized
      --  behavior): standalone alternative to SPLIT_TRAINING/RATIO. The
      --  regression fits only the leading Coverage fraction of the training
      --  interval [First, Cov_Last]; the resulting fit is then scored over
      --  the WHOLE training interval [First, Last] and that whole-interval
      --  score directly gates each thread's own accept/reject (unlike
      --  SPLIT_TRAINING's CorrCoeffTest / RATIO's OOB, which only affect
      --  cross-thread selection, never an individual thread's own search).
      --  A parameter set that only fits the coverage slice by chasing its
      --  noise will predict the withheld remainder poorly and get rejected
      --  right here, not just deprioritized after the fact. The true
      --  out-of-band interval (outside [First,Last] entirely) is untouched
      --  either way and still reported separately for final validation.
      Coverage : constant Long_Float := GEM.Getenv ("COVERAGE", 1.0);
      Cov_Last : constant Integer :=
        Integer'Min
          (Last, First + Integer (Long_Float (Last - First) * Coverage));

      function Coverage_Region (D : Data_Pairs) return Data_Pairs is
      begin
         return D (First .. Cov_Last);
      end Coverage_Region;

      --  ENCLOSING (default False): standalone alternative to COVERAGE,
      --  for the opposite geometry — training data that BRACKETS the
      --  [TRAIN_START,TRAIN_END] gap on both sides rather than a window
      --  the gap sits outside of. No fraction parameter is needed (unlike
      --  COVERAGE) because the split is already fully determined by
      --  First/Last: fit on LOWER only (Data'First .. First — the data
      --  before the gap); accept/reject is driven by ENCLOSED (below);
      --  and the gap itself, [First,Last], is always reported as the true
      --  held-out validation — untouched by both fitting and accept/
      --  reject in this mode, unlike Exclude_Metric's LOWER+UPPER
      --  combination which is no longer "held out" once ENCLOSING is
      --  using both flanks for the search. Takes precedence over
      --  EXCLUDE/SPLIT_TRAINING the same way COVERAGE does — the two
      --  standalone mechanisms assume opposite fit geometries and aren't
      --  meant to combine, so COVERAGE (checked first) wins if both are
      --  somehow set.
      Enclosing : constant Boolean := GEM.Getenv ("ENCLOSING", False);

      --  ENCLOSED (default "UL"): which region(s) drive ENCLOSING's
      --  accept/reject, mirroring SPLIT_TRAINING/SPLIT_LOW's pattern of a
      --  companion selector next to the activating flag.
      --    "UL" (default) — UPPER + LOWER combined (Exclude_Metric): the
      --      search is rewarded for fitting LOWER *and* predicting UPPER
      --      together, so it can't drift away from LOWER's pattern
      --      without being penalized — matches how COVERAGE scores the
      --      whole interval, not just the withheld remainder.
      --    "U" — UPPER only: a stricter, pure-generalization signal (the
      --      original behavior) — the search is judged solely on
      --      predicting UPPER, with no direct penalty for drifting away
      --      from LOWER's own pattern (the regression's own best-fit
      --      mechanics are the only thing keeping it anchored to LOWER).
      Enclosed_Upper_Only : constant Boolean :=
        GEM.Getenv ("ENCLOSED", "UL") = "U";

      function Lower_Region (D : Data_Pairs) return Data_Pairs is
      begin
         return D (D'First .. First);
      end Lower_Region;

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

      --  VALIDATE lockbox: the genuinely-untouched-by-training region --
      --  mirrors, independently (never by reusing CorrCoeffP itself, so
      --  this can't leak into Monitor.Check's OOB-weighted M), whichever
      --  region each mode's own CorrCoeffP-final assignment above already
      --  treats as the true held-out set. Precedence matches that same
      --  per-iteration branching exactly (COVERAGE first, unaffected by
      --  EXCLUDE; then ENCLOSING; then SPLIT_TRAINING; then plain).
      --
      --  CORRECTNESS NOTE (bug fixed here): the first version of this
      --  function always scored [Last, D'Last], on the unstated
      --  assumption that region is always held out. That's only true
      --  when EXCLUDE=FALSE. Under EXCLUDE=TRUE (and ENCLOSING),
      --  training itself fits the concatenation of BOTH outer flanks --
      --  [D'First,First] AND [Last,D'Last] -- so [Last,D'Last] is
      --  in-sample there, and the genuinely held-out region is the
      --  [First,Last] gap instead. Scoring the wrong region silently
      --  reports in-sample fit quality as if it were a validation score
      --  (confirmed live on a real EXCLUDE=TRUE run: "validate" tracked
      --  the true held-out gap score far more closely than it should
      --  have if it were really an independent, untouched check).
      --
      --  Takes the model to score explicitly (the thread's own KEPT best
      --  candidate, KeepModel -- NOT the live Model, which may reflect a
      --  since-rejected perturbation) rather than defaulting to Model.
      --  Below this many points, a correlation coefficient is not just
      --  noisy but MATHEMATICALLY DEGENERATE: Pearson CC on exactly 2
      --  points is a line through 2 points, which is ALWAYS exactly +-1
      --  (or 0.0 in CC's own zero-variance guard case) regardless of fit
      --  quality -- confirmed live: with no TRAIN_START/TRAIN_END set
      --  (First=D'First, Last=D'Last, the common default), Outer_Flanks'
      --  own "two 1-point flanks concatenated" is exactly this 2-point
      --  case, and hit exactly 1.0 on the very first report of a run.
      --  Harmless for Exclude_Metric's OWN existing use (CorrCoeffP is
      --  freshly recomputed every iteration, so a stray +-1.0 just
      --  flickers), but genuinely corrupting for Validate_Metric: its
      --  score feeds a running MAXIMUM (Monitor's Best_Validate_Live) and
      --  the once-only lockbox winner comparison (Report_Final) -- a
      --  single lucky +-1.0 hit gets locked in forever as "the best ever"
      --  in the former, and can win the save outright in the latter,
      --  with zero relation to real generalization.
      Min_Validate_Segment : constant := 24; -- 2 years of monthly data

      --  REVERTED to the ORIGINAL spec ("uses the 2nd outer as a
      --  validation set... instead of concatenating everything, the 2nd
      --  outer interval is used"): ALWAYS score [Last, D'Last] alone --
      --  not mode-dependent (no Coverage/Enclosing/Split_Training/Exclude
      --  branching, no Outer_Flanks concatenation). A prior revision made
      --  this Exclude-aware (switching to the [First,Last] gap under
      --  EXCLUDE=TRUE), reasoning that [Last,D'Last] is technically part
      --  of Exclude=TRUE's own combined training objective there -- true,
      --  but it silently overrode the explicit original design: a fixed,
      --  always-the-same-region check on the most RECENT real data,
      --  precisely so a strong post-training-window score there (even one
      --  the joint training objective also benefits from) is something
      --  the user can see and act on. Confirmed directly: the display was
      --  correctly reporting the [First,Last] EXCLUDE gap's score, not a
      --  display/store bug -- just the wrong region for what was wanted.
      function Validate_Metric (Scored_Model : Data_Pairs) return Long_Float is
      begin
         if Scored_Model'Last - Last + 1 < Min_Validate_Segment then
            return Long_Float'First; -- too degenerate to trust; see above
         end if;
         return Metric
           (Scored_Model (Last .. Scored_Model'Last),
            Data_Records (Last .. Data_Records'Last),
            Forcing (Last .. Forcing'Last));
      end Validate_Metric;

      function Excluded (D : Data_Pairs) return Data_Pairs is
      begin
         if Exclude then
            return D (D'First .. First) & D (Last .. D'Last);
         else
            return D (First .. Last);
         end if;
      end Excluded;

      --  The region Regression_Factors is handed (same choice as its
      --  call site); used to report the fit span. Under EXCLUDE it is the
      --  record minus the TRAIN_START..TRAIN_END test gap.
      --  DEBUG_KEEP (env, default off): full-precision fingerprints of a
      --  candidate, printed at every accept and again at save, to check
      --  that what is saved is what was kept (2026-10-01 investigation).
      Debug_Keep : constant Boolean := GEM.Getenv ("DEBUG_KEEP", False);

      procedure Fingerprint
        (Tag : String; P : GEM.LTE.Primitives.Shared.Param_S;
         Mdl, Frc : Data_Pairs; Mv : Modulations; H : Ns;
         T, Ae, CCv : Long_Float; Who : Integer)
      is
         SM, SF, SK, SL : Long_Float := 0.0;
         function HS (X : Ns) return String is
           (if X'Length = 0 then ""
            else X (X'First)'Img & HS (X (X'First + 1 .. X'Last)));
      begin
         if not Debug_Keep then
            return;
         end if;
         for I in Mdl'Range loop
            SM := SM + Mdl (I).Value;
            SF := SF + Frc (I).Value;
         end loop;
         for I in Mv'Range loop
            SK := SK + Mv (I);
         end loop;
         for I in P.B.LPAP'Range loop
            SL := SL + P.B.LPAP (I).Amplitude + P.B.LPAP (I).Phase;
         end loop;
         Text_IO.Put_Line
           ("KEEPDBG " & Tag & Who'Img & " cc" & CCv'Img & " sumModel" & SM'Img
            & " sumForcing" & SF'Img & " sumM" & SK'Img & " trend" & T'Img
            & " aero" & Ae'Img & " level" & P.A.level'Img & " k0" & P.A.k0'Img
            & " init" & P.B.init'Img & " IR" & P.B.IR'Img & " impC" & P.B.ImpC'Img
            & " mP" & P.B.mP'Img & " year" & P.B.Year'Img & " sumLPAP" & SL'Img
            & " harm" & HS (H));
      end Fingerprint;

      function Fit_Region (D : Data_Pairs) return Data_Pairs is
        (if Coverage < 1.0 then Coverage_Region (D)
         elsif Enclosing then Lower_Region (D)
         else Excluded (D));

      Max_Harmonics : Positive := GEM.Getenv ("MAXH", 1_000);

      ------------------------------------------------------------------------
   -- Array overlay for random descent optimization (Type-Safe Version)
   -- The Param_B record is overlaid with an array to allow Walker.Markov
   -- to randomly modify parameters. This is documented and verified via
   -- the Param_B_Overlay package for maximum safety.
      ------------------------------------------------------------------------
      NM : constant Integer := GEM.Getenv ("NM", N_Modulations);
      Minimum_Modulation : constant Long_Float := GEM.Getenv ("MIN_LT", 0.01);
      
      -- Full overlay size for verification (entire D.B structure)
      Full_Size : constant Positive :=
        GEM.LTE.Primitives.Param_B_Overlay.Overlay_Size (D.B.NLP, D.B.NLT);
      
      -- Reduced size for Walker search based on NM (only search used LT entries)
      -- Size = 25 scalars + (NLP * 2) LPAP + NM LT entries
      Size_Shared : constant Positive :=
        GEM.Getenv ("DSIZE", 25 + (D.B.NLP * 2) + NM);

      package Walker is new GEM.Random_Descent
        (Fixed => Is_Fixed, Set_Range => Size_Shared,
         Harmonic_Range => Max_Harmonics);

      Set, Keep, Set0 : Walker.LF_Array (1 .. Size_Shared);

      --  ltep occupies exactly the last NM entries of Set (matches
      --  Param_B_Overlay.First_LT_Index (D.B.NLP) = 25 + NLP*2 =
      --  Size_Shared - NM + 1, confirmed against that package). Only
      --  used when Lock_Winding is set -- see the population loop right
      --  after this procedure's `begin` and the snap call after each
      --  Markov step in the search loop below.
      LT_Lo : constant Positive := Size_Shared - NM + 1;
      LT_Hi : constant Positive := Size_Shared;
      --  Both signs of backbone*harmonic and backbone/subharmonic, since
      --  ltep can be negative (a phase-convention flip, confirmed seen
      --  in real fits this session) and Nearest_In_Set has no notion of
      --  sign-preserving comparison on its own.
      Canon_N : constant Positive := Canon_Harmonics'Length + 1;
      Canon_Candidates : Walker.LF_Array (1 .. 2 * Canon_N);

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
      --  Date Regression_Factors centred Accel on (its own return value,
      --  never re-derived here); LTE and the windings JSON use it as-is.
      Accel_Ref : Long_Float := 0.0;
      Aero : Long_Float := 0.0;  -- AERO coefficient (0 when AERO unset)
      --  ALPHA (default 0.0 = off): winding terms are doubled while
      --  F - M_Min <= Alpha*|M_Min|. M_Min is the minimum of the WHOLE
      --  manifold for the candidate being evaluated (not the training rows
      --  only), computed once and handed to both Regression_Factors and
      --  LTE so the fit and the evaluation use the same envelope.
      Alpha : constant Long_Float := GEM.Getenv ("ALPHA", 0.0);
      M_Min : Long_Float := 0.0;

      function Manifold_Min (F : Data_Pairs) return Long_Float is
         V : Long_Float := Alpha_Signal (F, F'First);
      begin
         for I in F'Range loop
            V := Long_Float'Min (V, Alpha_Signal (F, I));
         end loop;
         return V;
      end Manifold_Min;
      --  The accepted (kept) candidate's regression outputs. Model,
      --  Forcing and D are already kept at accept time (KeepModel,
      --  KeepForcing, DKeep); these were not, so after the search loop
      --  the reported trend/accel/windings came from the LAST candidate
      --  evaluated, which may have been rejected. Restored after the loop.
      Keep_Trend, Keep_Accel, Keep_Aero : Long_Float := 0.0;
      --  Training and held-out scores of the kept candidate (SELECT_TOP)
      Keep_CC, Keep_P : Long_Float := 0.0;
      Own_Score, Keep_Own : Long_Float := 0.0;  -- this leg's own score
      --  LEG_SCORE FAST work space
      MAP_Use : Modulations_Amp_Phase (1 .. NM + NH);
      Model_Low : Data_Pairs := Data_Records;
      Model_Full : Data_Pairs := Data_Records;
      --  The one-sided lag mask of the full model, reused for the pass
      --  without the fast part so that the fast part is a LINEAR piece of
      --  the model. (With a mask of its own, that pass differed from the
      --  full model by lagged slow terms as well, and a search on the fast
      --  score exploited it by driving IR to -2 and -3.)
      Full_Side : Side_Mask (Data_Records'Range) := (others => True);
      Save_Model : Data_Pairs := Data_Records;
      Save_Data : Data_Pairs := Data_Records;
      Final_Reported : Boolean := False;
      Keep_M : Modulations (1 .. NM + NH);
      Keep_MAP : Modulations_Amp_Phase (1 .. NM + NH);
      Annual_Cycle : Annual_Harmonics;
      Keep_Initial_Value, Init_Value0 : Long_Float := 0.0;

      function Modulations_Are_Valid return Boolean is
      begin
         if Minimum_Modulation <= 0.0 then
            raise Constraint_Error with "MIN_LT must be positive";
         end if;
         for I in 1 .. NM loop
            if abs (D.B.LT (I)) < Minimum_Modulation then
               return False;
            end if;
         end loop;
         return True;
      end Modulations_Are_Valid;

      function Jerked_Tidal_Factors return Long_Periods_Amp_Phase is
         Result : Long_Periods_Amp_Phase := D.B.LPAP;
         Reference_Period : constant Long_Float := 13.660_830_77;
         Normal_Weight, Derivative_Weight : Long_Float;
      begin
         if Jerk < 0.0 or else Jerk > 1.0 then
            raise Constraint_Error with "JERK must be between 0.0 and 1.0";
         end if;
         if Jerk = 0.0 then
            return Result;
         end if;

         for I in Result'Range loop
            if D.A.LP (I) = 0.0 then
               raise Constraint_Error with "Cannot differentiate a zero-period tide";
            end if;
            Normal_Weight := 1.0 - Jerk;
            Derivative_Weight := Jerk * Reference_Period / D.A.LP (I);
            -- Tide_Sum uses cosine phases here; d/dt cos(theta) is
            -- proportional to cos(theta + pi / 2).
            Result (I).Amplitude :=
              Result (I).Amplitude *
                LEF.Sqrt (Normal_Weight**2 + Derivative_Weight**2);
            Result (I).Phase := Result (I).Phase +
              LEF.Arctan (Derivative_Weight, Normal_Weight);
         end loop;
         return Result;
      end Jerked_Tidal_Factors;

      --  mB (2026-10-10): second-order shaping of the impulse response. The
      --  impulses pass through a second real pole R2 before the integrator,
      --     G(i) = (1 - R2) * Raw(i) + R2 * G(i-1),
      --  so together with the integrator's own pole 1 - mA the response to an
      --  impulse is  (1 - R2) / ((1 - (1-mA) z^-1) (1 - R2 z^-1)):  the plateau
      --  is reached along an exponential of time constant -1 / ln (mB) months
      --  instead of in one step, and then declines at the mA / mP rate as
      --  before. Both poles are real, so it cannot ring; mB = 1 - mA is the
      --  critically damped case (a double pole). The factor 1 - R2 keeps the
      --  plateau height what it was. mB is per MONTH like mA and mP
      --  (R2 = mB ** (12 / SAMPLING)); |mB| is used and held below 1.
      --  mB = 0 (the default, and a value the search never leaves): no change.
      function Shape_Impulses (Raw : Data_Pairs) return Data_Pairs is
         use Ada.Numerics.Long_Elementary_Functions;
      begin
         if D.B.mB = 0.0 then
            return Raw;
         end if;
         declare
            B : constant Long_Float := Long_Float'Min (abs D.B.mB, 0.9999);
            R2 : constant Long_Float := (if Monthly then B else B ** Step_Scale);
            G : Data_Pairs := Raw;
         begin
            G (G'First).Value := (1.0 - R2) * Raw (Raw'First).Value;
            for I in G'First + 1 .. G'Last loop
               G (I).Value := (1.0 - R2) * Raw (I).Value + R2 * G (I - 1).Value;
            end loop;
            return G;
         end;
      end Shape_Impulses;

      --  IMPULSE_EXACT: a copy of the template in which every sample that
      --  carries an impulse has its date moved to the exact impulse date
      --  (at most half a sample away); Restore_Dates puts the sample dates
      --  back on the tide values computed from it.
      function Exact_Template (T : Data_Pairs) return Data_Pairs is
      begin
         if not Impulse_Exact or else GEM.Zonal.On or else not Impulse_Only then
            return T;
         end if;
         declare
            R : Data_Pairs := T;
            Tau : constant Long_Float :=
              abs D.B.DelB - Long_Float'Floor (abs D.B.DelB);
            F, D1, D2 : Long_Float;
         begin
            for I in T'Range loop
               if Impulse (T (I).Date) /= 0.0 then
                  F := T (I).Date - Long_Float'Floor (T (I).Date);
                  D1 := Tau - F;
                  D1 := D1 - Long_Float'Rounding (D1);
                  D2 := Tau + 0.5 - F;
                  D2 := D2 - Long_Float'Rounding (D2);
                  R (I).Date := T (I).Date + (if abs D1 <= abs D2 then D1 else D2);
               end if;
            end loop;
            return R;
         end;
      end Exact_Template;

      function Restore_Dates (R, T : Data_Pairs) return Data_Pairs is
         X : Data_Pairs := R;
      begin
         if Impulse_Exact then
            for I in X'Range loop
               X (I).Date := T (I).Date;
            end loop;
         end if;
         return X;
      end Restore_Dates;

      function Calc_Forcing return Data_Pairs is
         F : Data_Pairs := Forcing;
      begin
         der := 1.0 - D.B.mA; -- keeps the integrator stable

         if Strict_IDate then
            --  Extend the (pure-date) template back to IDATE before
            --  Tide_Sum/Impulse_Amplify/IIR ever run, instead of letting
            --  IIR's own Start_Index search clamp to Data_Records'First.
            --  Tide_Sum and Impulse_Amplify (with Offset=Ramp=0.0, as
            --  called here) are both pure per-element functions of .Date
            --  alone -- see their own bodies -- so extending the template
            --  needs no real data for the synthetic rows. The extended
            --  result's real-data portion sits at exactly
            --  Data_Records'Range (Extend_Backward only ever prepends),
            --  so slicing back is a direct index copy, not a
            --  date-matching search.
            --
            --  REVISED (Option A, see project notes): seeded at
            --  Initial_Conditions_Date (IDATE), NOT Init_Date. D.B.init
            --  now ALWAYS means "value at IDATE" -- the same convention
            --  the non-Strict branch below (and ws.py's own iir(),
            --  independently) already uses -- rather than "value at
            --  Init_Date", which required a human to keep `init`
            --  numerically re-synced by hand every time Init_Date was
            --  changed (found the hard way: changing a resp's INIT_DATE
            --  alone, without also updating `init`, either silently
            --  mis-anchors the whole trajectory or drives
            --  Regression_Factors singular -- deadlocking every thread,
            --  caught only because of the Ever_Accepted guard).
            --
            --  Init_Date is no longer a seeding input at all here: with
            --  init pinned to IDATE, the computed trajectory is now
            --  IDENTICAL regardless of what Init_Date says -- declaring a
            --  different Init_Date (e.g. for a companion file's own
            --  cross-file reuse) changes nothing about THIS run's own
            --  fit, only what a human would read off the resulting
            --  Forcing curve at that date to seed elsewhere. Reporting
            --  that value automatically (reading Forcing at Init_Date
            --  when it falls inside Data_Records'Range) is left for a
            --  follow-on diagnostic, not implemented here.
            declare
               Extended_Template : constant Data_Pairs :=
                 Extend_Backward
                   (Data_Records, Initial_Conditions_Date, Sampling_Per_Year);
               Extended_Impulses : constant Data_Pairs :=
                 Impulse_Amplify
                   (Raw => Restore_Dates (
                      Tide_Sum
                        (Template => Exact_Template (Extended_Template),
                         Constituents => Jerked_Tidal_Factors,
                         Periods => D.A.LP, Ref_Time => 0.0, Scaling => 0.0,
                         Year_Len => Year_Length (D.B.Year), Integ => D.B.ShiftT,
                         Lod => Lod_Term, Rho => D.B.rho,
                         ZTau => D.B.ztau, ZDist => D.B.zdist, ZTilt => D.B.ztilt,
                         SsaA => D.B.ssaA, SsaP => D.B.ssaP), Extended_Template),
                    Offset => 0.0, Ramp => 0.0,
                    Start => Extended_Template (Extended_Template'First).Date);
               --  Start_Index searches for the first Date STRICTLY
               --  GREATER than Start -- passing Initial_Conditions_Date
               --  itself would land one row LATE whenever it exactly
               --  equals a real sample's own date (the common case:
               --  IDATE defaults to Data_Records'First's own date),
               --  seeding `init` a full sample late and throwing off
               --  everything downstream (caught during development: this
               --  exact bug produced discrepancies up to ~0.29 against
               --  the non-strict path on the SAME file, which this
               --  seeding is supposed to reproduce exactly whenever
               --  Initial_Conditions_Date already matches
               --  Data_Records'First). A tenth of a sample step is
               --  comfortably inside the gap to the next real row for
               --  any reasonable Sampling_Per_Year, so this can't skip
               --  back an extra row.
               Extended_F : constant Data_Pairs :=
                 IIR
                   (Raw => Shape_Impulses (Extended_Impulses),
                    lagA => (if Monthly then der
                             else 1.0 - (1.0 - der) * Step_Scale),
                    lagC => (if Monthly then D.B.mP else D.B.mP * Step_Scale),
                    iA => D.B.init,
                    Start => Initial_Conditions_Date - 0.1 / Sampling_Per_Year);
            begin
               for I in Data_Records'Range loop
                  F (I).Value := Extended_F (I).Value;
               end loop;
            end;
         else
            Impulses :=
              Impulse_Amplify
                (Raw => Restore_Dates (
                   Tide_Sum
                     (Template => Exact_Template (Data_Records), Constituents => Jerked_Tidal_Factors,
                      Periods => D.A.LP, Ref_Time => 0.0, Scaling => 0.0,
                      Year_Len => Year_Length (D.B.Year), Integ => D.B.ShiftT,
                         Lod => Lod_Term, Rho => D.B.rho,
                         ZTau => D.B.ztau, ZDist => D.B.zdist, ZTilt => D.B.ztilt,
                         SsaA => D.B.ssaA, SsaP => D.B.ssaP), Data_Records),
                 Offset => 0.0, Ramp => 0.0,
                 Start => Data_Records (Data_Records'First).Date);

               F :=
                 IIR
                   (Raw => Shape_Impulses (Impulses),
                    lagA => (if Monthly then der
                             else 1.0 - (1.0 - der) * Step_Scale),
                    lagC => (if Monthly then D.B.mP else D.B.mP * Step_Scale),
                    iA => D.B.init, -- iB => 0.0, iC => 0.0,
                    Start => Initial_Conditions_Date);
         end if;

         return F;
      end Calc_Forcing;

   begin
      if Data_Records'Length = 0 then
         Text_IO.Put_Line (File_Name & " empty or not found");
         GNAT.OS_Lib.OS_Exit (0);
      end if;

      --  Populate the canonical winding candidate set (harmless to build
      --  even when Lock_Winding is False -- only ever read when it's
      --  True). Both signs of each backbone*harmonic, plus both signs of
      --  the backbone/subharmonic term.
      for I in Canon_Harmonics'Range loop
         Canon_Candidates (2 * I - 1) :=
           Canon_Backbone * Long_Float (Canon_Harmonics (I));
         Canon_Candidates (2 * I) := -Canon_Candidates (2 * I - 1);
      end loop;
      Canon_Candidates (2 * Canon_N - 1) :=
        Canon_Backbone / Long_Float (Canon_Subharmonic);
      Canon_Candidates (2 * Canon_N) := -Canon_Candidates (2 * Canon_N - 1);

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

         --  LEGS: the follower evaluates the candidate its leader posts, and
         --  judges it against the leader's own current best.
         if Follower then
            Exchanges (Pair).Fetch (Leg)
              (D, Harms, Pair_Threshold, Follower_Quit);
            exit when Follower_Quit;
            Old_CC := Pair_Threshold;
         end if;
         if Has_Leg_IR then
            D.B.IR := Leg_IR;
         end if;

         GEM.LTE.Year_Adjustment (D.B.Year, D.A.LP);
         if not Modulations_Are_Valid then
            raise Constraint_Error with "Loaded LT modulation is below MIN_LT";
         end if;

         -- Tidal constituents summed, amplified by impulse, and LTE modulated
         --  (a pair leader does this further down, after posting the
         --  candidate, so that both legs are evaluated at the same time)
         if not Leader then
            Forcing := Calc_Forcing;
         end if;

         M (1 .. NM) := D.B.LT (1 .. NM);
         if Lock_Freq then
            M (NM) := 1.0 / (Decay * D.B.mP);
         end if;
         MAP (1 .. NM) := D.A.LTAP (1 .. NM);
         for I in 1 .. NH loop
            M (NM + I) := Long_Float (Harms (I)) * M (NM);
         end loop;
         
         -- Check for duplicate periods that could cause singular regression
         -- matrix. Recover by regenerating the offending harmonic's own
         -- multiplier via Force_Harmonic — an UNCONDITIONAL redraw (unlike
         -- Random_Harmonic, which only mutates probabilistically, gated by
         -- FLIP/FIX; with the default FLIP=0.0 it would never actually
         -- change anything, which is exactly what let this collision go
         -- unresolved originally). Force_Harmonic only ever draws from
         -- Harmonic_Index'Range = 2 .. Harmonic_Range (gem-random_descent
         -- .adb), so it can never re-propose the trivial multiplier=1 case
         -- that would trivially equal M(NM) itself. If a clean value can't
         -- be found within a bounded number of tries (a persistent
         -- configuration problem, not just a random coincidence — should
         -- essentially never happen), fall through leaving M as-is:
         -- Regression_Factors' own Constraint_Error handler already
         -- recovers gracefully from a singular design matrix, and the
         -- existing Singular-triggered reset below then picks a fresh
         -- starting point — the same recovery path any other rejected
         -- trial already goes through, rather than a special case that
         -- kills this thread outright.
         declare
            Has_Duplicates : Boolean := False;
         begin
            for Attempt in 1 .. 25 loop
               Has_Duplicates := False;

               -- Check if harmonic periods match any base periods
               for I in 1 .. NH loop
                  for J in 1 .. NM loop
                     if abs (M (NM + I) - M (J)) < 1.0E-10 then
                        Has_Duplicates := True;
                        Walker.Force_Harmonic (Harms (I));
                        M (NM + I) := Long_Float (Harms (I)) * M (NM);
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
                           Has_Duplicates := True;
                           Walker.Force_Harmonic (Harms (I));
                           M (NM + I) := Long_Float (Harms (I)) * M (NM);
                           exit;
                        end if;
                     end loop;
                     exit when Has_Duplicates;
                  end loop;
               end if;

               exit when not Has_Duplicates;
            end loop;

            if Has_Duplicates then
               Text_IO.Put_Line
                 ("Warning:" & ID'Img &
                  " harmonic collision persisted after regeneration " &
                  "attempts, resetting");
            end if;
         end;

         if Leader then
            Exchanges (Pair).Propose (D, Harms, Old_CC);
            Forcing := Calc_Forcing;
         end if;

         if MLR_On or not (Forcing_Only or Is_Minimum_Entropy) then
            if Climate_Trend then
               Secular_Trend := 1.0;
            else
               Secular_Trend := 0.0;
            end if;
            DR := Data_Records;
            if D.B.IR /= 0.0 and then not Uncompensated then
               declare
                  Side : constant Side_Mask := IR_Mask (Data_Records);
               begin
                  for I in reverse Data_Records'First + SPY .. Data_Records'Last loop
                     if Side (I - SPY) then
                        DR (I).Value := Data_Records (I).Value + D.B.IR*Data_Records (I - SPY).Value;
                     end if;
                  end loop;
               end;
            end if;
            --  Uncompensated=True: DR stays equal to raw Data_Records --
            --  Regression_Factors folds the delay differential into its
            --  own basis columns instead (see IR parameter below).
            
--            DR := Annual_Add (DR, -1.0);
            if Dither_A /= 0.0 and then Dither_Stage = 0 then
               Forcing := Add_Dither (Forcing);
            end if;
            if Layer > 0 and then NM >= 2 then
               Forcing := Layered
                 (Bessel (Forcing, D.B.ImpA, D.B.ImpB,
                          (if Layer_K > 0.0 then Layer_K
                           elsif NM >= 3 then M (2)
                           elsif Lock_Freq then M (NM - 1) else M (NM)),
                          0.0, 0.0, 0.0),
                  D.B.Offset, D.B.bg, M (1),
                  W_Mix => (if Layer = 1 and then Layer_Mix then D.B.ImpC
                            else 0.0));
            elsif Has_Friction then
               Forcing := Frictional(Forcing, D.B.ImpA, D.B.ImpB, 0.0);
            elsif NM = 1 then
               Forcing := Bessel(Forcing, D.B.ImpA, D.B.ImpB, M(NM)*(1.0-D.B.BG), 0.0, D.B.Offset, D.B.bg);
            else
               if Lock_Freq then
                  Forcing := Bessel(Forcing, D.B.ImpA, D.B.ImpB, M(NM-1), M(NM), D.B.Offset, D.B.bg);
               else
                  Forcing := Bessel(Forcing, D.B.ImpA, D.B.ImpB, M(NM), M(NM-1), D.B.Offset, D.B.bg);
               end if;
            end if;
            if Dither_A /= 0.0 and then Dither_Stage /= 0 then
               Forcing := Add_Dither (Forcing);
            end if;
            if Alpha /= 0.0 then
               M_Min := Manifold_Min (Forcing);
            end if;
            declare
               --  LAYER: the lowest winding is not a regressor.
               Lo : constant Integer :=
                 (if Layer > 0 and then NM >= 3 and then Layer_K <= 0.0
                  then 3
                  elsif Layer > 0 and then NM >= 2 then 2 else 1);
               --  LEG_COMMON: Idx lists the windings this task regresses; the
               --  others (harmonics above the split, on a follower) take the
               --  leader's amplitudes. Without it Idx is simply Lo .. NM+NH.
               Use_Common : constant Boolean := Follower and then Common_Split > 0;
               Idx : array (1 .. NM + NH) of Integer := (others => 0);
               N_Reg : Natural := 0;
               Lead : Modulations_Amp_Phase (1 .. NM + NH) := MAP;
            begin
            if Use_Common then
               Exchanges (Pair).Get_Amps (Leg) (Lead);
            end if;
            for J in Lo .. NM + NH loop
               if Use_Common and then J > NM
                 and then Harms (J - NM) > Common_Split
               then
                  MAP (J) := (Amplitude => Leg_Scale * Lead (J).Amplitude,
                              Phase => Lead (J).Phase);
               else
                  N_Reg := N_Reg + 1;
                  Idx (N_Reg) := J;
               end if;
            end loop;
            --  With LEG_SCORE FAST the slow terms are regressed on the data as
            --  they are, NOT on the data less the leader's fast part: the
            --  slow model must not be able to react to the fast part, or a
            --  fast part of runaway amplitude (near-collinear windings on
            --  the leader's leg) is partly cancelled by it and then
            --  correlates with what is left by construction. Searches did
            --  exploit exactly that.
            if Use_Common and then N_Reg < NM + NH - Lo + 1
              and then not Fast_Score
            then
               --  Take the leader's fast part out of the regression target
               declare
                  NC : constant Integer := NM + NH - Lo + 1 - N_Reg;
                  CM : Modulations (1 .. NC);
                  CA : Modulations_Amp_Phase (1 .. NC);
                  K : Natural := 0;
               begin
                  for J in Lo .. NM + NH loop
                     if J > NM and then Harms (J - NM) > Common_Split then
                        K := K + 1;
                        CM (K) := M (J);
                        CA (K) := MAP (J);
                     end if;
                  end loop;
                  declare
                     HP : constant Data_Pairs :=
                       LTE (Forcing => Forcing, Wave_Numbers => CM,
                            Amp_Phase => CA, NonLin => NonLin,
                            Accel_Ref => Accel_Ref);
                  begin
                     for I in DR'Range loop
                        --  compensated: the target is the pre-lag series;
                        --  uncompensated: it is the data, i.e. after the lag
                        if Uncompensated and then I - SPY >= DR'First then
                           DR (I).Value := DR (I).Value -
                             (HP (I).Value - D.B.IR * HP (I - SPY).Value);
                        else
                           DR (I).Value := DR (I).Value - HP (I).Value;
                        end if;
                     end loop;
                  end;
               end;
            end if;
            declare
               RM : Modulations (1 .. N_Reg);
               RMAP : Modulations_Amp_Phase (1 .. N_Reg);
            begin
            for R in 1 .. N_Reg loop
               RM (R) := M (Idx (R));
               RMAP (R) := MAP (Idx (R));
            end loop;
            Regression_Factors
              (Data_Records =>
                 (if Coverage < 1.0 then Coverage_Region (DR)
                  elsif Enclosing then Lower_Region (DR)
                  else Excluded (DR)), -- Time series
               Forcing =>
                 (if Coverage < 1.0 then Coverage_Region (Forcing)
                  elsif Enclosing then Lower_Region (Forcing)
                  else Excluded (Forcing)),  -- Value @ Time
               NM => RM'Length, -- # modulations
               DBLT => RM, --D.B.LT,
               DALTAP => RMAP, --D.A.LTAP,
               DALEVEL => D.A.level,
               DAK0 => D.A.k0, 
               Secular_Trend => Secular_Trend, 
               Accel => Accel,
               Accel_Ref => Accel_Ref,
               Aero => Aero,
               Singular => Singular, 
               Annual => Annual_Cycle,
               Third => 0.0,
               IR => D.B.IR,
               NonLin => NonLin,
               Alpha => Alpha,
               M_Min => M_Min
            );
            for R in 1 .. N_Reg loop
               MAP (Idx (R)) := RMAP (R);
            end loop;
            end;
            --  Zero amplitude for the windings that are not regressors
            --  (layer, and Bessel with NM >= 3): LTE (given the full
            --  arrays) then adds nothing for them, and the windings JSON
            --  shows them as 0.
            for J in 1 .. Lo - 1 loop
               MAP (J) := (Amplitude => 0.0, Phase => 0.0);
            end loop;
            if Leader and then Common_Split > 0 then
               Exchanges (Pair).Post_Amps (MAP);
            end if;
            end;
            if Climate_Trend and then not Singular then
               D.B.Ann1 := Annual_Cycle.Ann1;
               D.B.Ann2 := Annual_Cycle.Ann2;
               D.B.Sem1 := Annual_Cycle.Semi1;
               D.B.Sem2 := Annual_Cycle.Semi2;
            end if;
         else
            Singular := False;
         end if;

         CorrCoeff := 0.0;
         --  One pass normally. With LEG_SCORE FAST pass 1 builds the full
         --  model and pass 2 the model without the leader's fast part.
         for Pass in 1 .. (if Fast_Score then 2 else 1) loop
            MAP_Use := MAP (1 .. NM + NH);
            if Pass = 2 then
               for J in NM + 1 .. NM + NH loop
                  if Harms (J - NM) > Common_Split then
                     MAP_Use (J) := (Amplitude => 0.0, Phase => 0.0);
                  end if;
               end loop;
            end if;
            if Forcing_Only or (Is_Minimum_Entropy and not MLR_On) then
               Model := Forcing;
            else
               Model :=
                 LTE
                   (Forcing => Forcing,
                    Wave_Numbers =>
                      M (1 .. NM+NH), --D.B.Lt(1..NM),
                    Amp_Phase => MAP_Use, --D.A.LTAP,
                    Offset => D.A.level, 
                    K0 => D.A.k0, 
                    Trend => Secular_Trend,
                    Accel => Accel,
                    NonLin => NonLin,
                    Annual => Annual_Cycle,
                    Third => 0.0,
                    --  Anchor Accel at the date Regression_Factors
                    --  itself reports centring it on. Forcing(First).Date
                    --  (used before) is only that date when EXCLUDE is
                    --  off; under EXCLUDE, First is the start of the
                    --  excluded test gap, which re-centred Accel on 2000
                    --  while the fit used 1950 -- a linear tilt of
                    --  100*Accel per year relative to the model's own MLR.
                    Accel_Ref => Accel_Ref,
                    Aero => Aero,
                    Alpha => Alpha,
                    M_Min => M_Min);

--               Model := Annual_Add (Model);
               if Dither_A /= 0.0 then
                  Model := Remove_Climatology (Model);
               end if;
               if D.B.ImpC /= 0.0
                 and then not (Layer = 1 and then Layer_Mix)
               then
                  for I in Model'Range loop
                     Model (I).Value :=
                       Model (I).Value + Annual_Impulse (Model (I).Date);
                  end loop;
               end if;

               -- Delay differential
               if D.B.IR /= 0.0 then
                  declare
                     Side : constant Side_Mask :=
                       (if Pass = 2 then Full_Side else IR_Mask (Model));
                  begin
                     if Fast_Score and then Pass = 1 then
                        Full_Side := Side;
                     end if;
                     for I in reverse Model'First + SPY .. Model'Last loop
                        if Side (I - SPY) then
                           Model (I).Value := Model (I).Value - D.B.IR*Model (I - SPY).Value;
                        end if;
                     end loop;
                  end;
               end if;


            end if;


            if Filter9Pt > 0 then
               for F in 1 .. Filter9Pt loop
                  Model := Filter9Point (Model);
               end loop;
            elsif Filter9Pt < 0 then
               --  F9 < 0: no filtering of any kind, on either side (Data
               --  is only filtered for F9 > 0, above). Added 2026-09-30 to
               --  expose sub-annual structure the default filters hide.
               null;
            else
               -- extra filtering, 2 equal-weighted 3-point box windows creating triangle
               Model := Median (Model);
               Model :=
                 FIR
                   (FIR (Model, Filter, 1.0 - 2.0 * Filter, Filter), Filter,
                    1.0 - 2.0 * Filter, Filter);
            end if;

            if Fast_Score and then Pass = 1 then
               Model_Full := Model;
            elsif Pass = 2 then
               Model_Low := Model;
               Model := Model_Full;
            end if;
         end loop;

            -- pragma Debug ( Dump(Model, Data_Records, Run_Time) );

            --  LEG_SCORE FAST: score the fast part against the data less
            --  the rest of the model (restored after the scoring below)
            if Fast_Score then
               Save_Model := Model;
               Save_Data := Data_Records;
               for I in Model'Range loop
                  Model (I).Value := Save_Model (I).Value - Model_Low (I).Value;
                  Data_Records (I).Value :=
                    Save_Data (I).Value - Model_Low (I).Value;
               end loop;
            end if;

            if Coverage < 1.0 then
               -- COVERAGE: fit was on Coverage_Region (the leading fraction
               -- above); score on the WHOLE training interval so a thread's
               -- own accept/reject (Old_CC/Keep, below) directly rejects
               -- parameter sets that don't predict the withheld remainder.
               CorrCoeff :=
                 Metric
                   (Model (First .. Last), Data_Records (First .. Last),
                    Forcing (First .. Last));
               CorrCoeffP := Exclude_Metric; -- true out-of-band, unaffected
            elsif Enclosing then
               -- ENCLOSING: fit was on LOWER (above); accept/reject is
               -- driven by ENCLOSED_UPPER_ONLY's choice of region(s) — see
               -- its declaration for the "U" vs "UL" tradeoff. CorrCoeffP
               -- reports the enclosed [First,Last] gap either way — the
               -- true held-out region in this mode, since both flanks are
               -- now spent on fitting/selection, not Exclude_Metric's
               -- LOWER+UPPER (no longer held out here at all).
               if Enclosed_Upper_Only then
                  CorrCoeff :=
                    Metric
                      (Model (Last .. Model'Last),
                       Data_Records (Last .. Data_Records'Last),
                       Forcing (Last .. Forcing'Last));
               else
                  CorrCoeff := Exclude_Metric; -- LOWER + UPPER combined
               end if;
               CorrCoeffP :=
                 Metric
                   (Model (First .. Last), Data_Records (First .. Last),
                    Forcing (First .. Last));
            elsif Split_Training then
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

         -- BUG FIX (was unconditional): CorrCoeffP is only assigned in the
         -- "else" branch above (Split_Training's branch sets CorrCoeff/
         -- CorrCoeffTest directly and never touches CorrCoeffP), so an
         -- unconditional "CorrCoeff := CorrCoeffP" here clobbered
         -- Split_Training's correctly-computed CorrCoeff with an
         -- uninitialized value every iteration — which also broke each
         -- thread's own accept/reject below (Old_CC never improved past 0),
         -- silently degrading the search into undirected random sampling
         -- whenever SPLIT_TRAINING=true. COVERAGE's and ENCLOSING's
         -- branches above already set both CorrCoeff and CorrCoeffP
         -- directly, so both are excluded here the same way Split_Training
         -- is — this guard is only for the plain EXCLUDE/no-modifier case.
         if Coverage >= 1.0 and then not Enclosing
           and then not Split_Training
         then
            CorrCoeff := CorrCoeffP;
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

         if Fast_Score then
            Model := Save_Model;
            Data_Records := Save_Data;
         end if;

         Own_Score := CorrCoeff;

         --  LEGS: the score from here on is the pair's composite.
         if Leader then
            declare
               C, P : Long_Float;
            begin
               Exchanges (Pair).Result (CorrCoeff, CorrCoeffP, Leg_Weight, C, P);
               CorrCoeff := C;
               CorrCoeffP := P;
            end;
         elsif Follower then
            Exchanges (Pair).Answer (Leg, CorrCoeff, CorrCoeffP, Leg_Weight);
            Exchanges (Pair).Verdict (Leg) (CorrCoeff);
         end if;

         --  BUG FIX (ordering): this accept/reject block used to run
         --  AFTER "register the results with a monitor" below, so on
         --  every iteration Live_Validate_Score was computed from
         --  KeepModel as it stood BEFORE this iteration's own update --
         --  one full iteration stale. Harmless most of the time (KeepModel
         --  changes slowly), but catastrophic on the very FIRST iteration
         --  of a run: KeepModel's own declaration-time initializer is
         --  `:= Data_Records` (see its declaration), so on iteration 1,
         --  before any accept has ever happened, Validate_Metric(KeepModel)
         --  was really comparing (approximately) Data against Data --
         --  confirmed live via DEBUG_CC (KeepModel's stats in the scored
         --  region matched Data's variance, not Model's). That spurious
         --  near-self-correlation (0.99+) then gets reported to Check on
         --  iteration 1, and because Best_Validate_Live is a running
         --  MAXIMUM, gets permanently locked in as "the best ever seen"
         --  before the search has done any real work -- exactly the
         --  "Validate stuck near 1.0 while training keeps improving"
         --  symptom. Moved before the monitor registration so KeepModel is
         --  always this iteration's own freshly-decided state.
         if CorrCoeff > Old_CC then
            Old_CC := CorrCoeff;
            Ever_Accepted := True;
            Keep := Set;
            KeepModel := Model;
            KeepForcing := Forcing;
            --  The harmonic multipliers live in the local Harms array;
            --  D.C is only read once at start-up. Without copying them
            --  here, DKeep carried the ORIGINAL file's harmonics, and the
            --  save-time rebuild of M from DKeep.C (2026-09-29) wrote
            --  those back to lt.exe.p -- so a search that accepted a
            --  harmonic change saved a file that did not reproduce its
            --  own kept model (reload scored lower than displayed).
            for I in 1 .. NH loop
               exit when I > D.C'Last;
               D.C (I) := Harms (I);
            end loop;
            DKeep := D;
            Keep_Trend := Secular_Trend;
            Keep_Accel := Accel;
            Keep_Aero := Aero;
            Keep_CC := CorrCoeff;
            Keep_Own := Own_Score;
            Keep_P := CorrCoeffP;
            Keep_M := M;
            Keep_MAP := MAP;
            if Catchup and then not Follower then  -- save it for other threads to reset from
               GEM.LTE.Primitives.Shared.Put (D);
            end if;
            Harms_Keep := Harms;
            Fingerprint ("accept", D, Model, Forcing, M, Harms,
                         Secular_Trend, Aero, CorrCoeff, ID);
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

         -- Register the results with a monitor
         declare
            --  Live, display-only validate readout of this thread's own
            --  best-kept candidate (KeepModel, not the live/just-tried
            --  Model) -- see Monitor.Check's own comment: never fed back
            --  into Best/Best_Client above, purely for Status's progress
            --  line. Skipped (kept at 0.0) when VALIDATE is off, to avoid
            --  the extra Metric call on every iteration for nothing.
            --
            --  Also 0.0, not Validate_Metric(KeepModel), whenever this
            --  thread hasn't accepted anything yet (Ever_Accepted False):
            --  KeepModel is still literally raw Data at that point (see
            --  the deadlock guard below), and Validate_Metric would be
            --  comparing Data against itself -- confirmed live, this
            --  leaked a spurious V:0.95 onto the display even after
            --  Report_Final's own copy of this same mistake was fixed.
            Live_Validate_Score : constant Long_Float :=
              (if Validate and then Ever_Accepted
               then Validate_Metric (KeepModel) else 0.0);
         begin
            if Follower then
               Best := False;
               Percentage := 100;
               Monitor.Leader (Best_Client);
            elsif Split_Training then
               Monitor.Check
                 (CorrCoeffTest, CorrCoeff, ID, Counter, Live_Validate_Score,
                  Best, Best_Client, Percentage);
            else
               Monitor.Check
                 (CorrCoeff, CorrCoeffP, ID, Counter, Live_Validate_Score,
                  Best, Best_Client, Percentage);
            end if;
         end;

         if ID = Best_Client then
            Counter := 1; -- no use penalizing thread in the lead
         end if;

         if not Follower and then
           (Singular or
            ((not Best) and Counter > Maximum_Loops and Percentage < Threshold))
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

         --  A follower stops only when its leader has (Fetch reports it)
         exit when Halted and then not Follower;

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
         elsif Follower then
            null;  -- the next candidate comes from the leader
         elsif Vary_Initial then
            Walker.Markov (D.B.init, Keep_Initial_Value, Spread, Init_Value0);
         else
            Sv_ZTau := D.B.ztau; Sv_ZDist := D.B.zdist; Sv_ZTilt := D.B.ztilt;
            Sv_SsaA := D.B.ssaA; Sv_SsaP := D.B.ssaP;
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
            if not Free_ZTau then D.B.ztau := Sv_ZTau; end if;
            if not Free_ZDist then D.B.zdist := Sv_ZDist; end if;
            if not Free_ZTilt then D.B.ztilt := Sv_ZTilt; end if;
            if not Free_Ssa then D.B.ssaA := Sv_SsaA; D.B.ssaP := Sv_SsaP; end if;
            if Lock_Winding then
               --  Markov perturbs exactly one Set index per call; this
               --  finds whichever one changed (if any) within the ltep
               --  slice and snaps it to the nearest canonical harmonic/
               --  subharmonic, so the search can only choose AMONG the
               --  canonical set for its winding value(s), never drift to
               --  an arbitrary alias -- amplitude/phase expression of
               --  whichever candidate gets picked is still fully free,
               --  via the ordinary optimization of every other
               --  parameter (impA/impB/delA/delB/asym/etc.).
               for I in LT_Lo .. LT_Hi loop
                  if Set (I) /= Keep (I) then
                     Set (I) := Walker.Nearest_In_Set (Set (I), Canon_Candidates);
                  end if;
               end loop;
            end if;
            Walker.Random_Harmonic (Harms, Harms_Keep);
            if not Modulations_Are_Valid then
               Set := Keep;
            end if;
         end if;

      end loop;
      Monitor.Stop;
      if Leader then
         Exchanges (Pair).Close;
      elsif Follower then
         Monitor.Leader (Best_Client);
      end if;
      if Ever_Accepted then
         --  Report the KEPT candidate's regression outputs, matching
         --  KeepModel/DKeep (see Keep_Trend's declaration).
         Secular_Trend := Keep_Trend;
         Accel := Keep_Accel;
         Aero := Keep_Aero;
         M := Keep_M;
         MAP := Keep_MAP;
      end if;

      --  DEADLOCK GUARD: KeepModel's/DKeep's own declaration-time
      --  initializers (`:= Data_Records` / `:= D`) are a plain
      --  convenience to give them the right bounds/shape before any real
      --  candidate exists -- never meant to be a standalone fallback
      --  result. If this thread's search never once accepted a
      --  candidate (CorrCoeff > Old_CC, starting at 0.0), KeepModel is
      --  still that exact, untouched copy of the real data -- silently
      --  saving it would write raw Data into the "Model" column,
      --  indistinguishable from a suspiciously perfect fit. Confirmed
      --  live: kN040_E130_1856's lte_results.csv had Model == Data ==
      --  the raw .dat file, bit for bit, traced to exactly this path.
      --  A thread that never accepts anything across its whole run is
      --  deadlocked, not merely unlucky -- most likely because
      --  Regression_Factors keeps hitting its own Singular fallback
      --  (e.g. from near-collinear, closely-spaced windings), so no
      --  perturbation ever scores above the Old_CC=0.0 starting floor.
      --
      --  IMPORTANT: this disqualifies only THIS thread's own candidate
      --  from ever being saved or winning the VALIDATE lockbox below --
      --  it does NOT halt the process. One thread's regression being
      --  chronically singular says nothing about the other threads'
      --  independent random-walk trajectories, which may be finding
      --  perfectly good candidates; killing the whole run over one
      --  thread's bad luck would throw those away too. (An earlier
      --  version of this guard called GNAT.OS_Lib.OS_Exit here,
      --  unconditionally halting every thread -- corrected per explicit
      --  feedback: "the fact that a single thread deadlocks doesn't mean
      --  that other threads can't continue.")
      if not Ever_Accepted then
         Text_IO.Put_Line
           ("WARNING: " & File_Name & " thread" & ID'Img &
            " never accepted a single candidate (deadlocked -- "
            & "Regression_Factors likely singular throughout this "
            & "thread's run). Excluding this thread's candidate from "
            & "saving/winning; other threads continue normally.");
      end if;

      --  VALIDATE lockbox: every thread reports its own best-kept
      --  candidate (DKeep/KeepModel, never the live/possibly-just-
      --  rejected D/Model) and that candidate's score on the held-out
      --  validation segment; only the last thread to report is told to
      --  save -- using whichever candidate (potentially from a different
      --  thread) scored best on validation. Default VALIDATE=False keeps
      --  today's exact Best_Client-based decision, untouched.
      declare
         Save_Now : Boolean;
      begin
         if Validate and not Test_Only and not Paired then
            declare
               Am_I_Last : Boolean;
            begin
               Monitor.Report_Final
                 (D => DKeep, Model => KeepModel, Forcing => KeepForcing,
                  M => M, MAP => MAP,
                  Trend => Secular_Trend, Accel => Accel, Aero => Aero,
                  --  A deadlocked thread's own KeepModel is untouched raw
                  --  Data (see the guard above) -- Validate_Metric on
                  --  Data-vs-itself can score deceptively HIGH, exactly
                  --  the failure this whole investigation started from.
                  --  Report the same "can never win" sentinel Winner
                  --  itself starts from, so this thread's candidate is
                  --  never selected unless literally every thread in the
                  --  run is equally deadlocked (in which case there is no
                  --  real candidate to prefer anyway).
                  Validate_Score =>
                    (if Ever_Accepted then Validate_Metric (KeepModel)
                     else Long_Float'First),
                  Am_I_Last => Am_I_Last);
               if Am_I_Last then
                  Monitor.Winner
                    (D => DKeep, Model => KeepModel, Forcing => KeepForcing,
                     M => M, MAP => MAP,
                     Trend => Secular_Trend, Accel => Accel, Aero => Aero);
                  --  Keep D.A/D.B AND the live Model in sync w/ the
                  --  winning DKeep/KeepModel: the reporting below (CorrCoeff/
                  --  CorrCoeffP, Exclude_Metric) reads Model directly, so
                  --  without this it would print numbers from this thread's
                  --  OWN last-tried candidate instead of the one actually
                  --  being saved. KeepForcing is now part of the same
                  --  cross-thread swap (2026-09-29 fix) -- Winner may hand
                  --  back a DIFFERENT thread's D/Model/M/MAP than this
                  --  thread's own, and Forcing must come from that SAME
                  --  winning thread's own accept-time snapshot, not this
                  --  thread's local KeepForcing left over from before the
                  --  swap.
                  D := DKeep;
                  Model := KeepModel;
               end if;
               Save_Now := Am_I_Last;
            end;
         else
            --  Same disqualification for the VALIDATE=FALSE/legacy path:
            --  a deadlocked thread's own KeepModel is untouched raw Data,
            --  so it must never be the one saved even if (in some
            --  Ratio-weighted edge case) it still ended up as Best_Client.
            --  With legs the choice is always made centrally (SELECT_TOP 0
            --  acts as 1 = best training score), so that the leader and the
            --  follower of one and the same pair save: each task deciding
            --  for itself could leave the two leg files from different pairs.
            if (Select_Top > 0 or else Paired) and then not Test_Only then
               if not Follower then
                  Final_Reported := True;
                  Monitor.Final_Report
                    (ID, (if Ever_Accepted then Keep_CC else Long_Float'First),
                     Keep_P);
               end if;
               Monitor.Selected (Leader_ID, Save_Now);
               Save_Now := Save_Now and then Ever_Accepted;
            else
               Save_Now := Ever_Accepted and then
                 (Test_Only or Best_Client = Leader_ID);
            end if;
            if Ever_Accepted then
               --  As in the VALIDATE branch above: report and save the
               --  KEPT candidate, not the live last-tried one (the final
               --  score line below is computed from Model).
               D := DKeep;
               Model := KeepModel;
            end if;
         end if;

      --  Forcing is recomputed unconditionally every iteration regardless
      --  of accept/reject, so by the time Save fires (possibly many
      --  rejected perturbations after this thread's -- or, under
      --  VALIDATE=TRUE, even a DIFFERENT thread's -- own last accept) the
      --  live Forcing can silently disagree with the KeepModel/DKeep
      --  actually being saved, in BOTH the Validate and legacy paths
      --  above (neither resyncs Forcing itself). Save()'s Forcing column
      --  (lte_results.csv's column 4, the "manifold") could then
      --  disagree with its own Model/DKeep/lt.exe.p -- found 2026-09-29
      --  while cross-checking experiments/Feb2026/lte_forward.py's
      --  reimplementation column-by-column: a from-scratch Python replay
      --  of a saved lt.exe.p matched a live, single-threaded, first-
      --  iteration Ada evaluation of that same lt.exe.p to 1e-13, yet
      --  disagreed with that same run's own saved lte_results.csv
      --  Forcing column by ~1e-4 relative -- exactly this gap. KeepForcing
      --  is cached at the same accept-time assignment as KeepModel/DKeep
      --  (see "Keep := Set; KeepModel := Model;" above), so it is always
      --  valid whenever Ever_Accepted is True.
      if Ever_Accepted then
         Forcing := KeepForcing;
      end if;

      if Save_Now then
         Save_Mutex.Seize (ID);

         --  M (windings/frequencies) is recomputed unconditionally every
         --  iteration regardless of accept/reject (same root cause as the
         --  Forcing fix above), so by the time Save fires it can differ
         --  from DKeep's own accepted windings -- and the report loop
         --  below writes DKeep.C back FROM M
         --  (`DKeep.C(I-NM) := Integer(M(I)/M(NM))`), which
         --  Shared.Save(DKeep) then serializes as lt.exe.p's "harm"
         --  field. Unlike MAP/Save_Windings (amplitude/phase, write-only
         --  into windings.json, never read back by anything), lt.exe.p
         --  IS reloaded by future runs and by lte_forward.py's own
         --  harm-array fix, so this one specific derivation needs to be
         --  self-consistent with DKeep, not the live thread state.
         --  Rebuild M from DKeep itself first so that write-back becomes
         --  a harmless round-trip (2026-09-29).
         M (1 .. NM) := DKeep.B.LT (1 .. NM);
         if Lock_Freq then
            M (NM) := 1.0 / (Decay * DKeep.B.mP);
         end if;
         for I in 1 .. NH loop
            M (NM + I) := Long_Float (DKeep.C (I)) * M (NM);
         end loop;

         -- Text_IO.Put_Line("### " & File_Name);
         -- Walker.Dump(Keep); -- Print results of last best evaluation,
         GEM.LTE.Primitives.Shared.Dump (DKeep);

         Text_IO.Put_Line ("---- LTE ----");
         Put (Ridge_Lambda, " :ridge:", NL);
         Put (Coverage, " :coverage:", NL);
         Text_IO.Put_Line (" " & Boolean'Image (Enclosing) & " :enclosing:");
         Text_IO.Put_Line
           (" " & (if Enclosed_Upper_Only then "U" else "UL") &
            " :enclosed:");
         Put (Secular_Trend, " :trend:", NL);
         Put (Accel, " :accel:", NL);
         if GEM.LTE.Primitives.Aero_On then
            Put (Aero, " :aero:", NL);
         end if;
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

         declare
            Fit_Array : constant Data_Pairs := Fit_Region (Data_Records);
         begin
            GEM.LTE.Primitives.Shared.Save_Windings
              (Trend => Secular_Trend, Accel => Accel, K0 => D.A.k0,
               Level => D.A.level, IR => D.B.IR, M => M, MAP => MAP,
               NM => NM, NH => NH, B => D.B,
               Ctx =>
                 (Known          => True,
                  Accel_Ref      => Accel_Ref,
                  Fit_Start      => Fit_Array (Fit_Array'First).Date,
                  Fit_End        => Fit_Array (Fit_Array'Last).Date,
                  Interval_Start => Forcing (First).Date,
                  Interval_End   => Forcing (Last).Date,
                  Record_Start   => Data_Records (Data_Records'First).Date,
                  Record_End     => Data_Records (Data_Records'Last).Date,
                  Exclude        => Exclude,
                  Enclosing      => Enclosing,
                  Coverage       => Coverage,
                  Accel_Enabled  => GEM.Getenv ("ACCEL", True),
                  Aero_Enabled   => GEM.LTE.Primitives.Aero_On,
                  Aero_Coef      => Aero,
                  Aero_Lag       => GEM.LTE.Primitives.Aero_Lag,
                  Alpha          => Alpha,
                  M_Min          => (if Alpha /= 0.0 then Manifold_Min (Forcing) else 0.0)));
         end;

         Fingerprint ("save  ", DKeep, KeepModel, Forcing, M, Harms_Keep,
                      Secular_Trend, Aero, Old_CC, ID);
         if not Follower then  -- one parameter file, written by the leader
            GEM.LTE.Primitives.Shared.Save (DKeep);
         end if;
         Save (KeepModel, Data_Records, Forcing, File_Name => Results_File,
               IR => D.B.IR);    -- saves to file

         if Validate then
            --  VALIDATE=TRUE: print exactly the same three numbers the
            --  live Status display last showed (Best_Metric/Best_OOB/
            --  Best_Validate_Live, via the same non-blocking Monitor.
            --  Current peek Status itself reads) -- i.e. a copy of the
            --  running display, for the leader -- rather than freshly
            --  recomputing train/test from Model/KeepModel with mode-
            --  specific swap logic that can (and did) drift out of sync
            --  with what was actually shown throughout the run.
            declare
               Final_Train, Final_Test, Final_Validate : Long_Float;
            begin
               Monitor.Current (Final_Train, Final_Test, Final_Validate);
               Put_CC (Final_Train, Final_Test, Counter, ID, Final_Validate);
            end;
         elsif Split_Training then
            Put_CC (CorrCoeff, CorrCoeffTest, Counter, ID);
         elsif Enclosing then
            -- Mirror the per-iteration ENCLOSING branch above: report
            -- whichever region(s) drove accept/reject (per ENCLOSED_
            -- UPPER_ONLY) and the enclosed gap (the true held-out region
            -- in this mode).
            if Enclosed_Upper_Only then
               CorrCoeff :=
                 Metric
                   (Model (Last .. Model'Last),
                    Data_Records (Last .. Data_Records'Last),
                    Forcing (Last .. Forcing'Last));
            else
               CorrCoeff := Exclude_Metric;
            end if;
            CorrCoeffP :=
              Metric
                (Model (First .. Last), Data_Records (First .. Last),
                 Forcing (First .. Last));
            Put_CC (CorrCoeff, CorrCoeffP, Counter, ID);
         else
            --  VALIDATE=FALSE (default): original, long-standing,
            --  Exclude-oblivious definition, unchanged from before this
            --  session -- CorrCoeffP is UNCONDITIONALLY the concatenation
            --  AROUND the excluded [First,Last] area (Exclude_Metric),
            --  CorrCoeff is UNCONDITIONALLY the [First,Last] window
            --  itself, regardless of the Exclude flag.
            CorrCoeffP := Exclude_Metric;
            CorrCoeff :=
              Metric
                (Model (First .. Last), Data_Records (First .. Last),
                 Forcing (First .. Last));
            Put_CC (CorrCoeff, CorrCoeffP, Counter, ID);
         end if;

         if Paired then
            --  the score the pair was accepted on (both legs combined)
            Text_IO.Put ("COMPOSITE");
            Ada.Long_Float_Text_IO.Put (Old_CC, Fore => 4, Aft => 10, Exp => 0);
            Text_IO.Put ("  leg" & Leg'Img & " pair" & Pair'Img & "  own" &
                         (if Fast_Score then " (fast part)" else "") & ":");
            Ada.Long_Float_Text_IO.Put (Keep_Own, Fore => 3, Aft => 6, Exp => 0);
            Text_IO.New_Line;
         end if;

         CorrCoeff :=
           CompareRef (DKeep.A.LP, DKeep.B.LPAP, DKeep.B.Year, DKeep.B.rho,
                       DKeep.B.ztau, DKeep.B.zdist, DKeep.B.ztilt,
                       DKeep.B.ssaA, DKeep.B.ssaP);
         Put (CorrCoeff, ":dLOD:   ");
         Put (Year_Length (DKeep.B.Year), ":Yr: " & File_Name);
         Save_Mutex.Release (ID);

      else
         null; -- Text_IO.Put_Line("Exited " & ID'Img);
      end if;
      end;
      if Test_Only then
         Monitor.Test_Done;
      end if;

   exception
      when E : others =>
         Save_Mutex.Release (ID);
         if Leader then
            Exchanges (Pair).Close;
         elsif Follower then
            Exchanges (Pair).Follower_Gone (Leg);
         end if;
         if (Select_Top > 0 or else Paired) and then not Test_Only
           and then not Follower and then not Final_Reported
         then
            Monitor.Final_Report (ID, Long_Float'First, 0.0);
         end if;
         if Test_Only then
            Monitor.Test_Done;
         end if;
         Text_IO.Put_Line
           ("Solution err: " & Ada.Exceptions.Exception_Information (E));
         -- The following may need a debug-specifi compiler switch to activate
         Text_IO.Put_Line (GNAT.Traceback.Symbolic.Symbolic_Traceback (E));

         --  VALIDATE safety net: a thread dying here never reaches its
         --  normal Report_Final call above, which would otherwise leave
         --  Threads_Reported permanently short of Validate_Thread_Count --
         --  hanging every surviving thread's Am_I_Last check forever.
         --  Report a deliberately unbeatable-low score so this thread can
         --  never become Winner; the nested handler guards against a
         --  second failure (e.g. DKeep/KeepModel never having been fully
         --  assigned before the exception) turning a hang into a crash.
         if Validate then
            begin
               declare
                  Dummy_Last : Boolean;
               begin
                  Monitor.Report_Final
                    (D => DKeep, Model => KeepModel, Forcing => KeepForcing,
                     M => M, MAP => MAP,
                     Trend => Secular_Trend, Accel => Accel, Aero => Aero,
                     Validate_Score => Long_Float'First,
                     Am_I_Last => Dummy_Last);
               end;
            exception
               when others =>
                  null;
            end;
         end if;

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
