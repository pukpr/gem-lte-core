--  ENSO Optimization Driver
--
--  Main procedure for optimizing Laplace's Tidal Equation (LTE) model parameters
--  to fit climate indices (e.g., ENSO/Niño data) using parallel search tasks.
--
--  The program:
--  1. Loads initial parameters and dLOD (day length) reference data
--  2. Spawns parallel worker tasks to search parameter space
--  3. Monitors progress and handles user interrupts (q/s/x keys, timeout)
--  4. Saves best-fit parameters when optimization completes

with System.Task_Info;
with GEM.LTE.Primitives.Solution;
with GEM.LTE.Primitives.Refine;
with GEM.LTE.Primitives.Shared;
with Text_IO;
with GEM.dLOD;
with GEM.Tidal_Poly;
with GNAT.OS_Lib;
with Ada.Calendar;
with Ada.Directories;

procedure ENSO_Opt is  -- gprbuild lte.gpr enso_opt -largs -Wl, --stack=40000000

   --  Configuration from environment variables
   N : Positive :=  --  Number of parallel worker threads
     GEM.Getenv
       ("NUMBER_OF_PROCESSORS", System.Task_Info.Number_Of_Processors);
   dLOD_dat : String := GEM.Getenv ("DLOD_DAT", "../dlod3.dat");
   Ch : Character := ' ';  --  User keyboard input
   --  A test run takes no keyboard commands: reading the keyboard blocks
   --  when standard input is an idle pipe, which kept the run from ending.
   Test_Mode : constant Boolean := GEM.Getenv ("TEST_ONLY", False);
   Avail : Boolean;  --  Keyboard input available flag
   T : Ada.Calendar.Time;  --  Timeout tracking
   Cycle : Duration :=  --  Maximum optimization duration
     Duration (GEM.Getenv ("TIMEOUT", Long_Float (Duration'Last) / 2.0));
   dLOD_Scale : constant Long_Float := GEM.Getenv ("DLOD_SCALE", 0.0);
   Expect : constant Boolean := GEM.Getenv ("EXPECT", False);
   use type Ada.Calendar.Time;

   --  DLOD_REF: TRUE/REF (default), FALSE/CAL, or the name of a
   --  wisdom_lpap/<name>.lt.exe.p to seed lpap from; case-insensitive.
   Dlod_Ref_Value : constant String := GEM.Getenv ("DLOD_REF", "TRUE");
   type Dlod_Ref_Mode is (Ref, Cal, Named);
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
   Dlod_Up : constant String := Upper (Dlod_Ref_Value);
   Dlod_Mode : constant Dlod_Ref_Mode :=
     (if Dlod_Up = "TRUE" or else Dlod_Up = "T" or else Dlod_Up = "REF"
        or else Dlod_Up = "YES" or else Dlod_Up = "1" or else Dlod_Up = ""
      then Ref
      elsif Dlod_Up = "FALSE" or else Dlod_Up = "F" or else Dlod_Up = "CAL"
        or else Dlod_Up = "NO" or else Dlod_Up = "0"
      then Cal
      else Named);

   --  D = Complete model parameter set (tidal periods, amplitudes, phases,
   --  LTE modulation parameters, impulse response coefficients)
   D : GEM.LTE.Primitives.Shared.Param_S :=
     (NLP => GEM.LTE.LP'Length, NLT => GEM.LTE.LTM'Length,
      A =>
        (k0 => 0.0, level => 0.0, NLP => GEM.LTE.LP'Length,
         NLT => GEM.LTE.LTM'Length, LP => GEM.LTE.LP, LTAP => GEM.LTE.LTAP),
      B =>
        (NLP => GEM.LTE.LP'Length, NLT => GEM.LTE.LTM'Length,
         LPAP => GEM.LTE.LPAP, LT => GEM.LTE.LTM, Offset => 0.0, bg => 0.0,
         ImpA => 1.0, ImpB => 0.0, ImpC => 0.0, DelA => 0.0, DelB => 0.0,
         Asym => 0.0, Ann1 => 0.0, Ann2 => 0.0, Sem1 => 0.0, Sem2 => 0.0,
         Year => 0.0, IR => 0.0, mA => 0.0, mP => 0.0, shiftT => 0.000_00,
         init => 0.006_3, rho => 1.0, ztau => 1.0, zdist => 1.0, ztilt => 0.0,
         ssaA => 0.0, ssaP => 0.0, mB => 0.0, bg2 => 0.0, delC => 0.0,
         besK => 0.0),
      C => (others => 0), CW1 => 0.0, CW2 => 0.0);

begin
   --  Stack size limit is now raised by the launcher (ulimit -s unlimited)
   --  before this executable is invoked, so the setrlimit block below is
   --  no longer needed here.
   --
   --  --  Remove the main thread stack size limit (equivalent to ulimit -s unlimited).
   --  --  The default ~8 MB limit causes STORAGE_ERROR with large data sets.
   --  declare
   --     RLIMIT_STACK : constant Interfaces.C.int := 3;
   --     type rlim_t is mod 2 ** 64;
   --     type Rlimit is record
   --        rlim_cur, rlim_max : rlim_t;
   --     end record;
   --     pragma Convention (C, Rlimit);
   --     function setrlimit
   --       (resource : Interfaces.C.int; rlp : access Rlimit)
   --        return Interfaces.C.int;
   --     pragma Import (C, setrlimit, "setrlimit");
   --     R : aliased Rlimit := (rlim_t'Last, rlim_t'Last);
   --     use type Interfaces.C.int;
   --  begin
   --     if setrlimit (RLIMIT_STACK, R'Access) /= 0 then
   --        Text_IO.Put_Line ("Warning: could not remove stack size limit");
   --     end if;
   --  end;

   declare
      --  AP = Amplitude/Phase data from day-length-of-day (dLOD) measurements
      --  Used as reference forcing for the tidal model
      AP : GEM.LTE.Long_Periods_Amp_Phase (GEM.LTE.Doodson_Args'Range);
   begin

      Text_IO.Put_Line (N'Img & " processors available, timeout=" & Cycle'Img);

      --  For TIDES=POLY, set the polynomial basis periods BEFORE Shared.Load
      --  so saved .p files match their polynomial triplets by period.
      if GEM.Tidal_Poly.On then
         for I in D.A.LP'Range loop
            if I <= GEM.Tidal_Poly.Active_Count then
               D.A.LP (I) := GEM.Tidal_Poly.Term_Period_Days (I);
            else
               D.A.LP (I) := 0.0;
            end if;
         end loop;
      end if;

      --  Load previously saved parameters if available (warm start)
      GEM.LTE.Primitives.Shared.Load (D);

      --  Set the startup candidate periods before calibrating dLOD phases.
      if not GEM.Tidal_Poly.On then
         GEM.LTE.Year_Adjustment (D.B.Year, D.A.LP);
      end if;
      AP := GEM.dLOD (dLOD_dat, D.B.Year);

      Text_IO.Put_Line ("YA=" & D.B.Year'Img);

      --  Initialize tidal forcing amplitudes and phases from dLOD data.
      --  DLOD_REF (case-insensitive):
      --    TRUE / REF (default) : lpap as loaded; dLOD kept as reference
      --    FALSE / CAL          : lpap set to the dLOD calibration
      --    <name>               : lpap seeded from <name>.lt.exe.p in the
      --                           nearest experiments/wisdom_lpap (its lpap
      --                           rows only); dLOD kept as reference
      if GEM.Command_Line_Option_Exists ("r") or Dlod_Mode = Cal then
         --  Fresh start: use dLOD data directly with optional scaling
         Text_IO.Put_Line ("Loading dLOD");
         for I in D.B.LPAP'Range loop
            GEM.LTE.LPRef (I).Amplitude := AP (I).Amplitude;
            GEM.LTE.LPRef (I).Phase := AP (I).Phase;
            D.B.LPAP (I).Amplitude :=
              AP (I).Amplitude * (1.0 + dLOD_Scale * D.A.LP (I));
            D.B.LPAP (I).Phase := AP (I).Phase;
         end loop;
      else
         --  Reference mode: store dLOD but don't overwrite current parameters
         Text_IO.Put_Line ("Referencing dLOD");
         for I in D.B.LPAP'Range loop
            GEM.LTE.LPRef (I).Amplitude := AP (I).Amplitude;
            GEM.LTE.LPRef (I).Phase := AP (I).Phase;
         end loop;
         if Dlod_Mode = Named then
            declare
               File : constant String :=
                 GEM.LTE.Primitives.Shared.Wisdom_Named (Dlod_Ref_Value);
               Seeded, Skipped : Natural;
            begin
               if File = "" then
                  Text_IO.Put_Line
                    ("DLOD_REF " & Dlod_Ref_Value & ": not TRUE/REF, FALSE/CAL " &
                     "or a file in wisdom_lpap; found " &
                     GEM.LTE.Primitives.Shared.Wisdom_Names);
                  GNAT.OS_Lib.OS_Exit (1);
               end if;
               GEM.LTE.Primitives.Shared.Seed_LPAP (D, File, Seeded, Skipped);
               Text_IO.Put_Line
                 ("DLOD_REF " & Dlod_Ref_Value & ": seeding lpap from " & File &
                  ":" & Seeded'Img & " lines seeded," & Skipped'Img &
                  " rows left as loaded");
            end;
         end if;
      end if;

      --  WISDOM_REF TRUE: seed lpap from the crowd-median tidal set
      --  (experiments/wisdom_lpap/wisdom_raw.lt.exe.p, or WISDOM_P) -- the
      --  counterpart of DLOD_REF FALSE above, and applied after it when
      --  both are set. LPRef, the dLOD reference, is left as it is.
      if GEM.Getenv ("WISDOM_REF", False) then
         declare
            File : constant String := GEM.LTE.Primitives.Shared.Wisdom_File;
            Seeded, Skipped : Natural;
         begin
            if File = "" or else not Ada.Directories.Exists (File) then
               Text_IO.Put_Line
                 ("WISDOM_REF: no wisdom file (" & File &
                  "); set WISDOM_P to a wisdom_raw.lt.exe.p");
               GNAT.OS_Lib.OS_Exit (1);
            end if;
            GEM.LTE.Primitives.Shared.Seed_LPAP (D, File, Seeded, Skipped);
            Text_IO.Put_Line
              ("Seeding lpap from wisdom " & File & ":" & Seeded'Img &
               " lines seeded," & Skipped'Img & " rows left as loaded");
         end;
      end if;

      --  REFINE=TRUE: run deterministic in-process GN/homotopy refinement and exit
      if GEM.Getenv ("REFINE", False) then
         GEM.LTE.Primitives.Refine.Run (D);
         return;
      end if;

      --  Store initialized parameters for worker tasks to access
      GEM.LTE.Primitives.Shared.Put (D);

   end;

   T := Ada.Calendar.Clock;

   --  Launch N parallel worker tasks to search parameter space
   --  Each task performs random descent optimization independently
   GEM.LTE.Primitives.Solution.Start (D.NLP, D.NLT, N);

   --  Main monitoring loop: check status, handle user input, enforce timeout
   for I in 1 .. Integer'Last loop

      --  Get status from worker tasks (blocking call)
      declare
         S : String := GEM.LTE.Primitives.Solution.Status;
      begin
         --  Periodic progress reporting
         if I mod GEM.LTE.Primitives.Solution.Check_Every_N_Loops = 0 then
            Text_IO.Put_Line (S & " #" & I'Img);
         end if;
         exit when GEM.LTE.Primitives.Halted;
      end;

      --  Check for user keyboard commands
      if Test_Mode then
         Avail := False;
      elsif Expect then
         Text_IO.Get (Ch);  --  Blocking for scripted use
         Avail := True;
      else
         begin
            Text_IO.Get_Immediate (Ch, Avail);  --  Non-blocking for interactive
         exception
            when Text_IO.End_Error =>
               Avail := False;  --  no keyboard (input redirected): carry on
         end;
      end if;

      if Avail then
         if Ch = 'q' or Ch = 's' then
            --  'q' or 's': Stop optimization and save results
            GEM.LTE.Primitives.Stop;
         elsif Ch = 'x' then
            --  'x': Exit immediately without saving
            Text_IO.Put_Line ("Exiting, no save");
            GNAT.OS_Lib.OS_Exit (0);
         elsif Ch in '1' .. '9' then
            --  '1'-'9': Adjust correlation threshold trigger dynamically
            GEM.LTE.Primitives.Solution.Set_Trigger
              (Character'Pos (Ch) - Character'Pos ('1') + 1);
         end if;
      end if;

      if Ch = 'q' then
         exit when GEM.LTE.Primitives.Halted;
      end if;

      --  Enforce timeout limit
      if Ada.Calendar.Clock > T + Cycle then
         GEM.LTE.Primitives.Stop;
         T := Ada.Calendar.Clock;
         exit when GEM.LTE.Primitives.Halted;
      end if;

   end loop;

   --  COMMENTED CODE: Debug message for main thread exit
   --  TODO: Can remove - leftover from debugging parallel task shutdown
   --Text_IO.Put_Line("Main exiting, flushing other tasks");

   --  Wait for worker tasks to complete their final output writes
   delay 5.0;

end ENSO_Opt;
