--  ============================================================================
--  GEM.LTE.Primitives.Shared - Parameter Storage and I/O Management
--  ============================================================================
--
--  PURPOSE:
--    Manages the shared parameter state for multi-threaded GEM-LTE optimization.
--    Provides thread-safe storage, persistence (read/write), and serialization
--    (text .par files and JSON) of model parameters.
--
--  KEY RESPONSIBILITIES:
--    1. Protected Server: Thread-safe parameter sharing across worker threads
--    2. Parameter Persistence: Save/load from .par files (text format)
--    3. JSON Serialization: Load parameters from JSON format
--    4. Parameter Validation: Read with format checking and error handling
--
--  DATA STRUCTURES:
--    - Param_S: Main parameter record (defined in .ads)
--      - B: Variable parameters (offsets, impulses, LTE coefficients)
--      - A: Constant parameters (tidal periods, amplitudes, phases)
--      - C: Harmonic configuration array
--    - Protected Server: Manages single shared instance with mutex access
--
--  FILE FORMATS:
--    - .par files: Human-readable text format (4-char tags + values)
--      Example: "offs  0.12345", "impA  1.23456", "ltep  2.71828"
--    - .json files: JSON format for easier configuration management
--      Example: {"offs": 0.12345, "impA": 1.23456, "ltep": [2.71828, ...]}
--    - .parms files: Binary Direct_IO format (legacy, fast but not portable)
--
--  NAMING CONVENTION:
--    Files are named after the executable:
--    - enso_opt.exe.par (general parameters)
--    - enso_opt.exe.nino4.dat.par (climate-index-specific parameters)
--    - enso_opt.exe.json (JSON format)
--
--  THREAD SAFETY:
--    Server protected object ensures only one thread can modify shared
--    parameters at a time. Get() blocks until Put() has been called at
--    least once.
--
--  ============================================================================

with Ada.Command_Line;
with Ada.Text_IO;
with GNAT.OS_Lib;
with GNATCOLL.JSON;
with Ada.Strings.Unbounded;
with Ada.Characters.Latin_1;
with Ada.IO_Exceptions;
with Ada.Directories;

package body GEM.LTE.Primitives.Shared is

   --  Access type required for protected object storage
   --  (cannot directly store discriminated record in protected object)
   type Param_P is access all Param_S;

   --  =========================================================================
   --  Server: Protected object for thread-safe parameter sharing
   --
   --  PATTERN:
   --    Single-writer-multiple-reader with blocking Get(). First Put()
   --    allocates storage and sets Available flag, subsequent Put() calls
   --    reuse the allocated space. Get() blocks until Available is true.
   --
   --  USAGE:
   --    - Put() called by optimization thread when better parameters found
   --    - Get() called by worker threads to initialize or update their state
   --    - One-time allocation on first Put() avoids repeated heap allocation
   --  =========================================================================
   protected Server is
      procedure Put (P : in Param_S);
      entry Get (P : in out Param_S);
   private
      Available : Boolean := False;
      Params : Param_P;
   end Server;

   ------------
   -- Server --
   ------------

   protected body Server is

      --  Store parameters (allocates on first call, updates thereafter)
      procedure Put (P : in Param_S) is
      begin
         if not Available then -- only do once to allocate
            Params := new Param_S'(P);
         else  -- use allocated space
            Params.all := P;
         end if;
         Available := True;
      end Put;

      --  Retrieve parameters (blocks until first Put() completes)
      entry Get (P : in out Param_S) when Available is
      begin
         P := Params.all;
      end Get;

   end Server;

   --  Public wrappers for Server access
   procedure Put (P : in Param_S) is
   begin
      Server.Put (P);
   end Put;

   function Get (N_Tides, N_Modulations : in Integer) return Param_S is
      P : Param_S (N_Tides, N_Modulations);
   begin
      Server.Get (P);
      return P;
   end Get;


   --  =========================================================================
   --  Write_JSON: Save parameters to JSON .p files
   --
   --  Generates two files:
   --    1. <executable>.p - General parameters for all climate indices
   --    2. <executable>.<climate_index>.p - Index-specific parameters
   --
   --  Format: JSON with same structure as Read_JSON expects
   --  =========================================================================
   
   function File_To_String (Name : in String) return String;  -- Forward declaration
   
   --  Detect NM and NH from JSON file array lengths
   --  Overrides default values to match actual JSON array sizes
   procedure Detect_NM_NH_From_JSON (Filename : in String) is
      use GNATCOLL.JSON;
      Result : Read_Result;
      Data : JSON_Value;
      Arr : JSON_Array;
      NM_Count, NH_Count : Natural;
   begin
      -- Read JSON file
      declare
         Text : constant String := File_To_String (Filename);
      begin
         Result := Read (Text);
      exception
         when Ada.Text_IO.Name_Error | Ada.IO_Exceptions.Name_Error =>
            return;  -- File doesn't exist, skip detection
         when others =>
            return;
      end;
      
      if not Result.Success then
         return;  -- Parse error, skip detection
      end if;
      
      Data := Result.Value;
      
      -- Count ltep array length and override NM
      if Kind (Data) = JSON_Object_Type and then Has_Field (Data, "ltep") then
         Arr := Get (Data, "ltep");
         NM_Count := Length (Arr);
         GEM.Setenv ("NM", Integer'Image (NM_Count));
         --  Silently detect NM from JSON (uncomment for debug)
         --  Ada.Text_IO.Put_Line ("Detected NM=" & Integer'Image (NM_Count) & " from JSON ltep array");
      end if;
      
      -- Count harm array length and override NH
      if Kind (Data) = JSON_Object_Type and then Has_Field (Data, "harm") then
         Arr := Get (Data, "harm");
         NH_Count := Length (Arr);
         -- Build NH string as "1 1 1 ..." with NH_Count entries
         declare
            NH_Str : String (1 .. NH_Count * 2 - 1) := (others => ' ');
            Pos : Positive := 1;
         begin
            for I in 1 .. NH_Count loop
               NH_Str (Pos) := '1';
               if I < NH_Count then
                  Pos := Pos + 2;  -- Skip space
               end if;
            end loop;
            GEM.Setenv ("NH", NH_Str);
            --  Silently detect NH from JSON (uncomment for debug)
            --  Ada.Text_IO.Put_Line ("Detected NH with" & Natural'Image (NH_Count) & " entries from JSON harm array");
         end;
      end if;
   exception
      when others =>
         null;  -- Silently ignore any errors during detection
   end Detect_NM_NH_From_JSON;
   
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

   procedure Write_JSON (D : in Param_S) is
      use GNATCOLL.JSON;
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & ".p";

      -- Use NH/NM from environment (.resp file takes precedence)
      NM : constant Integer := GEM.Getenv ("NM", D.B.LT'Length);
      Harms : constant Ns := Parse_NH (GEM.Getenv ("NH", ""));
      NH : constant Integer := Harms'Length;
      
      procedure Write_To_File (Filename : String; Include_LPAP : Boolean);
      
      procedure Write_To_File (Filename : String; Include_LPAP : Boolean) is
         Obj : constant JSON_Value := Create_Object;
         LPAP_Arr : JSON_Array := Empty_Array;
         LT_Arr : JSON_Array := Empty_Array;
         Harm_Arr : JSON_Array := Empty_Array;
      begin
         -- Scalar parameters
         Set_Field (Obj, "offs", Create (D.B.Offset));
         Set_Field (Obj, "bg", Create (D.B.bg));
         Set_Field (Obj, "impA", Create (D.B.ImpA));
         Set_Field (Obj, "impB", Create (D.B.ImpB));
         Set_Field (Obj, "impC", Create (D.B.ImpC));
         Set_Field (Obj, "delA", Create (D.B.DelA));
         Set_Field (Obj, "delB", Create (D.B.DelB));
         Set_Field (Obj, "asym", Create (D.B.Asym));
         Set_Field (Obj, "ann1", Create (D.B.Ann1));
         Set_Field (Obj, "ann2", Create (D.B.Ann2));
         Set_Field (Obj, "sem1", Create (D.B.Sem1));
         Set_Field (Obj, "sem2", Create (D.B.Sem2));
         Set_Field (Obj, "cw1", Create (D.CW1));
         Set_Field (Obj, "cw2", Create (D.CW2));
         Set_Field (Obj, "year", Create (D.B.Year));
         Set_Field (Obj, "IR", Create (D.B.IR));
         Set_Field (Obj, "ma", Create (D.B.mA));
         Set_Field (Obj, "mp", Create (D.B.mP));
         Set_Field (Obj, "shfT", Create (D.B.shiftT));
         Set_Field (Obj, "init", Create (D.B.init));
         Set_Field (Obj, "rho", Create (D.B.rho));
         Set_Field (Obj, "ztau", Create (D.B.ztau));
         Set_Field (Obj, "zdist", Create (D.B.zdist));
         Set_Field (Obj, "ztilt", Create (D.B.ztilt));
         Set_Field (Obj, "ssaA", Create (D.B.ssaA));
         Set_Field (Obj, "ssaP", Create (D.B.ssaP));
         Set_Field (Obj, "mb", Create (D.B.mB));
         Set_Field (Obj, "bg2", Create (D.B.bg2));
         Set_Field (Obj, "delC", Create (D.B.delC));
         Set_Field (Obj, "besk", Create (D.B.besK));
         
         -- LPAP array (tidal constituents) - only in primary file
         if Include_LPAP then
            for I in D.B.LPAP'Range loop
               declare
                  Triplet : JSON_Array := Empty_Array;
               begin
                  Append (Triplet, Create (D.A.LP (I)));
                  Append (Triplet, Create (D.B.LPAP (I).Amplitude));
                  Append (Triplet, Create (D.B.LPAP (I).Phase));
                  Append (LPAP_Arr, Create (Triplet));
               end;
            end loop;
            Set_Field (Obj, "lpap", Create (LPAP_Arr));
         end if;
         
         -- ltep array (LT modulation periods) - only first NM terms
         for I in 1 .. NM loop
            Append (LT_Arr, Create (D.B.LT (I)));
         end loop;
         Set_Field (Obj, "ltep", Create (LT_Arr));
         
         -- harm array (harmonics) - only first NH values if NH specified
         if NH > 0 then
            for I in 1 .. NH loop
               Append (Harm_Arr, Create (Long_Float (D.C (I))));
            end loop;
         else
            -- Fallback: only non-zero values (original behavior)
            for I in D.C'Range loop
               exit when D.C (I) = 0;
               Append (Harm_Arr, Create (Long_Float (D.C (I))));
            end loop;
         end if;
         Set_Field (Obj, "harm", Create (Harm_Arr));
         
         -- Write to file
         declare
            JSON_Text : constant String := Write (Obj, Compact => False);
            FT : Ada.Text_IO.File_Type;
         begin
            Ada.Text_IO.Create (FT, Ada.Text_IO.Out_File, Filename);
            Ada.Text_IO.Put_Line (FT, JSON_Text);
            Ada.Text_IO.Close (FT);
         end;
      end Write_To_File;
      
   begin
      --  The per-CLIMATE_INDEX secondary file (FN2, "<exe>.<index>.p") used
      --  to also be written here -- removed alongside the matching read in
      --  Read (D) above: it duplicated every field except lpap, and its
      --  read side silently overwrote the primary's scalars/ltep/harm with
      --  whatever it happened to hold, which could go stale relative to
      --  lt.exe.p (confirmed live on "pdo", 2026-09-29). lt.exe.p alone is
      --  sufficient -- see that comment for the full explanation.
      Write_To_File (FN, Include_LPAP => True);
   end Write_JSON;

   --  =========================================================================
   --  Save_Windings: Persist the per-winding amplitude/phase table as its
   --  own JSON file, named "<exe>.windings.json". This mirrors the
   --  "---- LTE ----" stdout block emitted at the end of a run so the same
   --  amplitude/phase data can be collected from disk instead of scraped
   --  from stdout, for cross-index/cross-region phase analysis.
   --  Additive only: no other Save/Load/Dump behavior is touched.
   --
   --  Includes a "manifold" sub-object of the forcing-shaping parameters
   --  (DelA/DelB/Asym/mA/mP/ShiftT/ImpA/ImpB/Offset/Bg/Year). This is NOT
   --  needed to use the winding table on its own -- it exists so a
   --  cross-index/cross-region phase comparison can first check whether
   --  two files' manifolds actually agree closely enough to compare
   --  phases directly. Verified this session (two separate cases:
   --  Baltic MSL vs Kaplan SST, and kN_baltic vs kN_northsea) that two
   --  independently-optimized manifolds correlating as high as 99.6%
   --  overall can still decorrelate to near-zero or even NEGATIVE once
   --  multiplied by a winding number M inside sin(2*pi*M*Forcing) --
   --  manufacturing a spurious "antiphase" reading with no counterpart
   --  in the real data. Comparing this block first (mismatched values
   --  are a warning to rebuild a shared reference forcing rather than
   --  trust a direct phase comparison) is the fix.
   --  =========================================================================
   procedure Save_Windings
     (Trend, Accel, K0, Level, IR : in Long_Float;
      M                           : in Modulations;
      MAP                         : in Modulations_Amp_Phase;
      NM, NH                      : in Integer;
      B                           : in Param_B;
      Ctx                         : in Secular_Context := No_Secular_Context)
   is
      use GNATCOLL.JSON;
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) &
        ".windings" & GEM.Leg_Tag & ".json";
      Obj      : constant JSON_Value := Create_Object;
      K_Arr    : JSON_Array          := Empty_Array;
      Manifold : constant JSON_Value := Create_Object;

      function Span (A, B : Long_Float) return JSON_Array is
         Arr : JSON_Array := Empty_Array;
      begin
         Append (Arr, Create (A));
         Append (Arr, Create (B));
         return Arr;
      end Span;

      --  Local warming rate of the secular term as LTE evaluates it.
      function Rate (T : Long_Float) return Long_Float is
        (Trend + 2.0 * Accel * (T - Ctx.Accel_Ref));

      procedure Add_Secular_Context is
         Sec    : constant JSON_Value := Create_Object;
         Low    : constant JSON_Value := Create_Object;
         Lowest : Integer             := 1;
      begin
         Set_Field (Sec, "form",
           "model += level + k0*F(t) + trend*t"
           & " + accel*(t - accel_ref)**2"
           & " + aero*A(t - aero_lag)"
           & " + G*sum amp*sin(2*pi*k*F(t) + phase), G = ALPHA_GAIN if"
           & " F(t)-m_min <= alpha*|m_min| and the sum > 0, else 1"
           & "; t = decimal year, F = manifold (lte_results.csv col 4)");
         Set_Field (Sec, "accel_enabled", Create (Ctx.Accel_Enabled));
         Set_Field (Sec, "aero_enabled", Create (Ctx.Aero_Enabled));
         Set_Field (Sec, "aero_file", GEM.Getenv ("AERO", ""));
         Set_Field (Sec, "aero_lag", Create (Ctx.Aero_Lag));
         Set_Field (Sec, "alpha", Create (Ctx.Alpha));
         Set_Field (Sec, "m_min", Create (Ctx.M_Min));
         --  ALPHA_SLOPE > 0: the window is tested on the manifold's signed
         --  slope over that many months (negated when alpha_sign > 0), and
         --  m_min is the minimum of that quantity, not of F.
         Set_Field (Sec, "tides", GEM.Getenv ("TIDES", "TABLE"));
         Set_Field (Sec, "annual_dither",
                    Create (Long_Float'(GEM.Getenv ("ANNUAL_DITHER", 0.0))));
         Set_Field (Sec, "annual_dither_phase",
                    Create (Long_Float'(GEM.Getenv ("ANNUAL_DITHER_PHASE", 0.0))));
         Set_Field (Sec, "annual_dither_stage",
                    Create (Integer'(GEM.Getenv ("ANNUAL_DITHER_STAGE", 1))));
         Set_Field (Sec, "ir_side",
                    Create (Integer'(GEM.Getenv ("IR_SIDE", 0))));
         Set_Field (Sec, "ir_base_months",
                    Create (Integer'(GEM.Getenv ("IR_BASE", 121))));
         Set_Field (Sec, "alpha_slope_months",
                    Create (Integer'(GEM.Getenv ("ALPHA_SLOPE", 0))));
         Set_Field (Sec, "alpha_sign",
                    Create (Long_Float'(GEM.Getenv ("ALPHA_SIGN", -1.0))));
         Set_Field (Sec, "alpha_gain",
                    Create (Long_Float'(GEM.Getenv ("ALPHA_GAIN", 2.0))));
         Set_Field (Sec, "accel_ref", Create (Ctx.Accel_Ref));
         Set_Field (Sec, "accel_ref_note",
           "first date of the array the MLR was fitted on; the same date"
           & " LTE evaluates accel around. Warming rate at year t is"
           & " trend + 2*accel*(t - accel_ref).");
         Set_Field (Sec, "fit_span",
           Create (Span (Ctx.Fit_Start, Ctx.Fit_End)));
         Set_Field (Sec, "train_start_end",
           Create (Span (Ctx.Interval_Start, Ctx.Interval_End)));
         Set_Field (Sec, "train_start_end_role",
           (if Ctx.Coverage < 1.0 then "coverage: fit on leading fraction"
            elsif Ctx.Enclosing then "enclosing: interval is held out"
            elsif Ctx.Exclude
            then "excluded TEST interval (EXCLUDE=TRUE; training is"
                 & " the record outside it)"
            else "training interval"));
         Set_Field (Sec, "exclude", Create (Ctx.Exclude));
         Set_Field (Sec, "enclosing", Create (Ctx.Enclosing));
         Set_Field (Sec, "coverage", Create (Ctx.Coverage));
         Set_Field (Sec, "record_span",
           Create (Span (Ctx.Record_Start, Ctx.Record_End)));
         Set_Field (Sec, "rate_at_record_start",
           Create (Rate (Ctx.Record_Start)));
         Set_Field (Sec, "rate_at_record_end",
           Create (Rate (Ctx.Record_End)));
         Set_Field (Sec, "rise_over_record",
           Create (Trend * (Ctx.Record_End - Ctx.Record_Start)
                   + Accel * ((Ctx.Record_End - Ctx.Accel_Ref)**2
                              - (Ctx.Record_Start - Ctx.Accel_Ref)**2)));

         --  The lowest winding is a slow sinusoid in F(t); over a short
         --  record it bends like a quadratic and the regression can
         --  trade curvature between it and accel.
         for I in 2 .. NM + NH loop
            if abs M (I) < abs M (Lowest) then
               Lowest := I;
            end if;
         end loop;
         Set_Field (Low, "index", Create (Lowest));
         Set_Field (Low, "k", Create (M (Lowest)));
         Set_Field (Low, "amp", Create (MAP (Lowest).Amplitude));
         Set_Field (Low, "phase", Create (MAP (Lowest).Phase));
         Set_Field (Sec, "lowest_winding", Low);
         Set_Field (Obj, "secular_context", Sec);
      end Add_Secular_Context;
   begin
      Set_Field (Obj, "trend", Create (Trend));
      Set_Field (Obj, "accel", Create (Accel));
      Set_Field (Obj, "k0", Create (K0));
      Set_Field (Obj, "level", Create (Level));
      Set_Field (Obj, "IR", Create (IR));

      Set_Field (Manifold, "delA", Create (B.DelA));
      Set_Field (Manifold, "delB", Create (B.DelB));
      Set_Field (Manifold, "asym", Create (B.Asym));
      Set_Field (Manifold, "ma", Create (B.mA));
      Set_Field (Manifold, "mp", Create (B.mP));
      Set_Field (Manifold, "shiftT", Create (B.shiftT));
      Set_Field (Manifold, "rho", Create (B.rho));
      Set_Field (Manifold, "ztau", Create (B.ztau));
      Set_Field (Manifold, "zdist", Create (B.zdist));
      Set_Field (Manifold, "ztilt", Create (B.ztilt));
      Set_Field (Manifold, "ssaA", Create (B.ssaA));
      Set_Field (Manifold, "ssaP", Create (B.ssaP));
      Set_Field (Manifold, "mb", Create (B.mB));
      Set_Field (Manifold, "bg2", Create (B.bg2));
      Set_Field (Manifold, "delC", Create (B.delC));
      Set_Field (Manifold, "besk", Create (B.besK));
      Set_Field (Manifold, "impA", Create (B.ImpA));
      Set_Field (Manifold, "impB", Create (B.ImpB));
      Set_Field (Manifold, "offs", Create (B.Offset));
      Set_Field (Manifold, "bg", Create (B.bg));
      Set_Field (Manifold, "year", Create (B.Year));
      Set_Field (Obj, "manifold", Manifold);

      for I in 1 .. NM + NH loop
         declare
            Triplet : JSON_Array := Empty_Array;
         begin
            Append (Triplet, Create (M (I)));
            Append (Triplet, Create (MAP (I).Amplitude));
            Append (Triplet, Create (MAP (I).Phase));
            Append (K_Arr, Create (Triplet));
         end;
      end loop;
      Set_Field (Obj, "k_amp_phase", Create (K_Arr));
      if Ctx.Known then
         Set_Field (Obj, "aero", Create (Ctx.Aero_Coef));
      end if;
      if Ctx.Known and then Ctx.Tess_Freq /= 0.0 then
         declare
            Tess : constant JSON_Value := Create_Object;
         begin
            Set_Field (Tess, "freq", Create (Ctx.Tess_Freq));
            Set_Field (Tess, "sin", Create (Ctx.Tess_Sin));
            Set_Field (Tess, "cos", Create (Ctx.Tess_Cos));
            Set_Field (Obj, "tess", Tess);
         end;
      end if;
      if Ctx.Known and then NM + NH >= 1 then
         Add_Secular_Context;
      end if;

      declare
         JSON_Text : constant String := Write (Obj, Compact => False);
         FT        : Ada.Text_IO.File_Type;
      begin
         Ada.Text_IO.Create (FT, Ada.Text_IO.Out_File, FN);
         Ada.Text_IO.Put_Line (FT, JSON_Text);
         Ada.Text_IO.Close (FT);
      end;
   exception
      when others =>
         -- Never let this best-effort side output disturb a real run.
         null;
   end Save_Windings;

   --  =========================================================================
   --  Save: Persist parameters to binary .parms file, text .par, and JSON .p files
   --
   --  Uses Ada.Direct_IO for fast binary serialization (legacy format).
   --  Also writes human-readable .par files and JSON .p files to keep all formats in sync.
   --  =========================================================================
   procedure Save (P : in Param_S) is
   begin
      Write_JSON (P);  -- Write JSON .p files (keep in sync with .par)
   end Save;

   --  =========================================================================
   --  Input Parsing Helpers - Read parameters from text files
   --  =========================================================================

   function File_To_String (Name : in String) return String is
      FT : Ada.Text_IO.File_Type;
      Line : String (1 .. 1_024);
      Last : Natural;
      Contents : Ada.Strings.Unbounded.Unbounded_String :=
        Ada.Strings.Unbounded.Null_Unbounded_String;
   begin
      Ada.Text_IO.Open (FT, Ada.Text_IO.In_File, Name);
      begin
         while not Ada.Text_IO.End_Of_File (FT) loop
            Ada.Text_IO.Get_Line (FT, Line, Last);
            Ada.Strings.Unbounded.Append (Contents, Line (1 .. Last));
            Ada.Strings.Unbounded.Append (Contents, Ada.Characters.Latin_1.LF);
         end loop;
         Ada.Text_IO.Close (FT);
      exception
         when others =>
            if Ada.Text_IO.Is_Open (FT) then
               Ada.Text_IO.Close (FT);
            end if;
            raise;
      end;
      return Ada.Strings.Unbounded.To_String (Contents);
   end File_To_String;

   function Get_Number (Value : in GNATCOLL.JSON.JSON_Value) return Long_Float
   is
      use GNATCOLL.JSON;
   begin
      case Kind (Value) is
         when JSON_Float_Type =>
            return Get_Long_Float (Value);
         when JSON_Int_Type =>
            declare
               Int_Value : constant Long_Long_Integer := Get (Value);
            begin
               return Long_Float (Int_Value);
            end;
         when others =>
            Ada.Text_IO.Put_Line ("JSON number error:" & Kind (Value)'Img);
            GNAT.OS_Lib.OS_Exit (0);
      end case;
   end Get_Number;

   function Get_Integer (Value : in GNATCOLL.JSON.JSON_Value) return Integer is
      use GNATCOLL.JSON;
   begin
      case Kind (Value) is
         when JSON_Int_Type =>
            declare
               Int_Value : constant Long_Long_Integer := Get (Value);
            begin
               return Integer (Int_Value);
            end;
         when JSON_Float_Type =>
            return Integer (Get_Long_Float (Value));
         when others =>
            Ada.Text_IO.Put_Line ("JSON integer error:" & Kind (Value)'Img);
            GNAT.OS_Lib.OS_Exit (0);
      end case;
   end Get_Integer;

   function Get_Float_Field
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String;
      Default : in Long_Float) return Long_Float
   is
      use GNATCOLL.JSON;
   begin
      if Kind (Data) = JSON_Object_Type and then Has_Field (Data, Name) then
         return Get_Number (Get (Data, Name));
      end if;
      return Default;
   end Get_Float_Field;

   procedure Read_JSON_Float_Array
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String;
      Target : in out Modulations)
   is
      use GNATCOLL.JSON;
      Arr : JSON_Array;
      Count : Natural;
   begin
      if Kind (Data) /= JSON_Object_Type or else not Has_Field (Data, Name)
      then
         return;
      end if;
      Arr := Get (Data, Name);
      Count := Length (Arr);
      declare
         J : Positive := 1;
      begin
         for I in Target'Range loop
            exit when J > Count;
            Target (I) := Get_Number (Get (Arr, J));
            J := J + 1;
         end loop;
      end;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON array error " & Name);
         GNAT.OS_Lib.OS_Exit (0);
   end Read_JSON_Float_Array;

   procedure Read_JSON_Int_Array
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String; Target : in out Ns)
   is
      use GNATCOLL.JSON;
      Arr : JSON_Array;
      Count : Natural;
   begin
      if Kind (Data) /= JSON_Object_Type or else not Has_Field (Data, Name)
      then
         return;
      end if;
      Arr := Get (Data, Name);
      Count := Length (Arr);
      declare
         J : Positive := 1;
      begin
         for I in Target'Range loop
            exit when J > Count;
            Target (I) := Get_Integer (Get (Arr, J));
            J := J + 1;
         end loop;
      end;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON array error " & Name);
         GNAT.OS_Lib.OS_Exit (0);
   end Read_JSON_Int_Array;

   --  =========================================================================
   --  JSON Reading Functions - Modern configuration format
   --
   --  Provides flexible JSON-based parameter loading with fallback to
   --  text .par format. JSON format is easier to edit and validate.
   --  =========================================================================

   --  Read JSON array of Long_Float triplets (period, amplitude, phase)
   --  Used for LPAP (Long Period Amplitude Phase) tidal constituent data
   procedure Read_JSON_LPAP
     (Data : in GNATCOLL.JSON.JSON_Value; Name : in String; D : in out Param_S)
   is
      use GNATCOLL.JSON;
      Arr : JSON_Array;
      Count : Natural;
      Matches_Found : Natural := 0;
      Matched : array (D.B.LPAP'Range) of Boolean := (others => False);
      
      procedure Apply (Period, Amp, Phase : in Long_Float) is
         Found : Boolean := False;
         Best_Index : Integer := 0;
         Best_Diff : Long_Float := Long_Float'Last;
      begin
         -- Find the closest unmatched period
         for I in D.B.LPAP'Range loop
            if not Matched (I) and then abs (D.A.LP (I)) > 0.0 then
               declare
                  Diff : constant Long_Float := abs (abs (Period) - abs (D.A.LP (I)));
                  Tol : constant Long_Float := abs (D.A.LP (I)) * 0.01;
               begin
                  if Diff <= Tol and then Diff < Best_Diff then
                     Best_Diff := Diff;
                     Best_Index := I;
                     Found := True;
                  end if;
               end;
            end if;
         end loop;
         
         if Found and Best_Index > 0 then
            Ada.Text_IO.Put_Line 
              ("  MATCH: JSON period" & Long_Float'Image (Period) &
               " -> Index" & Integer'Image (Best_Index) &
               " (LP=" & Long_Float'Image (D.A.LP (Best_Index)) &
               ", diff=" & Long_Float'Image (Best_Diff) &
               ", Amp=" & Long_Float'Image (Amp) &
               ", Phase=" & Long_Float'Image (Phase) & ")");
            D.B.LPAP (Best_Index).Amplitude := Amp;
            D.B.LPAP (Best_Index).Phase := Phase;
            Matched (Best_Index) := True;
            Matches_Found := Matches_Found + 1;
         else
            Ada.Text_IO.Put_Line 
              ("  NO MATCH: JSON period" & Long_Float'Image (Period) &
               " (Amp=" & Long_Float'Image (Amp) &
               ", Phase=" & Long_Float'Image (Phase) & ")");
         end if;
      end Apply;
   begin
      Ada.Text_IO.Put_Line ("=== Reading JSON LPAP from: " & Name & " ===");
      if Kind (Data) /= JSON_Object_Type or else not Has_Field (Data, Name)
      then
         Ada.Text_IO.Put_Line ("  Field '" & Name & "' not found in JSON");
         return;
      end if;
      Arr := Get (Data, Name);
      Count := Length (Arr);
      Ada.Text_IO.Put_Line ("  Array contains" & Natural'Image (Count) & " triplets");
      if Count = 0 then
         return;
      end if;
      if Kind (Get (Arr, 1)) = JSON_Array_Type then
         for I in 1 .. Count loop
            declare
               Triplet : constant JSON_Array := Get (Get (Arr, I));
            begin
               if Length (Triplet) >= 3 then
                  Apply
                    (Get_Number (Get (Triplet, 1)),
                     Get_Number (Get (Triplet, 2)),
                     Get_Number (Get (Triplet, 3)));
               end if;
            end;
         end loop;
      else
         declare
            Index : Positive := 1;
         begin
            while Index + 2 <= Count loop
               Apply
                 (Get_Number (Get (Arr, Index)),
                  Get_Number (Get (Arr, Index + 1)),
                  Get_Number (Get (Arr, Index + 2)));
               Index := Index + 3;
            end loop;
         end;
      end if;
      Ada.Text_IO.Put_Line 
        ("=== LPAP Summary: Matched" & Natural'Image (Matches_Found) &
         " out of" & Natural'Image (Count) & " JSON triplets ===" & ASCII.LF);
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON array error " & Name);
         GNAT.OS_Lib.OS_Exit (0);
   end Read_JSON_LPAP;

   --  Read complete parameter set from JSON file
   --  Returns True on success, False if file not found or parse error
   --  Include_LPAP flag controls whether tidal constituent data is loaded
   function Read_JSON
     (Name : in String; D : in out Param_S; Include_LPAP : in Boolean)
      return Boolean
   is
      use GNATCOLL.JSON;
      Result : Read_Result;
      Data : JSON_Value;
   begin
      declare
         Text : constant String := File_To_String (Name);
      begin
         Result := Read (Text);
      exception
         when Ada.Text_IO.Name_Error | Ada.IO_Exceptions.Name_Error =>
            return False;
         when others =>
            return False;
      end;
      if not Result.Success then
         Ada.Text_IO.Put_Line (Format_Parsing_Error (Result.Error));
         GNAT.OS_Lib.OS_Exit (0);
      end if;
      Data := Result.Value;
      D.B.Offset := Get_Float_Field (Data, "offs", D.B.Offset);
      D.B.bg := Get_Float_Field (Data, "bg", D.B.bg);
      D.B.ImpA := Get_Float_Field (Data, "impA", D.B.ImpA);
      D.B.ImpB := Get_Float_Field (Data, "impB", D.B.ImpB);
      D.B.ImpC := Get_Float_Field (Data, "impC", D.B.ImpC);
      D.B.DelA := Get_Float_Field (Data, "delA", D.B.DelA);
      D.B.DelB := Get_Float_Field (Data, "delB", D.B.DelB);
      D.B.Asym := Get_Float_Field (Data, "asym", D.B.Asym);
      D.B.Ann1 := Get_Float_Field (Data, "ann1", D.B.Ann1);
      D.B.Ann2 := Get_Float_Field (Data, "ann2", D.B.Ann2);
      D.B.Sem1 := Get_Float_Field (Data, "sem1", D.B.Sem1);
      D.B.Sem2 := Get_Float_Field (Data, "sem2", D.B.Sem2);
      --  CW alias pair, absent in older files: 0 = off
      D.CW1 := Get_Float_Field (Data, "cw1", 0.0);
      D.CW2 := Get_Float_Field (Data, "cw2", 0.0);
      D.B.Year := Get_Float_Field (Data, "year", D.B.Year);
      D.B.IR :=
        Get_Float_Field (Data, "IR", Get_Float_Field (Data, "ir", D.B.IR));
      D.B.mA := Get_Float_Field (Data, "ma", D.B.mA);
      D.B.mP := Get_Float_Field (Data, "mp", D.B.mP);
      D.B.shiftT :=
        Get_Float_Field
          (Data, "shfT", Get_Float_Field (Data, "shft", D.B.shiftT));
      D.B.init := Get_Float_Field (Data, "init", D.B.init);
      --  absent in older files: 1.0 = pure Meeus generator
      D.B.rho := Get_Float_Field (Data, "rho", 1.0);
      --  ZONAL adjustables, absent in older files: neutral values
      D.B.ztau := Get_Float_Field (Data, "ztau", 1.0);
      D.B.zdist := Get_Float_Field (Data, "zdist", 1.0);
      D.B.ztilt := Get_Float_Field (Data, "ztilt", 0.0);
      D.B.ssaA := Get_Float_Field (Data, "ssaA", 0.0);
      D.B.ssaP := Get_Float_Field (Data, "ssaP", 0.0);
      --  second pole of the impulse response, absent in older files: 0 = off
      D.B.mB := Get_Float_Field (Data, "mb", 0.0);
      --  quarter-year impulse copy, absent in older files: 0 = off
      D.B.delC := Get_Float_Field (Data, "delC", 0.0);
      --  NM=1 Bessel 2k cosine coefficient, split from bg (2026-11-02). Absent in
      --  older files: take bg, which is what that term used, so they reproduce.
      D.B.bg2 := Get_Float_Field (Data, "bg2", D.B.bg);
      if Include_LPAP then
         Read_JSON_LPAP (Data, "LPAP", D);
         Read_JSON_LPAP (Data, "lpap", D);
      end if;
      Read_JSON_Float_Array (Data, "ltep", D.B.LT);
      Read_JSON_Float_Array (Data, "LT", D.B.LT);
      Read_JSON_Float_Array (Data, "lt", D.B.LT);
      Read_JSON_Int_Array (Data, "harm", D.C);
      Read_JSON_Int_Array (Data, "harmonics", D.C);
      --  besk (2026-11-17): the Bessel correction's wavenumber, a manifold
      --  parameter of its own. Before, the correction used ltep (NM) -- or
      --  ltep (NM - 1) with LOCKF -- so a file without besk is seeded with
      --  that value and reproduces exactly.
      if Kind (Data) = JSON_Object_Type and then Has_Field (Data, "besk") then
         D.B.besK := Get_Number (Get (Data, "besk"));
      else
         declare
            NM : constant Integer :=
              Integer'Max (1, Integer'Min (D.B.LT'Last, GEM.Getenv ("NM", D.B.LT'Length)));
            Lock : constant Boolean := GEM.Getenv ("LOCKF", False);
         begin
            D.B.besK := D.B.LT (if Lock and then NM >= 2 then NM - 1 else NM);
         end;
      end if;
      for I in D.A.LP'Range loop
         D.A.LP (I) := GEM.LTE.LP (I);
      end loop;
      return True;
   exception
      when others =>
         Ada.Text_IO.Put_Line ("JSON read error " & Name);
         return False;
         -- GNAT.OS_Lib.Os_Exit(0);
   end Read_JSON;

   
   procedure Read (D : in out Param_S) is
      Exec : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name);
      Base : constant String :=
        (if Exec'Length > 0 and then Exec (Exec'Last) = '.' then
           Exec (Exec'First .. Exec'Last - 1)
         else Exec);
      FN_JSON : constant String := Base & ".p";
      FT : Ada.Text_IO.File_Type;
   begin
      -- Initialize D.A.LP with canonical Doodson-calculated periods BEFORE any JSON reading
      for I in D.A.LP'Range loop
         D.A.LP (I) := GEM.LTE.LP (I);
      end loop;
      
      -- Detect NM/NH from JSON to override defaults with actual JSON array sizes
      Detect_NM_NH_From_JSON (FN_JSON);

         Ada.Text_IO.Put_Line ("*** JSON-ONLY MODE ***");
         Ada.Text_IO.Put_Line ("Reading primary file: " & FN_JSON);
         if not Read_JSON (FN_JSON, D, True) then
            raise Ada.Text_IO.Name_Error with "JSON file not found: " & FN_JSON;
         end if;
         Ada.Text_IO.Put_Line ("Successfully loaded: " & FN_JSON);
         --  The per-CLIMATE_INDEX secondary file (FN2_JSON, "<exe>.<index>.p")
         --  used to be read here too, AFTER the primary -- but Read_JSON's
         --  Include_LPAP flag only ever gates whether LPAP gets parsed; every
         --  scalar, ltep, AND harm field gets applied unconditionally
         --  regardless of that flag, and the secondary file never contains
         --  an "lpap" key at all. So this second read did nothing but
         --  silently OVERWRITE the just-loaded primary's scalars/ltep/harm
         --  with the secondary file's OWN (frequently stale -- e.g. copied
         --  in from a different run, or simply never refreshed alongside a
         --  later primary-only update) values -- confirmed directly this
         --  session on a live cell (pdo): the secondary file held an older
         --  fit than lt.exe.p, and re-running silently regressed the loaded
         --  state to it (2026-09-29). Removed entirely per the user's own
         --  diagnosis -- lt.exe.p already contains everything needed.
         Ada.Text_IO.Put_Line ("*** JSON FILE LOADED SUCCESSFULLY ***");


   exception
      when others =>
         if Ada.Text_IO.Is_Open (FT) then
            Ada.Text_IO.Close (FT);
         end if;
   end Read;

   --  =========================================================================
   --  Load: Initialize parameters from file or legacy binary format
   --
   --  Command-line flags control behavior:
   --    -j : JSON-only mode - read only .p files, fail if not found (for debugging)
   --    -l : Load from legacy .parms binary file
   --    -p : Dump parameters and exit
   --    -w : Write parameters to .par files and exit
   --  =========================================================================
   function Leg_File (Leg : in Positive) return String is
      Exec : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name);
      Base : constant String :=
        (if Exec'Length > 0 and then Exec (Exec'Last) = '.' then
           Exec (Exec'First .. Exec'Last - 1)
         else Exec);
      L : constant String := Positive'Image (Leg);
   begin
      return Base & ".p." & L (L'First + 1 .. L'Last);
   end Leg_File;

   --  The leg's own scalars beyond IR, under their .p keys (one list for
   --  Read_Leg and Write_Leg). Own_Params is every Param_B scalar but IR and
   --  year (Param_B_Overlay's 28 less those 2), plus cw1/cw2 from Param_S.
   procedure Own_Fields
     (D : in out Param_S; Own : in Leg_Own_Set;
      Apply : not null access procedure (Key : in String; V : in out Long_Float)) is
   begin
      if Own = Own_None then
         return;
      end if;
      Apply ("impA", D.B.ImpA);
      Apply ("impB", D.B.ImpB);
      Apply ("besk", D.B.besK);
      Apply ("offs", D.B.Offset);
      Apply ("bg", D.B.bg);
      Apply ("bg2", D.B.bg2);
      if Own = Own_Params then
         Apply ("impC", D.B.ImpC);
         Apply ("delA", D.B.DelA);
         Apply ("delB", D.B.DelB);
         Apply ("asym", D.B.Asym);
         Apply ("ann1", D.B.Ann1);
         Apply ("ann2", D.B.Ann2);
         Apply ("sem1", D.B.Sem1);
         Apply ("sem2", D.B.Sem2);
         Apply ("ma", D.B.mA);
         Apply ("mp", D.B.mP);
         Apply ("shfT", D.B.shiftT);
         Apply ("init", D.B.init);
         Apply ("rho", D.B.rho);
         Apply ("ztau", D.B.ztau);
         Apply ("zdist", D.B.zdist);
         Apply ("ztilt", D.B.ztilt);
         Apply ("ssaA", D.B.ssaA);
         Apply ("ssaP", D.B.ssaP);
         Apply ("mb", D.B.mB);
         Apply ("delC", D.B.delC);
         Apply ("cw1", D.CW1);
         Apply ("cw2", D.CW2);
      end if;
   end Own_Fields;

   procedure Read_Leg (Leg : in Positive; D : in out Param_S; Found : out Boolean;
                       Own : in Leg_Own_Set := Own_None)
   is
      use GNATCOLL.JSON;
      Name : constant String := Leg_File (Leg);
      Result : Read_Result;
      procedure Get (Key : in String; V : in out Long_Float) is
      begin
         V := Get_Float_Field (Result.Value, Key, V);
      end Get;
   begin
      Found := False;
      if not Ada.Directories.Exists (Name) then
         return;
      end if;
      Result := Read (File_To_String (Name));
      if not Result.Success then
         Ada.Text_IO.Put_Line (Name & ": " & Format_Parsing_Error (Result.Error));
         GNAT.OS_Lib.OS_Exit (0);
      end if;
      D.B.IR := Get_Float_Field (Result.Value, "IR", D.B.IR);
      Read_JSON_Float_Array (Result.Value, "ltep", D.B.LT);
      Read_JSON_Int_Array (Result.Value, "harm", D.C);
      Own_Fields (D, Own, Get'Access);
      Found := True;
   end Read_Leg;

   procedure Write_Leg (Leg : in Positive; D : in Param_S;
                        Own : in Leg_Own_Set := Own_None) is
      use GNATCOLL.JSON;
      NM : constant Integer := GEM.Getenv ("NM", D.B.LT'Length);
      Harms : constant Ns := Parse_NH (GEM.Getenv ("NH", ""));
      Obj : constant JSON_Value := Create_Object;
      LT_Arr, Harm_Arr : JSON_Array := Empty_Array;
      FT : Ada.Text_IO.File_Type;
      DD : Param_S := D;
      procedure Put (Key : in String; V : in out Long_Float) is
      begin
         Set_Field (Obj, Key, Create (V));
      end Put;
   begin
      Own_Fields (DD, Own, Put'Access);
      Set_Field (Obj, "IR", Create (D.B.IR));
      for I in 1 .. NM loop
         Append (LT_Arr, Create (D.B.LT (I)));
      end loop;
      Set_Field (Obj, "ltep", Create (LT_Arr));
      for I in 1 .. Harms'Length loop
         exit when I > D.C'Last;
         Append (Harm_Arr, Create (Long_Float (D.C (I))));
      end loop;
      Set_Field (Obj, "harm", Create (Harm_Arr));
      Ada.Text_IO.Create (FT, Ada.Text_IO.Out_File, Leg_File (Leg));
      Ada.Text_IO.Put_Line (FT, Write (Obj, Compact => False));
      Ada.Text_IO.Close (FT);
   end Write_Leg;

   function Wisdom_File return String is
      Env : constant String := GEM.Getenv ("WISDOM_P", "");
      Rel : constant String := "experiments/wisdom_lpap/wisdom_raw.lt.exe.p";
   begin
      if Env /= "" then
         return Env;
      end if;
      declare
         Dir : Ada.Strings.Unbounded.Unbounded_String :=
           Ada.Strings.Unbounded.To_Unbounded_String
             (Ada.Directories.Current_Directory);
      begin
         for Level in 0 .. 10 loop
            declare
               Here : constant String := Ada.Strings.Unbounded.To_String (Dir);
               Name : constant String := Here & "/" & Rel;
            begin
               if Ada.Directories.Exists (Name) then
                  return Name;
               end if;
               exit when Here = "/" or else Here = "";
               Dir := Ada.Strings.Unbounded.To_Unbounded_String
                 (Ada.Directories.Containing_Directory (Here));
            end;
         end loop;
      exception
         when Ada.Directories.Use_Error | Ada.Directories.Name_Error =>
            null;
      end;
      return "";
   end Wisdom_File;

   --  The first experiments/wisdom_lpap at or above the run folder ("" none)
   function Wisdom_Dir return String is
      Dir : Ada.Strings.Unbounded.Unbounded_String :=
        Ada.Strings.Unbounded.To_Unbounded_String
          (Ada.Directories.Current_Directory);
   begin
      for Level in 0 .. 10 loop
         declare
            Here : constant String := Ada.Strings.Unbounded.To_String (Dir);
            Name : constant String := Here & "/experiments/wisdom_lpap";
         begin
            if Ada.Directories.Exists (Name) then
               return Name;
            end if;
            exit when Here = "/" or else Here = "";
            Dir := Ada.Strings.Unbounded.To_Unbounded_String
              (Ada.Directories.Containing_Directory (Here));
         end;
      end loop;
      return "";
   exception
      when Ada.Directories.Use_Error | Ada.Directories.Name_Error =>
         return "";
   end Wisdom_Dir;

   function Lower (S : String) return String is
      R : String := S;
   begin
      for I in R'Range loop
         if R (I) in 'A' .. 'Z' then
            R (I) := Character'Val (Character'Pos (R (I)) + 32);
         end if;
      end loop;
      return R;
   end Lower;

   function Wisdom_Named (Name : String) return String is
      use Ada.Directories;
      Dir : constant String := Wisdom_Dir;
      Want : constant String := Lower (Name) & ".lt.exe.p";
      Found : Ada.Strings.Unbounded.Unbounded_String;
      S : Search_Type;
      E : Directory_Entry_Type;
   begin
      if Dir = "" then
         return "";
      end if;
      Start_Search (S, Dir, "*.lt.exe.p", (Ordinary_File => True, others => False));
      while More_Entries (S) loop
         Get_Next_Entry (S, E);
         if Lower (Simple_Name (E)) = Want then
            Found := Ada.Strings.Unbounded.To_Unbounded_String (Full_Name (E));
         end if;
      end loop;
      End_Search (S);
      return Ada.Strings.Unbounded.To_String (Found);
   exception
      when Use_Error | Name_Error =>
         return "";
   end Wisdom_Named;

   function Wisdom_Names return String is
      use Ada.Directories;
      Dir : constant String := Wisdom_Dir;
      List : Ada.Strings.Unbounded.Unbounded_String;
      S : Search_Type;
      E : Directory_Entry_Type;
   begin
      if Dir = "" then
         return "(no experiments/wisdom_lpap found)";
      end if;
      Start_Search (S, Dir, "*.lt.exe.p", (Ordinary_File => True, others => False));
      while More_Entries (S) loop
         Get_Next_Entry (S, E);
         declare
            N : constant String := Simple_Name (E);
         begin
            Ada.Strings.Unbounded.Append (List, " " & N (N'First .. N'Last - 9));
         end;
      end loop;
      End_Search (S);
      return Dir & ":" & Ada.Strings.Unbounded.To_String (List);
   exception
      when Use_Error | Name_Error =>
         return Dir;
   end Wisdom_Names;

   procedure Seed_LPAP (D : in out Param_S; File : in String;
                        Seeded, Skipped : out Natural)
   is
      use GNATCOLL.JSON;
      Result : Read_Result;
      Arr : JSON_Array;
      Used : array (D.B.LPAP'Range) of Boolean := (others => False);
   begin
      Seeded := 0;
      Skipped := 0;
      Result := Read (File_To_String (File));
      if not Result.Success then
         raise Constraint_Error with
           File & ": " & Format_Parsing_Error (Result.Error);
      end if;
      if not Has_Field (Result.Value, "lpap") then
         raise Constraint_Error with File & ": no lpap";
      end if;
      Arr := Get (Result.Value, "lpap");
      for I in 1 .. Length (Arr) loop
         declare
            T : constant JSON_Array := Get (Get (Arr, I));
            Period, Amp, Phase : Long_Float;
            Best : Integer := 0;
            Best_Diff : Long_Float := Long_Float'Last;
         begin
            if Length (T) >= 3 then
               Period := Get_Number (Get (T, 1));
               Amp := Get_Number (Get (T, 2));
               Phase := Get_Number (Get (T, 3));
               if Amp = 0.0 then
                  Skipped := Skipped + 1;
               else
                  for J in D.B.LPAP'Range loop
                     if not Used (J) and then D.A.LP (J) /= 0.0 then
                        declare
                           Diff : constant Long_Float :=
                             abs (abs Period - abs D.A.LP (J));
                        begin
                           if Diff <= 1.0E-4 * abs D.A.LP (J)
                             and then Diff < Best_Diff
                           then
                              Best := J;
                              Best_Diff := Diff;
                           end if;
                        end;
                     end if;
                  end loop;
                  if Best = 0 then
                     Ada.Text_IO.Put_Line
                       ("  WISDOM: no table period for" & Period'Img);
                     Skipped := Skipped + 1;
                  else
                     D.B.LPAP (Best) := (Amplitude => Amp, Phase => Phase);
                     Used (Best) := True;
                     Seeded := Seeded + 1;
                  end if;
               end if;
            end if;
         end;
      end loop;
   end Seed_LPAP;

   procedure Load (P : in out Param_S) is
   begin
      Read (P);

      if GEM.Command_Line_Option_Exists ("p") then
         Dump (P);
         GNAT.OS_Lib.OS_Exit (0);
      end if;
   end Load;

   --  =========================================================================
   --  Dump: Print all parameters in human-readable markdown format
   --
   --  Outputs parameters with percentage differences from reference values.
   --  Used with -p command-line flag for parameter inspection.
   --  =========================================================================
   procedure Dump (D : in Param_S) is
      function Percent (A, B : in Long_Float) return String is
      begin
         if A = 0.0 then
            return ", 0.0";
         else
            return ", " & Integer'Image (Integer ((B - A) / A * 100.0));
         end if;
      end Percent;
   begin
      Ada.Text_IO.Put_Line ("```");
      Put (D.B.Offset, " :offset:", NL);
      Put (D.B.bg, " :bg:", NL);
      Put (D.B.ImpA, " :impA:", NL);
      Put (D.B.ImpB, " :impB:", NL);
      Put (D.B.ImpC, " :impC:", NL);
      Put (D.B.DelA, ":delA:", NL);
      Put (D.B.DelB, ":delB:", NL);
      Put (D.B.Asym, ":asym:", NL);
      Put (D.B.Ann1, ":ann1:", NL);
      Put (D.B.Ann2, ":ann2:", NL);
      Put (D.B.Sem1, ":sem1:", NL);
      Put (D.B.Sem2, ":sem2:", NL);
      Put (D.CW1, ":cw1:", NL);
      Put (D.CW2, ":cw2:", NL);
      Put (D.B.Year, ":year:", NL);
      Put (D.B.IR, ":IR:", NL);
      Put (D.B.mA, " :mA:", NL);
      Put (D.B.mP, " :mP:", NL);
      Put (D.B.shiftT, " :shiftT:", NL);
      Put (D.B.init, " :init:", NL);
      Put (D.B.rho, " :rho:", NL);
      Put (D.B.ztau, " :ztau:", NL);
      Put (D.B.zdist, " :zdist:", NL);
      Put (D.B.ztilt, " :ztilt:", NL);
      Put (D.B.ssaA, " :ssaA:", NL);
      Put (D.B.ssaP, " :ssaP:", NL);
      Put (D.B.mB, " :mB:", NL);
      Put (D.B.bg2, " :bg2:", NL);
      Put (D.B.delC, " :delC:", NL);
      Put (D.B.besK, " :besk:", NL);
      Ada.Text_IO.Put_Line ("---- Tidal ----");
      for I in D.B.LPAP'Range loop
         Put (D.A.LP (I), ", ");
         Put (D.B.LPAP (I).Amplitude, ", ");
         Put
           (D.B.LPAP (I).Phase,
            ", " & I'Img &
            Percent (GEM.LTE.LPRef (I).Amplitude, D.B.LPAP (I).Amplitude) &
            ", " & GEM.LTE.LPRef (I).Amplitude'Img,
            NL);
      end loop;
      --  TODO: Can remove - Debug output for LTE static parameters was used
      --  during development to verify parameter loading. Now redundant as main
      --  Dump output shows all critical parameters.
      --Ada.Text_IO.Put_Line("---- LTE static ----");
      --for I in D.B.LT'Range loop
      --   Put(D.B.LT(I), "", NL);
      --end loop;
   end Dump;

end GEM.LTE.Primitives.Shared;
