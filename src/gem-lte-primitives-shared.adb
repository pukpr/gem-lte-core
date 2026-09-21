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

   CI : constant String := GEM.Getenv (Name => "CLIMATE_INDEX", Default => "");

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
      FN2 : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) & "." &
        CI & ".p";
      
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
         Set_Field (Obj, "year", Create (D.B.Year));
         Set_Field (Obj, "IR", Create (D.B.IR));
         Set_Field (Obj, "ma", Create (D.B.mA));
         Set_Field (Obj, "mp", Create (D.B.mP));
         Set_Field (Obj, "shfT", Create (D.B.shiftT));
         Set_Field (Obj, "init", Create (D.B.init));
         
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
      Write_To_File (FN, Include_LPAP => True);   -- Primary file with LPAP
      Write_To_File (FN2, Include_LPAP => False); -- Secondary file without LPAP
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
      B                           : in Param_B)
   is
      use GNATCOLL.JSON;
      FN : constant String :=
        Ada.Directories.Simple_Name (Ada.Command_Line.Command_Name) &
        ".windings.json";
      Obj      : constant JSON_Value := Create_Object;
      K_Arr    : JSON_Array          := Empty_Array;
      Manifold : constant JSON_Value := Create_Object;
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
      D.B.Year := Get_Float_Field (Data, "year", D.B.Year);
      D.B.IR :=
        Get_Float_Field (Data, "IR", Get_Float_Field (Data, "ir", D.B.IR));
      D.B.mA := Get_Float_Field (Data, "ma", D.B.mA);
      D.B.mP := Get_Float_Field (Data, "mp", D.B.mP);
      D.B.shiftT :=
        Get_Float_Field
          (Data, "shfT", Get_Float_Field (Data, "shft", D.B.shiftT));
      D.B.init := Get_Float_Field (Data, "init", D.B.init);
      if Include_LPAP then
         Read_JSON_LPAP (Data, "LPAP", D);
         Read_JSON_LPAP (Data, "lpap", D);
      end if;
      Read_JSON_Float_Array (Data, "ltep", D.B.LT);
      Read_JSON_Float_Array (Data, "LT", D.B.LT);
      Read_JSON_Float_Array (Data, "lt", D.B.LT);
      Read_JSON_Int_Array (Data, "harm", D.C);
      Read_JSON_Int_Array (Data, "harmonics", D.C);
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
      FN2_JSON : constant String := Base & "." & CI & ".p";
      FT : Ada.Text_IO.File_Type;
   begin
      -- Initialize D.A.LP with canonical Doodson-calculated periods BEFORE any JSON reading
      for I in D.A.LP'Range loop
         D.A.LP (I) := GEM.LTE.LP (I);
      end loop;
      
      -- Detect NM/NH from JSON to override defaults with actual JSON array sizes
      Detect_NM_NH_From_JSON (FN2_JSON);  -- Try secondary file first (more specific)
      Detect_NM_NH_From_JSON (FN_JSON);   -- Then primary file (may override if present)
      
         Ada.Text_IO.Put_Line ("*** JSON-ONLY MODE ***");
         Ada.Text_IO.Put_Line ("Reading primary file: " & FN_JSON);
         if not Read_JSON (FN_JSON, D, True) then
            raise Ada.Text_IO.Name_Error with "JSON file not found: " & FN_JSON;
         end if;
         Ada.Text_IO.Put_Line ("Successfully loaded: " & FN_JSON);
         
         Ada.Text_IO.Put_Line ("Reading secondary file: " & FN2_JSON);
         if not Read_JSON (FN2_JSON, D, False) then
            raise Ada.Text_IO.Name_Error with "JSON file not found: " & FN2_JSON;
         end if;
         Ada.Text_IO.Put_Line ("Successfully loaded: " & FN2_JSON);
         Ada.Text_IO.Put_Line ("*** JSON FILES LOADED SUCCESSFULLY ***");
      

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
      Put (D.B.Year, ":year:", NL);
      Put (D.B.IR, ":IR:", NL);
      Put (D.B.mA, " :mA:", NL);
      Put (D.B.mP, " :mP:", NL);
      Put (D.B.shiftT, " :shiftT:", NL);
      Put (D.B.init, " :init:", NL);
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
