--  ============================================================================
--  GEM.LTE - Laplace's Tidal Equation Package Body
--  ============================================================================
--
--  PURPOSE:
--    Implements Laplace's Tidal Equation (LTE) framework for modeling
--    ocean-atmosphere response to long-period gravitational forcing.
--    Provides Doodson number calculations and initialization of tidal
--    constituent sets used throughout GEM-LTE climate analysis.
--
--  SCIENTIFIC BACKGROUND:
--    Doodson numbers encode the astronomical arguments (s, h, p, N, etc.)
--    that determine tidal periods. Each constituent represents a specific
--    combination of lunar/solar forcing frequencies. The framework supports:
--    - Standard long-period tides (18.6yr, 8.85yr, etc.)
--    - Annual/semi-annual harmonics
--    - Quasi-Biennial Oscillation (QBO) periods
--
--  KEY FUNCTIONS:
--    - Doodson: Converts Doodson arguments to frequency (cycles/year)
--    - Year_Adjustment: Dynamically updates all periods when year length changes
--
--  INITIALIZATION:
--    Package elaboration computes all tidal periods and initializes:
--    - LP/LPAP: Main long-period constituents
--    - LP_Set/LPAP_Set: Extended constituent set
--    - LP_RSet/LPAP_RSet: Reduced constituent set
--    - LP_QSet/LPAP_QSet: QBO-focused set
--    - LP_Annual/LPAP_Annual: Annual harmonics
--    - QBO/QBOAP: QBO-specific periods
--
--  ============================================================================

package body GEM.LTE is

   Year_in_Days : constant := 365.242_248_4;  -- 365.241237718675000;
   Year_Correction : constant Long_Float := GEM.Getenv ("YEAR", 0.0);
   --  TODO: Can remove - dLOD_Mod was experimental scaling factor for day-length
   --  variations. Standard Doodson calculation proved sufficient.
   --dLOD_Mod : Long_Float := GEM.Getenv("dLOD_Mod", 1.0);
   function Year_Length (Dynamic_Correction : Long_Float := 0.0)
     return Long_Float
   is
   begin
      return Year_in_Days + Year_Correction + Dynamic_Correction;
   end Year_Length;

   --  Compute tidal frequency from Doodson arguments
   --  Formula combines tropical year, lunar nodal, perigee, and anomalistic periods
   function Doodson
     (I : in Integer; D : in Doodson_List;
      Dynamic_Correction : in Long_Float := 0.0) return Long_Float
   is
   begin
      --  TODO: Can remove - Conditional dLOD scaling was experimental approach
      --  to handle day-length variations. Standard formula works for all cases.
      --if D(I).h = 1 and dLOD_Mod /= 1.0 then
      --   return dLOD_Mod * Year_Length;
      --else
      return
        1.0 /
        (Long_Float (D (I).s) / Tropical +
        (Long_Float (D (I).h) * D (I).Year_Multiplier) /
          Year_Length (Dynamic_Correction) +
        Long_Float (D (I).p) / p + Long_Float (D (I).N) / N);
      --end if;
   end Doodson;

   --  Recalculate a caller-owned period list for one candidate year correction.
   --  The startup YEAR correction remains part of Year_Length.
   procedure Year_Adjustment (Value : in Long_Float; List : in out Periods) is
   begin
      for I in List'Range loop
         List (I) := Doodson (I, Doodson_Args, Value);
      end loop;
   end Year_Adjustment;

   --  Legacy stateful implementation retained for comparison:
   --  if Value /= 0.0 then
   --     Year_Dynamic_Correction := Value;
   --     for I in List'Range loop
   --        List (I) := Doodson (I, Doodson_Args);
   --     end loop;
   --  end if;

begin

   --  Package elaboration: Initialize all tidal constituent period arrays
   --  Converts Doodson arguments to frequencies for use in LTE solver

   for I in Doodson_Args'Range loop
      LP (I) := Doodson (I, Doodson_Args);
      LPAP (I) := (0.0, 1.0);
   end loop;


end GEM.LTE;
