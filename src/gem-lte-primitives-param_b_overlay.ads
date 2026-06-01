--  ============================================================================
--  GEM.LTE.Primitives.Param_B_Overlay - Safe Array Overlay for Param_B
--  ============================================================================
--
--  PURPOSE:
--    Provides type-safe wrapper around memory overlay used for random descent
--    optimization. Encapsulates the unsafe 'Address clause with runtime checks.
--
--  BACKGROUND:
--    The random descent algorithm (Walker.Markov) needs to randomly modify
--    parameters without knowledge of the Param_B record structure. Rather than
--    using unsafe C-style memory overlay, this package provides:
--    - Named constants for field positions in overlay array
--    - Runtime verification of layout assumptions
--    - Type-safe interface for creating overlays
--
--  COMPILE-TIME SAFETY:
--    Static assertions verify fundamental assumptions about data representation:
--    - Long_Float is 64 bits (8 bytes) as expected
--    - Ensures portability across platforms
--
--  FIELD INDICES IN OVERLAY ARRAY:
--    Scalars (8 fields):
--      1. Offset    5. mA
--      2. DelA      6. mP
--      3. DelB      7. shiftT
--      4. Asym      8. init
--
--    Arrays (variable size):
--      9..9+NLP*2-1: LPAP (each constituent = 2 floats: Amplitude, Phase)
--      9+NLP*2..end: LT modulations
--
--  USAGE:
--    Size := Param_B_Overlay.Overlay_Size (D.NLP, D.NLT);
--    Set : LF_Array (1 .. Size);
--    for Set'Address use D.B.Offset'Address;
--    Param_B_Overlay.Verify_Layout (D.B, Set);  -- Runtime check
--
--  ============================================================================

with GEM.LTE.Primitives.Shared;

package GEM.LTE.Primitives.Param_B_Overlay is

   use GEM.LTE.Primitives.Shared;
   
   --  ==========================================================================
   --  COMPILE-TIME SAFETY CHECKS
   --  ==========================================================================
   
   --  Verify Long_Float is 64 bits as expected
   --  This is fundamental to the overlay working correctly
   pragma Compile_Time_Error 
      (Long_Float'Size /= 64,
       "Long_Float'Size must be 64 bits for overlay to work");
   
   --  Verify Long_Float alignment is reasonable
   pragma Compile_Time_Error
      (Long_Float'Alignment > 8,
       "Long_Float'Alignment unexpectedly large");
   
   --  ==========================================================================
   
   --  Named constants for scalar field positions in overlay array
   Offset_Index  : constant := 1;
   DelA_Index    : constant := 2;
   DelB_Index    : constant := 3;
   Asym_Index    : constant := 4;
   MA_Index      : constant := 5;
   MP_Index      : constant := 6;
   ShiftT_Index  : constant := 7;
   Init_Index    : constant := 8;

   --  Overlay layout: 8 scalars + NLP*2 (LPAP amp+phase pairs) + NLT (LT periods)
   --  First_LPAP_Index = 9
   --  First_LT_Index(NLP) = 9 + NLP*2
   Scalar_Field_Count : constant := 8;

   --  Calculate overlay array size for given discriminants
   --  Formula: Scalar_Fields + (NLP * 2) + NLT
   function Overlay_Size (NLP, NLT : Integer) return Positive;
   
   --  Verify that overlay array correctly maps to Param_B record
   --  Raises Constraint_Error if layout assumptions are violated
   --  This should be called once at initialization for safety
   procedure Verify_Layout (
      P : in Param_B;
      Expected_Size : in Positive);
   
   --  Get index of first LPAP element (Amplitude of constituent 1)
   --  LPAP indices: First_LPAP_Index .. First_LPAP_Index + NLP*2 - 1
   function First_LPAP_Index return Positive;
   
   --  Get index of first LT element
   --  LT indices: First_LT_Index(NLP) .. First_LT_Index(NLP) + NLT - 1
   function First_LT_Index (NLP : Integer) return Positive;
   
   --  Get index for specific LPAP constituent's amplitude
   --  Constituent 1..NLP, returns index for Amplitude field
   function LPAP_Amplitude_Index (Constituent : Positive) return Positive;
   
   --  Get index for specific LPAP constituent's phase
   --  Constituent 1..NLP, returns index for Phase field
   function LPAP_Phase_Index (Constituent : Positive) return Positive;
   
   --  ==========================================================================
   --  DEBUG MODE - Field Mapping Verification
   --  ==========================================================================
   
   --  Long_Float array type for overlay debugging
   type LF_Array is array (Positive range <>) of Long_Float;
   
   --  Print detailed field-by-field mapping and verify overlay access
   --  Only runs if environment variable OVERLAY_DEBUG=1
   --  Useful for verifying overlay correctness when adding new fields
   procedure Debug_Print_Field_Mapping (
      P : in Param_B;
      Set : in LF_Array);

end GEM.LTE.Primitives.Param_B_Overlay;
