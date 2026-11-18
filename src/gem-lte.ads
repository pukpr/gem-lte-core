package GEM.LTE is
   -- Laplace's Tidal Equations model for equatorial waves -- Chap 11, Chap 12

   -- Top level holds the main Tidal Constituents and LTE modulation parameters

   pragma Elaborate_Body;

   -- The year value is close to the Tropical year of 365.24219 days
   -- which accumulates to ~1/2 day error over 100 years
   -- Year : constant Long_Float := 365.2412384; -- 365.246462;
--   Draconic : constant Long_Float := 27.21222082;
--   Tropical : constant Long_Float := 27.32163237; -- 27.32166155;
--   Anomalistic : constant Long_Float := 27.55454988;

   Draconic : constant Long_Float := 27.212_220_815;
   --Draconic : constant Long_Float := 27.212_0;
   --  The Moon's mean-longitude month used for the Doodson argument s.
   --    sidereal month 27.321_661_554 d (star to star)       -> Mf 13.660831 d
   --    tropical month 27.321_582_252 d (equinox to equinox) -> Mf 13.660791 d
   --  The value used up to now (under the name Tropical) is the SIDEREAL
   --  month, and it stays the default so existing fits are unchanged.
   --  TROPICAL=TRUE (or SIDEREAL=FALSE) selects the tropical month.
   Sidereal_Month : constant Long_Float := 27.321_661_554;
   Tropical_Month : constant Long_Float := 27.321_582_252;
   Use_Tropical : constant Boolean :=
     GEM.Getenv ("TROPICAL", False) or else not GEM.Getenv ("SIDEREAL", True);
   Tropical : constant Long_Float :=
     (if Use_Tropical then Tropical_Month else Sidereal_Month);
   Anomalistic : constant Long_Float := 27.554_549_886;
   Extra : constant Long_Float := GEM.Getenv ("EXTRA", 0.405);

   --  nodal and perigee periods, as exact static values for each month so
   --  the default reproduces earlier results bit for bit
   N_Sid : constant := 1.0 / (1.0 / Draconic - 1.0 / Sidereal_Month);
   p_Sid : constant := 1.0 / (1.0 / Sidereal_Month - 1.0 / Anomalistic);
   N_Trop : constant := 1.0 / (1.0 / Draconic - 1.0 / Tropical_Month);
   p_Trop : constant := 1.0 / (1.0 / Tropical_Month - 1.0 / Anomalistic);
   N : constant Long_Float := (if Use_Tropical then N_Trop else N_Sid);
   p : constant Long_Float := (if Use_Tropical then p_Trop else p_Sid);
   --Year_Length : constant := 365.2412384;  -- 365.241237718675000;

   function Year_Length (Dynamic_Correction : Long_Float := 0.0)
     return Long_Float;

   type Doodson_Argument is record
      s, h, p : Integer;
      N : Integer;
      Year_Multiplier : Long_Float;
      --  Day (2026-11-05): cycles per mean solar day, 1 for the diurnal
      --  (tesseral) lines, 2 for the semidiurnal (sectorial) lines;
      --  0 for every long-period line.
      Day : Integer;
   end record;

   type Doodson_List is array (Positive range <>) of Doodson_Argument;

   Doodson_Args : constant Doodson_List :=
     ((1, 0, 0, 0, 1.0, 0),
      (1, 0, 0, 1, 1.0, 0),
      (0, 0, 2, 2, 1.0, 0),
      (2, 0, 0, 1, 1.0, 0),
      (2, 0, 0, 0, 1.0, 0), 
      (2, 0, 0, 2, 1.0, 0),
      (1, 0, -1, 0, 1.0, 0),
      (2, 0, -2, 0, 1.0, 0),
      (0, 0, 0, 1, 1.0, 0),
      (0, 0, 2, 0, 1.0, 0), 
      (1,-2, 1, 0, 1.0, 0), 
      (0, 0, 2, 1, 1.0, 0), 
      (1, 0, -1, 1, 1.0, 0), 
      (1, 0, -1, -1, 1.0, 0),
      (0, 0, 1, 1, 1.0, 0),
      (1, 0, 1, 1, 1.0, 0), 
      (0, 0, 1, -1, 1.0, 0),
      (0, 0, -1, 0, 1.0, 0),
      (0, 0, -2, 1, 1.0, 0),
      (3, 0, -1, 0, 1.0, 0), 
      (3, 0, -1, 1, 1.0, 0),
      (3, 0, -1, 2, 1.0, 0), 
      (0, 0, 0, 2, 1.0, 0),
      (0, 0, 1, 2, 1.0, 0),

       (3,-2, 1, 0, 1.0, 0),
       (3, 0,-3, 0, 1.0, 0), 
       (3, -2, 1, 1, 1.0, 0),
       (4, -2, 0, 1, 1.0, 0),
       (4, 0,-2, 1, 1.0, 0),
       (4, 0,-2, 0, 1.0, 0),
       (4, -2, 0, 0, 1.0, 0),
       (5, -2, -1, 0, 1.0, 0),

       (2, -2, 0, 0, 1.0, 0), 
       (3, -2, -1, 0, 1.0, 0), 
       (1,  0, 1, 0, 1.0, 0),
       (5, -2, -1, 1, 1.0, 0),
       (5,  0, -3, 0, 1.0, 0),
       (2, -2, 0, -1, 1.0, 0), 
       (2, -2, 0,  1, 1.0, 0), 
       (2, -3, 0, 0, 1.0, 0), 
       (1,  2, -1, 0, 1.0, 0),
       (0, 2, 0, 0, 1.0, 0),
       --  OO1 (185.555, tau + 3s = 1 cpd + 2s + h): tesseral, aliases to ~59 yr
       (2,  1, 0, 0, 1.0, 1),
      --  Diurnal and semidiurnal lines (2026-11-05), written as
      --  Day cycles/solar day + (s, h, p) using tau = 1 cpd + h - s
      --  K1 (165.555, tau + s = 1 cpd + h)
      (0,  1, 0, 0, 1.0, 1),
      --  O1 (145.555, tau - s = 1 cpd - 2s + h)
      (-2, 1, 0, 0, 1.0, 1),
      --  M1 (155.555, tau = 1 cpd - s + h), mean lunar day
      (-1, 1, 0, 0, 1.0, 1),
      --  J1 (175.455, tau + 2s - p = 1 cpd + s + h - p)
      (1,  1, -1, 0, 1.0, 1),
      --  Q1 (135.655, tau - 2s + p = 1 cpd - 3s + h + p)
      (-3, 1, 1, 0, 1.0, 1),
      --  P1 (163.555, tau + s - 2h = 1 cpd - h)
      (0, -1, 0, 0, 1.0, 1),
      --  M2 (255.555, 2 tau = 2 cpd - 2s + 2h)
      (-2, 2, 0, 0, 1.0, 2),
      --  N2 (245.655, 2 tau - s + p = 2 cpd - 3s + 2h + p)
      (-3, 2, 1, 0, 1.0, 2),
      --  K2 (275.555, 2 tau + 2s = 2 cpd + 2h)
      (0,  2, 0, 0, 1.0, 2)

      );

   type Amp_Phase is record
      Amplitude, Phase : Long_Float;
   end record;
   subtype Period is Long_Float;

   type Amp_Phases is array (Positive range <>) of Amp_Phase;
   type Periods is array (Positive range <>) of Period;

   -- Tidal Factor types
   subtype Long_Periods_Amp_Phase is Amp_Phases;
   subtype Long_Periods is Periods;

   -- LTE modulation types
   subtype Modulations_Amp_Phase is Amp_Phases;
   subtype Modulations is Periods;

   type Period_Set (N : Positive) is record
      LP : Periods (1 .. N);
      AP : Amp_Phases (1 .. N);
   end record;

   -- Default values that produce a fit to NINO34 of cc = ~0.83
   -- The period values can be derived from the constants above
   -- but are left as values for easy identification.
   -- Also, the phases are not reduced to within (-PI .. PI)  boundaries

   --  SemiAnnual : constant := 182.623231; -- days
   --  Annual : constant := 365.246462;
   --  ThirdAnnual : constant := 121.748821;
   --  Msm : constant := 31.81209136;
   --  Mm_1 : constant := 27.66676713;
   --  Mm_f : constant := 13.77727494;  -- _f for fortnightly
   --  Nodical_f : constant := 13.60611041;
   --  Mf_f : constant := 13.66083077;
   --  Nodical : constant := 27.21222082;
   --  Mm : constant := 27.55454988;
   --  Mf : constant := 27.32166155;
   --  Mf_prime : constant := 13.63341568;
   --  Msqm : constant := 7.095810615; -- aka Msq
   --  Mqm : constant := 6.859402884;  -- Mq
   --  Mstm : constant := 9.556887401;
   --  Mfm : constant := 9.132950781;  -- aka Mt
   --  Mfm_1 : constant := 9.124566355;
   --  Msf_f : constant := 14.76532679;
   --  Msm_1 : constant := 27.09267692;
   --  Msf : constant := 29.53065358;
   --  Perigee_half : constant := SemiAnnual*8.85;
   --  Perigee : constant := Annual*8.85;
   --  Nodal_half : constant := SemiAnnual*18.6;
   --  Nodal : constant := Annual*18.6;
   --  Mf_Mfp : constant := 6.82355473;
   --  Mf_Mf : constant := 6.83041539;
   --  Mf_Nodical : constant := 6.81670784;
   --  Msp : constant := 5.643; -- not used
   --  M357 : constant := 1.0/(1.0/Nodal+1.0/Perigee_half);
   --  MTDD : constant := 1.0/(1.0/Mf+2.0/Nodical);
   --

   LPAP : Long_Periods_Amp_Phase (Doodson_Args'Range);
   LP : Long_Periods (Doodson_Args'Range);

   LPRef : Long_Periods_Amp_Phase (Doodson_Args'Range);

   LTM : constant Modulations :=
     (2.763_168_894, 42.387_030_79, 18.682_160_48, 100.0, 10.0, 200.0, 600.0,
      70.0, 150.0, 400.0, 0.5);

   LTAP : constant Modulations_Amp_Phase :=
     ((0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0),
      (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0));

   LT0 : Modulations (1 .. 0);


   procedure Year_Adjustment (Value : in Long_Float; List : in out Periods);

end GEM.LTE;
