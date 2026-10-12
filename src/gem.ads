package GEM is
   -- Top level directory for all GeoEnergyMath software routines
   -- Adapted from Mathematical GeoEnergy (Wiley/2009)
   -- All GEM child libraries inherit from this adaptation

   function Getenv (Name : in String; Default : in String) return String;
   function Getenv (Name : in String; Default : in Integer) return Integer;
   function Getenv
     (Name : in String; Default : in Long_Integer) return Long_Integer;
   function Getenv (Name : in String; Default : in Float) return Float;
   function Getenv
     (Name : in String; Default : in Long_Float) return Long_Float;
   function Getenv (Name : in String; Default : in Boolean) return Boolean;

   procedure Setenv (Name : in String; Value : in String);

   --  LEGS: a worker task can be given a leg number (1 or 2). Getenv, when
   --  called from that task, then looks first in the leg's own response file
   --  "<executable>.resp.<leg>", then in the common response file, then in
   --  the environment. Leg 0 (the default for every task) is the old lookup.
   procedure Set_Leg (L : in Natural);
   function Leg return Natural;
   function Leg_Tag return String;  -- "" for leg 0, else ".1" / ".2"
   procedure Clear (Name : in String);

   -- List of integers from "1 2 4 6" etc
   type Ns is array (Positive range <>) of Integer;
   function S_to_I (S : in String) return Ns;

   -- List of floats from "1.0 2.0 4.2 6.4" etc
   type Fs is array (Positive range <>) of Long_Float;
   function S_to_LF (S : in String) return Fs;

   function Command_Line_Option_Exists (Option : in String) return Boolean;

   type Annual_Harmonics is
      record
         Ann1, Ann2, Semi1, Semi2 : Long_Float;
      end record;

end GEM;
