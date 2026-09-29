--  Manifold_Regularizer_Test -- standalone self-test for the MANIFOLD
--  quadrature-flip regularizer added to GEM.LTE.Primitives.Solution's
--  local Metric function (2026-09-29, per user request): "regularize
--  the metric via a correlation coefficient that compares to a mean or
--  median manifold ... If CC~1.0, no penalty, but if CC~0, large
--  penalty."
--
--  The actual penalty logic lives in a function local to Dipole_Model
--  (not separately callable), so this test instead exercises the two
--  building blocks it's built from directly -- GEM.LTE.Primitives.CC
--  (already used throughout this project for every reported CorrCoeff)
--  and the Data_Pairs slicing convention Metric relies on (Manifold_Ref
--  indexed by the caller's own Z'First..Z'Last, since every Metric
--  call's Z argument is always some slice of Forcing sharing Data_Records'
--  own bounds) -- and checks the exact three properties the penalty
--  depends on:
--
--   1. A manifold compared against ITSELF gives CC ~= +1.0 (multiplier
--      ~= 1.0, i.e. no penalty) -- the "already in the right quadrature"
--      case.
--   2. A manifold compared against its own NEGATION gives CC ~= -1.0
--      (multiplier clamps to 0.0 via Long_Float'Max(0.0, CC) in the real
--      penalty code) -- the "wrong quadrature" case the regularizer
--      exists to suppress. (A genuine ~90 degree quadrature swap, the
--      case actually seen in this project's own fits, gives CC nearer
--      0 than -1; -1 is the more extreme 180 degree case, tested here
--      because it's exactly reproducible without needing a real
--      quadrature-ambiguous dataset -- both land in the same "clamps
--      toward/through zero" regime.)
--   3. Indexing the reference by a slice that runs past its own bounds
--      raises Constraint_Error -- confirming the guard Metric relies on
--      (an exception handler that leaves Raw_Score unpenalized rather
--      than propagating) is reachable for a real misaligned-reference
--      scenario, not just in principle.
with Text_IO;
with GEM.LTE.Primitives; use GEM.LTE.Primitives;

procedure Manifold_Regularizer_Test is
   N : constant := 200;
   Self_Ref, Flipped_Ref, Slice : Data_Pairs (1 .. N);
   Pass : Boolean := True;

   procedure Check (Label : String; Got, Want : Long_Float; Tol : Long_Float) is
   begin
      if abs (Got - Want) <= Tol then
         Text_IO.Put_Line ("PASS  " & Label & " =" & Got'Img);
      else
         Text_IO.Put_Line
           ("FAIL  " & Label & " =" & Got'Img & "  (want" & Want'Img &
            " +/-" & Tol'Img & ")");
         Pass := False;
      end if;
   end Check;

begin
   --  A synthetic "manifold" with real variation (not a flat line, which
   --  CC's own zero-variance guard would report as 0.0 regardless).
   for I in 1 .. N loop
      Self_Ref (I) :=
        (Date => 1880.0 + Long_Float (I) / 12.0,
         Value => Long_Float (I mod 17) - 8.0 + 0.3 * Long_Float (I));
      Flipped_Ref (I) := (Date => Self_Ref (I).Date, Value => -Self_Ref (I).Value);
   end loop;
   Slice := Self_Ref;

   Check ("self-vs-self CC", CC (Slice, Self_Ref), 1.0, 1.0e-9);
   Check ("self-vs-negated CC", CC (Slice, Flipped_Ref), -1.0, 1.0e-9);
   Check
     ("penalty multiplier on self match",
      Long_Float'Max (0.0, CC (Slice, Self_Ref)), 1.0, 1.0e-9);
   Check
     ("penalty multiplier on flipped match",
      Long_Float'Max (0.0, CC (Slice, Flipped_Ref)), 0.0, 1.0e-9);

   --  Mirrors Metric's own exact nesting (see gem-lte-primitives-
   --  solution.adb): an exception raised while ELABORATING a block's own
   --  declarative part is NOT caught by that same block's handler in
   --  Ada -- it propagates to the enclosing construct instead. The real
   --  Manifold_CC declaration has exactly this shape (`Manifold_Ref
   --  (Z'First .. Z'Last)` inside a `declare` block's own decl part), so
   --  the guard needs an OUTER bare block wrapping the inner `declare`
   --  as a single statement, not a handler on the `declare` block
   --  itself -- confirmed live: the single-block version of this test
   --  raised CONSTRAINT_ERROR uncaught instead of printing PASS below.
   begin
      declare
         Bad : Data_Pairs := Self_Ref (1 .. N + 50);
         pragma Unreferenced (Bad);
      begin
         Text_IO.Put_Line
           ("FAIL  out-of-bounds slice did not raise Constraint_Error");
         Pass := False;
      end;
   exception
      when Constraint_Error =>
         Text_IO.Put_Line
           ("PASS  out-of-bounds slice raised Constraint_Error " &
            "(the guard Metric's exception handler relies on)");
   end;

   Text_IO.New_Line;
   if Pass then
      Text_IO.Put_Line ("ALL CHECKS PASSED");
   else
      Text_IO.Put_Line ("SOME CHECKS FAILED");
   end if;
end Manifold_Regularizer_Test;
