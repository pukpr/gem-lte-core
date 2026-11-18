with GEM.LTE.Primitives.Shared;

--  INVERSION=TRUE (2026-11-08): deterministic refinement of the loaded
--  lt.exe.p in place of the random search. One task; every evaluation is
--  the production one (forcing, Bessel, Regression_Factors, LTE, IR, F9),
--  so all response-file options apply unchanged. Per iteration:
--    1. the current base parameters are evaluated (and accepted/saved as
--       an ordinary search would);
--    2. each active parameter -- Re and Im of z = A e^(i phase) for every
--       LPAP line with nonzero amplitude, plus ImpA and ImpB -- is probed
--       by a finite difference (one evaluation each);
--    3. a ridge-regularized step is solved and a grid of ridge strengths x
--       step sizes is evaluated; the candidate with the best TRAINING score
--       becomes the next base if it beats the current one, otherwise the
--       run stops. The held-out score is printed, never used to choose.
--  INVERSION_METHOD GN (default): Gauss-Newton on the training residual
--    Data - Model, the Jacobian taken through the full evaluation (the
--    regression is refit at every probe, i.e. variable projection).
--  INVERSION_METHOD VITERBI: a Viterbi path dF(t) of the manifold (prior
--    INVERSION_PRIOR and step INVERSION_STEP_SD, both in fastest-winding
--    periods) that best explains the training data with the current
--    regression held fixed, then the step that best reproduces dF from the
--    manifold Jacobian.
--  INVERSION_ITERS (default 4), INVERSION_RIDGE (default "0.001 0.01 0.1"),
--  INVERSION_STEPS (default "0.25 0.5 1").
--  LEGS (2026-11-09): one group; the leader runs the iteration, every leg
--  records its own rows (weighted by sqrt (LEG_WEIGHT) / its base residual
--  rms) and the step is solved from all of them; the candidates are scored
--  by the pair's composite. GN only (VITERBI falls back to GN).
--  Validation (after TRAIN_END, as VALIDATE scores it) and test (before
--  TRAIN_START) scores are printed with every base and candidate; they are
--  out-of-sample only with EXCLUDE=FALSE. INVERSION_STOP VALID (2026-11-09):
--  steps are still accepted on training, but the base with the best
--  validation score is kept, the run stops INVERSION_PATIENCE (default 2)
--  iterations after it, and that base is re-evaluated and saved.
--  Not available with COVERAGE, ENCLOSING or SPLIT_TRAINING.
package GEM.LTE.Primitives.Inversion is

   On : constant Boolean;

   type Mask is array (Integer range <>) of Boolean;

   --  The evaluation this leg (0 or 1 = the leader) just recorded was of the
   --  current base (to be accepted). Per leg: a follower asks after the
   --  composite is in, when the leader may already have moved on.
   function Is_Base (Leg : Natural := 0) return Boolean;

   --  Called after every evaluation. Pointwise (I, F) is the regression
   --  model at sample I with manifold value F (no lag, no filter);
   --  Center_Weight is the filter's weight on a sample's own value.
   procedure Record_Result
     (D : in Shared.Param_S; Model, Data, Forcing : in Data_Pairs;
      Train : in Mask; Score, Held : in Long_Float; Singular : in Boolean;
      Winding_Period, Center_Weight : in Long_Float;
      Pointwise : not null access function
        (I : Integer; F : Long_Float) return Long_Float;
      Leg : in Natural := 0; Leg_Weight : in Long_Float := 1.0;
      Valid, Test : in Long_Float := Long_Float'First);

   --  The evaluation this leg just recorded is the final best-validation
   --  base (INVERSION_STOP VALID): accept it whatever its training score
   function Force_Accept (Leg : Natural := 0) return Boolean;

   --  Sets D to the next parameters to evaluate; Done when finished
   procedure Next (D : in out Shared.Param_S; Done : out Boolean);

private

   On : constant Boolean := GEM.Getenv ("INVERSION", False);

end GEM.LTE.Primitives.Inversion;
