with GEM.LTE.Primitives.Shared;

--  GEM.LTE.Primitives.Refine: in-process deterministic refinement of lpap,
--  per-leg Bessel/ltep/IR parameters, and harmonic sets, mirroring
--  experiments/Oct2026_generator/legs/lpap_refine.py without subprocess or
--  disk overhead.
--
--  Activated by setting REFINE=TRUE when running enso_opt (lt.exe).
--  Environment variables:
--    REFINE_LAM           (default 1.0)    dLOD Tikhonov prior strength
--    REFINE_HOMOTOPY      (default FALSE)  coarse-to-fine |k| continuation
--    REFINE_HOMOTOPY_LAM  (default -1.0)   homotopy lam (<0 => REFINE_LAM)
--    REFINE_ROUNDS        (default 4)      alternating A/B/C rounds
--    REFINE_ITERS         (default 8)      GN iterations per round/stage
--    REFINE_STEP          (default 0.3)    initial lpap trust-region scale
--    REFINE_PER_LEG       (default TRUE)   Block B per-leg continuous GN
--    REFINE_LEG_STEP      (default 0.05)   initial Block B trust-region scale
--    REFINE_LEG_LAM       (default 1e-4)   Block B ridge damping
--    REFINE_HARM_SCREEN   (default TRUE)   Block C greedy 1-swap harmonics
--    REFINE_APPLY         (default TRUE)   write improved lt.exe.p[.k]
package GEM.LTE.Primitives.Refine is

   procedure Run (D : in out Shared.Param_S);

end GEM.LTE.Primitives.Refine;
