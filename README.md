# AA amyloidosis: SAA–HDL binding and renal biological age

MATLAB source code for the manuscript *Mathematical Modeling of AA
Amyloidosis: Coupling SAA–HDL Binding Dynamics with Path-Dependent Renal
Aging* (A. V. Kuznetsov). Base MATLAB only; developed and tested with R2023a.

Put all files in the same folder and run the scripts from that folder. Every
script is self-contained: open it and press Run. Each figure window it opens is
named after the figure and panel of the paper that it reproduces (for example
'Fig. 2a'), and the numerical values quoted in the text are printed in the
Command Window.

For Figs. 2–7 and S1–S26, panel (a) uses biologically relevant half-lives
(functional protein degradation) and panel (b) infinite half-lives (impaired
protein degradation). One run of the script produces both panels.

## Figures and tables → scripts

| figure | script | content |
|---|---|---|
| Fig. 1 | `AA_Kp.m` | SAA–HDL binding fractions |
| Fig. 2 | `AA_Htot.m` | combined damage, three [H]tot values |
| Fig. 3 | `AA_Htot.m` | combined damage: fibril removal and linearly increasing [B] |
| Fig. 4 | `AA_Htot.m` | biological age, three [H]tot values |
| Fig. 5 | `AA_Htot_sensitivity.m` | sensitivity to [H]tot |
| Fig. 6 | `AA_kcl.m` | biological age, three kcl values |
| Fig. 7 | `AA_kcl_sensitivity.m` | sensitivity to kcl |
| Fig. 8 | `AA_tipping_point.m` | autocatalytic threshold |
| Fig. 9 | `AA_flares.m` | intermittent flares vs sustained inflammation: (a) prescribed SAA; (b) free SAA; (c,d) renal oligomer concentration; (e,f) combined damage |
| Fig. 10 | `AA_treatment.m` | suppression of SAA at 1, 3, 5 years |
| Fig. 11 | `AA_validation.m` | comparison with clinical data |
| Table 3 | `*_sensitivity.m` | value at 6 years for the baseline case, printed for each panel; the impaired-degradation column is also printed by `AA_calibration_c.m` |
| Figs. S1–S6 | `AA_Htot.m` | [S], [A], [B], [I], fibril volume fraction, accumulated toxicity |
| Fig. S7 | `AA_Htot.m` | [B]: model vs linear increase |
| Figs. S8–S14 | `AA_kcl.m` | [S], [A], [B], [I], fibril volume fraction, accumulated toxicity, combined damage |
| Fig. S15 | `AA_hs.m` | biological age, three hS values |
| Fig. S16 | `AA_hs_sensitivity.m` | sensitivity to hS |
| Fig. S17 | `AA_theta12B.m` | biological age, three θ1/2,B values |
| Fig. S18 | `AA_theta12B_sensitivity.m` | sensitivity to θ1/2,B |
| Fig. S19 | `AA_Kd.m` | biological age, three Kd values |
| Fig. S20 | `AA_Kd_sensitivity.m` | sensitivity to Kd |
| Fig. S21 | `AA_Stot.m` | biological age, three Stot values |
| Fig. S22 | `AA_Stot_sensitivity.m` | sensitivity to Stot |
| Fig. S23 | `AA_k1.m` | biological age, three k1 values |
| Fig. S24 | `AA_k1_sensitivity.m` | sensitivity to k1 |
| Fig. S25 | `AA_k2.m` | biological age, three k2 values |
| Fig. S26 | `AA_k2_sensitivity.m` | sensitivity to k2 |
| Fig. S27 | `AA_flares.m` | early vs late inflammatory episode; (a) prescribed SAA |
| Fig. S28 | `AA_calibration_c.m` | calibration of c |

## Scripts

**`AA_Kp.m`** — Fig. 1. SAA–HDL binding equilibrium, Section 2.1.1. Standalone.

**`AA_Htot.m`, `AA_kcl.m`, `AA_hs.m`, `AA_theta12B.m`, `AA_Kd.m`, `AA_Stot.m`,
`AA_k1.m`, `AA_k2.m`** — Figs. 2–4, 6, S1–S15, S17, S19, S21, S23 and S25 (see
table above). Solve Eqs. (15), (19), (20) and (26) with initial conditions (32)
for three values of the parameter named in the file name, and plot the
numerical (`ode45`) and analytical solutions together. Damage measures from
Eqs. (45) and (47), biological age from Eq. (49). Require `funCDCMCB_AA.m`.

**`AA_<parameter>_sensitivity.m`** — Figs. 5, 7, S16, S18, S20, S22, S24 and
S26 (see table above). Normalized local sensitivity of biological age,
Eq. (53), computed by forward differences, Eq. (52), with a relative step of
1e-3. The accumulated toxicity, Eq. (45), is integrated together with the ODEs
so that the start-up transient is resolved. Require `funCDCMCB_AA.m`.

**`AA_tipping_point.m`** — Fig. 8. Dimensionless reduction of the aggregation
step, Eqs. (54)–(58); locates the transcritical threshold κ2 = 1 and the
critical HDL concentration for each k2, and verifies the closed form against
`ode45`. Also checks the linear stability analysis of Section 2.5: prints the
Jacobian eigenvalues of the monomer–oligomer subsystem on both steady-state
branches (κ1 = 0 and baseline κ1). Requires `AA_params.m`, `funCDCMCB_AA.m`,
`savefig_emf.m`.

**`AA_flares.m`** — Figs. 9 and S27. Intermittent inflammatory flares versus
sustained inflammation at the same time-averaged SAA concentration, and early
versus late inflammatory episodes of identical magnitude; panels (a) show the
prescribed total SAA concentration, and Figs. 9b–d the resulting free SAA and
renal oligomer concentrations. Requires `AA_params.m`, `funCDCMCB_AA_tv.m`,
`funAAsolve.m`, `funflare.m`, `savefig_emf.m`.

**`AA_treatment.m`** — Fig. 10. Suppression of circulating SAA below 10 mg/L
started at one, three and five years after onset. Requires `AA_params.m`,
`funCDCMCB_AA_tv.m`, `funAAsolve.m`, `savefig_emf.m`.

**`AA_validation.m`** — Fig. 11. Predicted time to the calibrated threshold
BA = 100 years and rate of renal biological ageing as functions of median
circulating SAA, compared with Gillmore et al. (2001) and Lachmann et al.
(2007). Uses the analytical steady state and verifies it against `ode45`.
Requires `AA_params.m`, `funAAsteady.m`, `funAAdamage.m`, `funAAsolve.m`,
`funCDCMCB_AA_tv.m`, `savefig_emf.m`.

**`AA_calibration_c.m`** — Fig. S28. Value of c implied by Eq. (51) over the
reported 4–8 year survival range, and the effect of that choice on the
dimensionless sensitivity coefficients. Requires `AA_params.m`,
`funCDCMCB_AA_tv.m`, `funAAsolve.m`, `savefig_emf.m`.

The scripts that call `savefig_emf.m` also write each figure to an `.emf` file
(for example `fig8a.emf`) in the current folder.

## Shared routines

| file | purpose |
|---|---|
| `funCDCMCB_AA.m` | right-hand side of Eqs. (15), (19), (20) and (26) |
| `funCDCMCB_AA_tv.m` | as above, but for a prescribed time-dependent free SAA concentration; used by the time-varying simulations |
| `AA_params.m` | baseline parameter values of Table 2 |
| `funAAsolve.m` | integrates the system for a prescribed total SAA history and returns the damage measures of Sections 2.2 and 2.3: the combined burden Qtot, Eq. (47); biological age formed from it, Eq. (49), which falls when fibrillar deposits regress; and the irreversible floor of biological age, Eq. (50), formed from the accumulated oligomer toxicity Q alone |
| `funAAsteady.m` | closed-form steady state written in terms of the dimensionless groups κ1 and κ2 (Section 2.5) |
| `funAAdamage.m` | damage measures from analytical solutions (27), (41), (46) |
| `funflare.m` | smoothed periodic inflammatory flare profile |
| `savefig_emf.m` | figure export helper |

## License

MIT; see `LICENSE`.
