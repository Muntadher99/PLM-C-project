# PLM-C: over-determined selective harmonic mitigation for three-level NPC inverters

MATLAB code, settings, seeds and results of the paper

> M. S. Sukkar and K. K. Abdalla, **"Over-Determined Selective Harmonic Mitigation: A Projected Levenberg-Marquardt Solver with Continuation for Three-Level NPC Inverters"**, under review.

PLM-C is an offline solver for selective harmonic mitigation (SHM) in three-level neutral-point-clamped (NPC) inverters. It computes the fifteen quarter-wave switching angles of a programmed-modulation pattern over the whole modulation range. The angles minimise low-order harmonic distortion, and a separate stage finds angles that keep every counted harmonic within its grid-code limit.

With fifteen angles, SHM is generically over-determined: the fundamental equality leaves fourteen degrees of freedom against sixteen counted harmonics (the odd, non-triplen orders 5 to 49), so in general no angle set cancels all of them. PLM-C therefore treats the design as constrained least squares. It works in gap coordinates, where the ordering conditions and the minimum pulse width define a simplex; a batched two-metric projected Levenberg-Marquardt corrector solves each operating point, and bidirectional continuation carries it across the modulation range.

The repository holds the solver, the five comparison methods, the matched-budget benchmark harness, the ablation, angle-count, compliance, device-level, loss and sensitivity studies, and every result file used in the paper. Every table of Sections 4 to 6 can be printed from the stored results, and every experiment can be rerun with the scripts, seeds and evaluation budgets that produced it. Section, table and equation numbers and bracketed references such as [2] refer to the paper.

---

## 1. Quick start

```matlab
cd PLM-C-project               % root of the repository
setup_paths                    % adds the code folders to the MATLAB path
reproduce_paper                % prints Tables 3-6, 8, 9 from the stored results (seconds)
reproduce_paper('all')         % reruns every experiment, then prints the tables again
reproduce_paper('witnesses')   % reruns one step
```

The steps are `benchmark`, `ablations`, `angle_count`, `witnesses`, `device`, `loss` and `sensitivity`, one per numbered folder (Section 7).

A rerun writes into `<folder>/reproduced/` under the same file names as the stored `<folder>/results/`, so the two can be compared file by file. `paper_tables('reproduced')` prints the tables from the reproduced files.

## 2. Requirements

| Software | Needed for |
| --- | --- |
| MATLAB R2024b (the release used for the paper) | everything |
| Optimization Toolbox (`fmincon`) | SQP-MS, SQP-C, GWO-SQP, MS-TB |
| Global Optimization Toolbox (`MultiStart`) | MS-TB |
| Simulink and Simscape Electrical | device-level validation (Section 6.2) |
| Parallel Computing Toolbox (`parsim`) | ideal-switch check (Section 6.2) |

PLM-C, PLM-MS, FSSA, the minimax witness search, the loss estimate and the sensitivity study use base MATLAB only. `paper_tables` also runs in GNU Octave. The method uses no trained model and no dataset.

## 3. Contents and where they appear in the paper

| Content | Folder | Paper |
| --- | --- | --- |
| PLM-C solver: harmonic model, gap coordinates, simplex projection, two-metric projected LM corrector, continuation | `1_proposed_method_PLM-C` | Sec. 3 |
| Comparison methods: GWO-SQP, MS-TB, FSSA (reimplementation of [2]) | `2_baselines` | Sec. 5.1, Table 2 |
| Matched-budget harness, common acceptance rule, settings of every method, NFE accounting and seed rules | `3_benchmark_protocol` | Sec. 4, Table 2 |
| Benchmark of the six methods at six budgets, ten seeds and 21 operating points | `3_benchmark_protocol` | Sec. 5.1-5.3, Tables 3, 4 |
| Achieved-NFE log of every benchmark run | `3_benchmark_protocol/results/nfe_log.csv` | Sec. 4 |
| FSSA with the weighted grid-code objective of [2] | `3_benchmark_protocol` (`bench_fssa_grid.csv`) | Sec. 5.6 |
| Local Gauss-Newton rate, corrector ablation at matched starts, corrector-by-strategy factorial | `4_ablations` | Sec. 2.2, 5.4, Tables 5, 6 |
| The same comparison at 11 and 19 angles | `5_angle_count` | Sec. 5.5 |
| Minimax compliance stage and every final angle set with its fundamental error, minimum gap, all normalised loadings q_n and max q_n | `6_compliance_witnesses` | Sec. 3.4, 5.7, Table 7 |
| Deployed design: distortion, compliance, switching and loss on one grid | `6_compliance_witnesses`, `8_loss_analysis` | Sec. 3.4, 5.7, Tables 7, 9 |
| Device-level Simulink validation with datasheet devices, dead time, DC-link imbalance and timer quantisation | `7_device_level_validation` | Sec. 6.2, Table 8 |
| Semiconductor loss estimate against carrier PWM | `8_loss_analysis` | Sec. 6.3, Table 9 |
| Sensitivity to dead time, neutral-point ripple, timer quantisation, modulation-index tolerance and minimum gap | `9_sensitivity` | Sec. 6.1 |

## 4. Package layout

Each numbered folder holds the code of one part of the study, a `run_*.m` script that reruns it with the settings of the paper, and a `results/` folder with the outputs used in the paper.

```
PLM-C-project/
  README.md
  setup_paths.m             add the code folders to the path
  reproduce_paper.m         single entry point (Section 1)
  paper_tables.m            Tables 3-6, 8, 9 and the statistics of Sections 4-6

  1_proposed_method_PLM-C/
    plm_core.m              Eq. (2), Jacobian, limits, gap coordinates, simplex projection,
                            batched two-metric projected LM corrector, NFE counter
    plm_sweep.m             sweep over the Ma grid: bidirectional continuation (PLM-C) or
                            multi-start (PLM-MS); PLM or fmincon-SQP corrector (SQP-C, SQP-MS)

  2_baselines/
    baseline_common.m       objective, constraints, SQP refinement and restoration shared by
                            GWO-SQP and MS-TB
    gwo_sqp_sweep.m         GWO-SQP
    mstb_sweep.m            MS-TB (MATLAB MultiStart)
    fssa_sweep.m            FSSA, reimplementation of the method of [2]

  3_benchmark_protocol/
    run_benchmark.m         matched-budget harness: arm registry, budgets, seeds, CSV writer
    evaluate_candidates.m   common restoration and acceptance rule applied to every arm
    summarize_benchmark.m   per-arm summary workbook of a benchmark CSV
    read_benchmark_csv.m, merge_benchmark_csv.m, csv_to_xlsx.m
    run_main_benchmark.m    -> results/bench_main.csv, bench_gwo_sqp.csv
    run_fssa_grid.m         FSSA with the objective of [2] (Section 5.6) -> results/bench_fssa_grid.csv
    export_nfe_log.m        achieved-NFE log of every benchmark CSV -> results/nfe_log.csv

  4_ablations/
    corrector_ablation.m    the local corrector at matched starts (Table 6)
    gauss_newton_rate.m     local rate rho of undamped Gauss-Newton (Section 2.2)
    run_ablations.m         -> results/bench_forward_continuation.csv, bench_sqp_continuation.csv,
                               corrector_ablation_production.csv, corrector_ablation_converged.csv,
                               gauss_newton_rate.csv

  5_angle_count/
    run_angle_count.m       -> results/bench_N11.csv, bench_N19.csv

  6_compliance_witnesses/
    minimax_witness_search.m   minimax compliance stage (Section 3.4)
    audit_angle_sets.m         independent audit of any angle set from Eq. (2) and Table 1
    build_witness_tables.m     selection of witnesses and designs, Table 7
    run_witnesses.m         -> results/minimax_search.csv, witnesses_all_audit.csv,
                               witnesses_selected.csv, thdmin_selected.csv, unified_table.csv

  7_device_level_validation/
    NPC_SHM.slx                 Simscape Electrical model of the NPC inverter
    device_parameters.m         IGBT and diode data from the datasheet [29]
    nonideal_simulink_runs.m    six cases x ten angle sets (Table 8)
    nonideal_switching_model.m  switching-function model of the same cases
    nonideal_model_init.m       sets the model variables for interactive use
    shmpwm3ph_15_deadtime.m     copy of the gate-logic block inside NPC_SHM.slx
    run_device_validation.m  -> results/nonideal_results.csv (input results/nonideal_cases.csv)
    ideal_switch_check/
      NPC_SHM_ideal.slx         the same inverter with ideal switches
      export_validation_cases.m angle sets of the six methods from the benchmark rows
      ideal_switch_validation.m simulation and THD49 comparison (script)

  8_loss_analysis/
    npc_loss_model.m        conduction and switching loss of the NPC inverter [30]
    loss_tradeoff_table.m   distortion, compliance, switching and loss on the grid
    run_loss_analysis.m     -> results/loss_table.csv (Table 9), loss_table_175C.csv

  9_sensitivity/
    sensitivity_study.m     dead time, neutral-point ripple, timer quantisation,
                            Ma interpolation, minimum gap (Section 6.1)
    sensitivity_summary.m   printed summary of a sensitivity run
    run_sensitivity.m       -> results/sensitivity*.csv/.xlsx, bench_mingap050.csv
```

## 5. Methods and the keys used in the result files

The `arm` column of every benchmark CSV carries a fixed key. The paper names are:

| Paper | `arm` key | Corrector | Start strategy | Code |
| --- | --- | --- | --- | --- |
| PLM-C (proposed) | `gn_hybrid` | projected two-metric LM | bidirectional continuation | `plm_sweep` |
| PLM-MS | `gn_multistart` | projected two-metric LM | multi-start | `plm_sweep` |
| SQP-MS | `sqp_multistart` | fmincon SQP | multi-start | `plm_sweep` |
| SQP-C | `sqp_hybrid` | fmincon SQP | bidirectional continuation | `plm_sweep` |
| Table 5, row 2 | `gn_hybrid_zeroth_fwd` | projected two-metric LM | forward continuation of [2] | `plm_sweep` |
| GWO-SQP | `gwo_sqp_box` | GWO, then fmincon SQP | box-initialised swarm [7] | `gwo_sqp_sweep` |
| MS-TB | `multistart_tb` | fmincon SQP | MATLAB MultiStart [28] | `mstb_sweep` |
| FSSA | `fssa` | simulated annealing | continuation in steps of 0.01 [2] | `fssa_sweep` |
| FSSA, objective of [2] (Sec. 5.6) | `fssa_grid` | simulated annealing | continuation in steps of 0.01 [2] | `fssa_sweep` (`OBJ = 'grid'`) |

The settings of every method (Table 2 of the paper) are the documented defaults at the top of each driver file and the overrides in the arm registry of `run_benchmark.m`.

## 6. Protocol

| Item | Setting |
| --- | --- |
| Operating points | Ma = 1.20, 1.15, ..., 0.20 (21 points), solved in this order |
| Angles, orders | N = 15 quarter-wave angles; THD49 over the 16 odd non-triplen orders 5 to 49 |
| Angle range, minimum gap | 0.5 deg <= theta_1 < ... < theta_15 <= 89 deg, consecutive angles at least 0.25 deg apart (0.50 deg in `bench_mingap050.csv`) |
| Budgets | 500, 1000, 2000, 4000, 8000 and 20000 NFE per operating point |
| Seeds | 1 to 10 for every benchmark; 1 to 5 for the minimax search |
| NFE | one per evaluation of Eq. (2) for one angle vector, Jacobian included; the external restoration of `evaluate_candidates` is not charged |
| Acceptance | after one common restoration onto b1 = Ma: \|b1 - Ma\| < 1e-8 and constraint violation < 1e-9 |
| Compliance | q_n = \|b_n\| / (L_n \|b1\|) <= 1 + 1e-4 for all 16 orders (limits of Table 1) |
| Exclusion | a comparison is not reported when the median achieved NFE of two methods differs by more than 15% |
| Statistics | per seed, the median THD49 over the points both methods accept; Hodges-Lehmann shift of the ten differences, exact two-sided Wilcoxon signed-rank p, Holm adjustment |

Seed rules (`s` is the seed, `i` the index of the operating point in solve order):

| Code | Random stream |
| --- | --- |
| PLM-C, PLM-MS, SQP-C, SQP-MS, forward continuation | `rng(1000*s + i, 'twister')`; reverse pass of the continuation: `1000*s + 500 + j` |
| GWO-SQP, MS-TB | `rng(s + i, 'twister')` |
| FSSA | `rng(s, 'twister')`, one stream per sweep |
| minimax witness search | `rng(1000*s + i)`, the same starts for its three objectives |
| corrector ablation | `20260930 + i`, the same 200 starts for every corrector |

Every random stream is set from these rules, so a rerun follows the same path: with MATLAB R2024b on the same hardware it returns the same angle sets and the same achieved NFE as the stored files. On other hardware, last-bit differences of the floating-point libraries can change individual angle sets. The wall-clock columns (`time_s`, `run_s`) depend on the machine.

## 7. Reproducing each table

| Paper item | Command | Reads | Writes |
| --- | --- | --- | --- |
| Tables 3, 4; Sec. 4, 5.2, 5.6 | `run_main_benchmark` | - | `bench_main.csv`, `bench_gwo_sqp.csv` |
| Sec. 5.6 (FSSA, objective of [2]) | `run_fssa_grid` | - | `bench_fssa_grid.csv` |
| achieved-NFE log | `export_nfe_log` (`export_nfe_log('reproduced')` after a rerun) | every benchmark CSV | `nfe_log.csv` |
| Table 5, Sec. 5.4 | `run_ablations` | `bench_main.csv` | `bench_forward_continuation.csv`, `bench_sqp_continuation.csv` |
| Table 6 | `run_ablations` | - | `corrector_ablation_production.csv`, `corrector_ablation_converged.csv` |
| Sec. 2.2 (rho) | `run_ablations` | `bench_main.csv` | `gauss_newton_rate.csv` |
| Sec. 5.5 | `run_angle_count` | - | `bench_N11.csv`, `bench_N19.csv` |
| Table 7, Sec. 5.7 | `run_witnesses` | `bench_main.csv` | `minimax_search.csv`, `witnesses_*.csv`, `thdmin_selected.csv`, `unified_table.csv`, `nonideal_cases.csv` |
| Table 8, Sec. 6.2 | `run_device_validation` | `nonideal_cases.csv`, benchmark CSVs | `nonideal_results.csv`, `validation_cases.xlsx`, `validation_cases_sim.xlsx` |
| Table 9, Sec. 6.3 | `run_loss_analysis` | `thdmin_selected.csv`, `witnesses_selected.csv` | `loss_table.csv`, `loss_table_175C.csv` |
| Sec. 6.1 | `run_sensitivity` | `bench_main.csv` | `sensitivity*.csv`, `sensitivity.xlsx`, `bench_mingap050.csv` |
| printed tables | `paper_tables` | all of the above | - |

Approximate run times on a desktop PC: main benchmark 1.4 h, ablations 0.5 h, angle counts 0.6 h, minimax search under 0.1 h, Simulink device study 0.7 h, sensitivity 1.3 h.

## 8. Result files

**Benchmark CSVs** (`bench_*.csv`, one row per method, budget, seed and operating point; `.xlsx` copy and `_summary.xlsx` alongside):
`arm, family, engine, strategy, use_spwm, nfe_target, nfe_achieved, nfe_ratio, seed, Ma, accepted, THD49_pct, THD99_pct, WTHD49_pct, WTHD999_pct, worst_harm_above49_pct, worst_harm_order, b3_pct, MAXQ, compliant, b1, b1_err, viol, n_over_limit, n_distinct_minima, time_s, theta01_deg ... thetaNN_deg`.
`nfe_achieved` is the achieved NFE per operating point: for PLM-C, PLM-MS, SQP-MS, SQP-C, the forward continuation and FSSA it is the evaluations charged over the whole sweep divided by its 21 points, and the same value is written on each of the 21 rows of the sweep; `time_s` is the wall time of the whole 21-point sweep. For GWO-SQP and MS-TB the evaluation counter is reset at every point and `nfe_achieved` is the count charged at the last point of the sweep.

**Achieved-NFE log** (`3_benchmark_protocol/results/nfe_log.csv`, written by `export_nfe_log` from the benchmark CSVs; one row per file, method, budget and seed):
`source, arm, method, nfe_target, seed, n_points, nfe_per_point, nfe_sweep_total, nfe_ratio, count_basis`.
`nfe_per_point` is the stored `nfe_achieved`; `nfe_sweep_total` is `n_points * nfe_per_point` where `count_basis` is `mean_over_points`, and is left empty where it is `last_point` (GWO-SQP, MS-TB).

**Audited angle sets** (`witnesses_all_audit.csv`, `witnesses_selected.csv`, `thdmin_selected.csv`; one row per angle set, recomputed from the angles alone by `audit_angle_sets`):
`label, Ma, theta01_deg..theta15_deg, b1, b1_err, THD49_pct, WTHD49_pct, q5..q49, q_max, binding_orders, min_gap_deg, first_margin_deg, last_margin_deg, viol_rad, feasible, compliant`.

**Table 7** (`unified_table.csv`): `Ma, THDmin_THD49_pct, THDmin_qmax, THDmin_compliant, seeds_compliant_of10, Qhat, Qhat_binding, design_THD49_pct, design_qmax, design_binding, dTHD49_pp, design_b1_err, design_min_gap_deg`.

The other files carry their column definitions in the header comment of the function that writes them.

## 9. Citation

The paper is under review. If you use this code or its results, please cite it; the full reference will be added here once it is published.
