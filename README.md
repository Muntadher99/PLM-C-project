# Over-Determined Selective Harmonic Mitigation: A Projected Levenberg-Marquardt Solver with Continuation for Three-Level NPC Inverters

**Muntadher S. Sukkar**, **Kasim K. Abdalla**
Department of Electrical Engineering, College of Engineering, University of Babylon, Hilla, Iraq

> **Paper under review.** Code, experiment scripts, results and Simulink models will be released upon acceptance.

## Overview

PLM-C is an offline solver for selective harmonic mitigation (SHM) in three-level neutral-point-clamped (NPC) inverters. It computes the fifteen quarter-wave switching angles of a programmed-modulation pattern over the whole modulation range. The angles minimise low-order harmonic distortion, and a separate stage finds angles that keep every counted harmonic within its grid-code limit.

With fifteen angles, SHM is generically over-determined. The fundamental equality leaves fourteen degrees of freedom against sixteen counted harmonics (the odd, non-triplen orders 5 to 49), so in general no angle set cancels all of them. The design is therefore a constrained least-squares problem rather than a root-finding one. Existing solvers either keep the square elimination formulation, which presupposes a root that does not exist at this angle count, or treat the design as a black box for a population metaheuristic. PLM-C instead uses the structure of the problem. It works in gap coordinates, where the ordering conditions and the minimum pulse width define a simplex. A batched two-metric projected Levenberg-Marquardt corrector solves each operating point, and bidirectional continuation carries it across the modulation range.

## Key Features

* **Gap-coordinate formulation** — the ordering conditions and the minimum pulse width become a simplex, so every iterate is admissible by construction and projection takes a single sort
* **Two-metric projected Levenberg-Marquardt corrector** — damped Gauss-Newton metric on the free gaps and identity on the gaps held at their bound, batched so that many starts advance together
* **Bidirectional continuation** — a forward pass carries solutions from one modulation index to the next, and a reverse pass reuses an archive pooled over the whole range; no lookup table or trained initialiser is needed
* **Minimax compliance stage** — Lawson reweighting on the same corrector finds a grid-code-compliant angle set at every operating point
* **Matched-evaluation benchmark** — six methods compared under the same per-point evaluation budget, one common acceptance rule and ten paired seeds, with ablations that isolate each component
* **Device-level and loss validation** — Simulink simulation with datasheet devices, dead time, DC-link imbalance and timer quantisation, and a semiconductor loss estimate against carrier PWM

## Code

Coming soon. The release will include:

* the MATLAB implementation of PLM-C and of every baseline method, with their settings
* the benchmark harness, the seed lists and the scripts that regenerate every table of the paper
* the per-seed results, including the number of function evaluations used in every run
* every final angle set, with its fundamental error, minimum gap and normalised harmonic loadings
* the Simulink models of the device-level validation

If you want to be notified, watch this repository.

## Citation

The citation will be added upon publication.
