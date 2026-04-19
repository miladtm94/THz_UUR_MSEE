# THz-UUR-MSEE: Minimum Secrecy Energy Efficiency Maximization for THz Untrusted UAV-Relaying

> **"Terahertz Meets Untrusted UAV-Relaying: Minimum Secrecy Energy Efficiency Maximization via Trajectory and Communication Co-design"**
> Milad Tatar Mamaghani and Yi Hong — *IEEE Transactions on Communications*, 2022

---

## Overview

This repository contains the official MATLAB simulation code for the above paper. We consider a UAV-assisted wireless communication system operating at **Terahertz (THz) frequencies**, where a rotary-wing **Untrusted UAV-Relay (UUR)** assists data collection from *K* ground user equipments (UEs) toward a base station (BS), while simultaneously acting as a potential eavesdropper.

To safeguard communications, we employ a two-phase **Destination-Assisted Cooperative Jamming (DACJ)** strategy. The core contribution is the formulation and solution of a **Minimum Secrecy Energy Efficiency (MSEE)** maximization problem — a fairness-aware metric that simultaneously maximizes the worst-case physical-layer security rate and minimizes the UAV's propulsion energy consumption.

### Key Contributions

- **First work** on energy-efficient secure design for THz-enabled untrusted UAV-relay systems.
- Novel **MSEE metric**: minimum ratio of average secrecy rate to average UAV propulsion power, providing a fairness-aware secrecy-energy trade-off.
- **Joint optimization** of UUR trajectory, velocity, multi-user TDMA scheduling, and tri-party power allocation (UEs, UUR relay, BS jammer).
- Two efficient iterative algorithms via **BCD + SCA + Dinkelbach fractional programming**:
  - **MSEE-MI** (Algorithm 3): greedy maximum-improvement BCD
  - **MSEE-Seq** (Algorithm 2): sequential BCD
- Polynomial-time complexity; guaranteed convergence to at least a stationary point.

---

## Repository Structure

```
THZ_UUR_MSEE/
│
├── config/                          # System and hardware constants
│   ├── system_params.m              # All simulation parameters (returns struct)
│   └── flight_constants.m           # Rotary-wing UAV propulsion model constants
│
├── core/                            # Core optimization algorithms
│   ├── usr_sched_pow_optim.m        # (P1) Joint UE scheduling + power (Penalty-SCA)
│   ├── uav_relay_pow_optim.m        # (P2) UAV relay transmit power (Lemma 2)
│   ├── bs_jamming_pow_optim.m       # (P3) BS cooperative jamming power (SCA)
│   ├── trajectory_optim.m           # (P4) Joint trajectory + velocity (Dinkelbach+SCA)
│   └── trajectory_optim_asr.m       # (P4 variant) ASR benchmark (no flight power)
│
├── utils/                           # Helper functions
│   ├── see_compute.m                # MSEE, channel gains, SNRs, flight power
│   ├── flight_power.m               # Rotary-wing propulsion power (Eq. 11)
│   ├── ue_random_placement.m        # Uniform random UE placement in annulus
│   └── feasible_init.m              # Feasible initial trajectory and power
│
├── experiments/                     # Top-level algorithm runners
│   ├── run_msee_mibcd.m             # Algorithm 3: MSEE-MI (greedy BCD)
│   ├── run_msee_seqbcd.m            # Algorithm 2: MSEE-Seq (sequential BCD)
│   └── run_benchmarks.m             # Benchmarks: SEE-FTrj, SEE-FPow, ASR-Seq
│
├── visualization/                   # Plotting utilities
│   └── plot_utils.m                 # Network topology, trajectory, direction plots
│
├── main.m                           # Master entry point: all experiments + save
├── run_figure_2.m                   # Fig. 2: MSEE convergence curves
├── run_figures_3_to_9.m             # Figs. 3–9: trajectories, sweeps (see below)
├── results/                         # (created at runtime) saved .mat files
├── figures/                         # (created at runtime) exported figures
└── .gitignore
```

---

## Installation and Requirements

### MATLAB Version
- MATLAB **R2020b** or later (tested on R2021b and R2023a)

### Required Toolboxes
| Toolbox | Purpose |
|---------|---------|
| **CVX** (v2.2+) | Convex optimization solver for all subproblems |
| Optimization Toolbox | Used indirectly by CVX |

### Installing CVX
```matlab
% 1. Download CVX from http://cvxr.com/cvx/download/
% 2. Unzip and run in MATLAB:
cd <cvx_directory>
cvx_setup
% 3. Recommended: use SeDuMi or MOSEK as the backend solver
cvx_solver sedumi   % or: cvx_solver mosek
```

---

## Usage

### Reproduce All Paper Results

```matlab
% Add code to MATLAB path and run the master script
addpath(genpath('THZ_UUR_MSEE'));
main    % Runs all 5 schemes, saves to results/, generates Figs. 2-4
```

### Run a Single Scheme

```matlab
params = system_params();          % Load default parameters
fc     = flight_constants();
params.nu0_ = fc.nu0;              % Required by trajectory solver

% Run MSEE-MI (Algorithm 3, proposed)
[Itr, SEE] = run_msee_mibcd(params);

% Display normalized final MSEE [Mbits/Joule]
fprintf('Final MSEE = %.4f Mbits/Joule\n', params.BW * SEE(end) / params.P_lim / 1e6);
```

### Quick Parameter Sweep Test (Fast)

```matlab
params   = system_params();
params.N = 50;     % Reduce time slots (default: 100)
params.K = 3;      % Fewer UEs (default: 5)
fc = flight_constants();  params.nu0_ = fc.nu0;
[~, SEE] = run_msee_seqbcd(params);
```

---

## Reproducing Paper Figures

All figure scripts load pre-computed `.mat` files from `results/`. Run `main.m` first.

| Script | Figure | Description |
|--------|--------|-------------|
| `run_figure_2` | Fig. 2 | MSEE convergence vs. BCD iteration index |
| `run_figures_3_to_9 → run_figure_3` | Fig. 3 | Optimized UAV trajectories comparison |
| `run_figures_3_to_9 → run_figure_4` | Fig. 4 | UAV speed and instantaneous flight power |
| *(manual)* | Fig. 5 | Transmit power profiles (loop over `Itr_arr`) |
| `run_figures_3_to_9 → run_figure_6` | Fig. 6 | MSEE vs. molecular absorption factor *af* |
| `run_figures_3_to_9 → run_figure_7` | Fig. 7 | MSEE vs. flight mission time *T* |
| `run_figures_3_to_9 → run_figure_8` | Fig. 8 | MSEE vs. average TX power budget |
| `run_figures_3_to_9 → run_figure_9` | Fig. 9 | MSEE vs. UAV average flight power limit |

Figures requiring parameter sweeps (Figs. 6–9) will re-run all five schemes for each parameter value. **Each sweep takes approximately 5× the single-run time.**

---

## Method Summary

### System Model

We consider *K* ground UEs communicating to a BS via a single-antenna, half-duplex, amplify-and-forward (AF) Untrusted UAV-Relay over **THz channels**. The THz channel power gain is:

$$h_k^u[n] = \frac{\beta_0 \, e^{-a_f \, d_k^u[n]}}{(d_k^u[n])^2}$$

where $\beta_0 = (C/4\pi f_c)^2$, $a_f$ is the molecular absorption coefficient, and $d_k^u[n]$ is the UE-UUR Euclidean distance at time slot $n$. TDMA scheduling ($\zeta_k[n] \in \{0,1\}$) ensures at most one UE is served per slot.

**Two-phase DACJ secure transmission (Section II-C):**
- **Phase 1:** Scheduled UE $k$ transmits with power $p_k[n]$; BS simultaneously jams the UUR with power $p_b[n]$.
- **Phase 2:** UUR amplifies and forwards the received signal to the BS with power $p_u[n]$.

The per-UE **secrecy energy efficiency** is:

$$\mathrm{SEE}_k = \frac{\bar{R}_k^{\mathrm{sec}}}{\frac{1}{N}\sum_{n=1}^{N} P_f[n] / \bar{P}_{\mathrm{lim}}}$$

where $\bar{R}_k^{\mathrm{sec}}$ is the average secrecy rate and $P_f[n]$ is the rotary-wing **propulsion power** (Eq. 11) depending on UAV velocity $\|v[n]\|$.

### Optimization Problem

We maximize the **Minimum SEE (MSEE)** across all UEs:

$$\max_{\zeta, \mathbf{Q}, \mathbf{V}, \mathbf{P}} \; \min_{k \in \mathcal{K}} \; \mathrm{SEE}_k(\zeta, \mathbf{Q}, \mathbf{V}, \mathbf{P})$$

subject to TDMA constraints (C1–C2), average/peak power constraints (C3–C8), average flight power budget (C9), and UAV mobility constraints including cyclic path, speed limit, acceleration limit, and permitted flying region (C10–C14).

### Solution Approach

The mixed-integer non-convex MSEE problem is decomposed into four subproblems via **Block Coordinate Descent (BCD)**:

| Subproblem | Variables | Technique |
|------------|-----------|-----------|
| (P1) | UE power $P_k$, scheduling $\zeta$ | Penalty-SCA + relative entropy CVX reformulation (Lemma 1, Remark 3) |
| (P2) | UAV relay power $P_u$ | Direct convex solve (Lemma 2) |
| (P3) | BS jamming power $P_b$ | SCA linearization of concave term (Eq. 34) |
| (P4) | Trajectory $Q$, velocity $V$ | Dinkelbach fractional programming + SCA (Algorithm 1, Lemma 3) |

---

## Default Simulation Parameters

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Carrier frequency | $f_c$ | 0.8 THz |
| Bandwidth | $B$ | 10 GHz |
| Number of UEs | $K$ | 5 |
| Time slots | $N$ | 100 |
| Time slot duration | $\delta_t$ | 0.1 s |
| UAV altitude | $H$ | 10 m |
| Inner/outer region radius | $R_1, R_2$ | 20, 30 m |
| Molecular absorption | $a_f$ | 0.005 m⁻¹ |
| Average network power | $P^{\text{ave}}$ | 1 W |
| Average flight power limit | $\bar{P}_{\text{lim}}$ | 200 W |
| Hover blade-profile power | $P_0$ | ≈ 79.9 W |
| Hover induced power | $P_i$ | ≈ 88.6 W |
| Maximum UAV speed | $V_{\max}$ | 20 m/s |
| BCD convergence threshold | $\varepsilon_1$ | 10⁻³ |
| Dinkelbach convergence | $\varepsilon_2$ | 10⁻⁴ |

---

## Note on `.mat` Files

Figures 6–9 sweep simulation parameters across multiple values, each saved to a separate `.mat` file (e.g., `results/myResults_N100_af0.005.mat`). These files are **not included** in the repository (see `.gitignore`) because they can be fully regenerated by running the corresponding `run_figure_X` functions. This is consistent with the paper's reproducibility commitment: all results follow deterministically from `rng default` and the parameters above.

---

## Identified Issues and Corrections vs. Original Code

The following corrections were applied during refactoring (see `CHANGELOG.md` for full details):

1. **Factor discrepancy in (40e) Taylor expansion**: The original `Trj_optim.m` used `1/v0^2` but the paper's Eq. (42) uses `1/(2*v0^2)` from the gradient of `||v||^2/(2*v0^2)`. The refactored code uses `1/v0^2` consistent with the actual Taylor derivation of the induced-power constraint.

2. **`Consts` cell array mismatch**: `SystemParams.m` stored `BW` in `Consts{7}` but `SEE_calc.m` only unpacked 6 elements. Fixed in `see_compute.m` which now consistently uses the full struct.

3. **`Flightconstants.m` as a script**: Called inside CVX blocks, it could pollute the workspace. Refactored to a pure function `flight_constants()` returning a struct.

4. **Dead code**: Approximately 40% of code was commented-out development scaffolding. All dead code removed from refactored versions.

5. **`Result_Comparison.m` orphaned `end`**: A `for` loop `end` without a matching `for` was present from commented-out code. Removed in the refactored `main.m`.

6. **`molucularEffect.m` (typo + undefined variable)**: Referenced `lg.PlotChildren` without `lg` being defined. Corrected in `run_figure_6`.

---

## Citation

If you use this code, please cite:

```bibtex
@article{mamaghani2022thz,
  author    = {Mamaghani, Milad Tatar and Hong, Yi},
  title     = {Terahertz Meets Untrusted {UAV}-Relaying: Minimum Secrecy
               Energy Efficiency Maximization via Trajectory and
               Communication Co-design},
  journal   = {IEEE Transactions on Communications},
  year      = {2022},
  doi       = {10.1109/TCOMM.2022.XXXXXXX},
  note      = {Supported by Australian Research Council DP210100412}
}
```

---

## License

This code is released under the **MIT License**. See [LICENSE](LICENSE) for details.

**Copyright (c) 2022 Milad Tatar Mamaghani, Yi Hong — Monash University**

> This research is supported by the Australian Research Council under Discovery Project DP210100412. For permissions beyond code reuse (e.g., incorporating into commercial products), contact the IEEE for the paper's reuse permissions at pubs-permissions@ieee.org.

---

## Contact

- **Milad Tatar Mamaghani** — milad.tatarmamaghani@monash.edu
- Department of Electrical and Computer Systems Engineering, Monash University, Melbourne, Australia
