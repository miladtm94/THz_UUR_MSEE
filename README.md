# Minimum Secrecy Energy Efficiency Maximization for THz Untrusted UAV-Relaying

**Companion code for:**

> Milad Tatar Mamaghani et el., "[Terahertz Meets Untrusted UAV-Relaying: Minimum Secrecy Energy Efficiency Maximization via Trajectory and Communication Co-design](https://doi.org/10.1109/TVT.2022.3150011)," in *IEEE Transactions on Vehicular Technology*, vol. 71, no. 5, pp. 4991-5006, May 2022.



[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue)](https://www.mathworks.com/)
[![CVX](https://img.shields.io/badge/CVX-2.2-orange)](http://cvxr.com/cvx/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

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
| MATLAB | Used indirectly by CVX |

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

## Note on `.mat` Files

Figures 6–9 sweep simulation parameters across multiple values, each saved to a separate `.mat` file (e.g., `results/myResults_N100_af0.005.mat`). These files are **not included** in the repository (see `.gitignore`) because they can be fully regenerated by running the corresponding `run_figure_X` functions. 

---


## Citation

If you use this code, please cite:

```bibtex
@article{mamaghani2022terahertz,
  title={Terahertz meets untrusted UAV-relaying: Minimum secrecy energy efficiency maximization via trajectory and communication co-design},
  author={Tatar Mamaghani, Milad and Hong, Yi},
  journal={IEEE Transactions on Vehicular Technology},
  volume={71},
  number={5},
  pages={4991--5006},
  year={2022},
}
```

---

## License

This project is released under the MIT License. See `LICENSE` for details. 