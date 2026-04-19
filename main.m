% MAIN  Master entry point: run all experiments and save results.
%
%   Run this script to reproduce all numerical results from:
%
%   M. T. Mamaghani and Y. Hong, "Terahertz Meets Untrusted UAV-Relaying:
%   Minimum Secrecy Energy Efficiency Maximization via Trajectory and
%   Communication Co-design," IEEE Transactions on Communications, 2022.
%
%   WORKFLOW:
%   1. Loads system parameters (default N=100, P_lim=200 W, af=0.005).
%   2. Runs all five schemes:
%       - MSEE-MI  (Algorithm 3, proposed)
%       - MSEE-Seq (Algorithm 2, proposed)
%       - SEE-FTrj (fixed trajectory benchmark)
%       - SEE-FPow (fixed power benchmark)
%       - ASR-Seq  (ASR maximization benchmark)
%   3. Saves results to results/myResults_N100_default.mat
%   4. Produces all paper figures (Figs. 2-9) via run_figure_X.m scripts.
%
%   REQUIREMENTS:
%       - MATLAB R2020b or later
%       - CVX toolbox (http://cvxr.com/cvx/) with SeDuMi or MOSEK solver
%
%   ESTIMATED RUNTIME: ~4-8 hours on a modern workstation (N=100, K=5).
%   For a quick test, set params.N = 50 in system_params.m.
%
%   NOTE ON SAVED .mat FILES:
%   Figures 6-9 require sweeping parameters (af, T, P_ave, P_lim).
%   Each sweep point saves its own .mat file. Run run_figure_X.m scripts
%   individually after main.m to generate the corresponding figures.

clc; clear; close all;

%% ======================================================================
%  Add all subdirectories to MATLAB path
% =======================================================================
addpath(genpath(fileparts(mfilename('fullpath'))));
if ~exist('results', 'dir'), mkdir('results'); end

%% ======================================================================
%  Load system parameters (default configuration)
% =======================================================================
params = system_params();
N      = params.N;

% Add nu0 field needed by run functions (flight constant)
fc         = flight_constants();
params.nu0_ = fc.nu0;

fprintf('=== THz-UUR MSEE Simulation ===\n');
fprintf('Parameters: N=%d, K=%d, P_lim=%d W, af=%.4f\n', ...
    N, params.K, params.P_lim, params.af);
fprintf('Expected runtime: several hours. Set N=50 for a quick test.\n\n');

filename = sprintf('results/myResults_N%d_default.mat', N);

%% ======================================================================
%  Scheme 1: MSEE-MI (Algorithm 3 — proposed, greedy BCD)
% =======================================================================
fprintf('\n====== MSEE-MI (Algorithm 3) ======\n');
[Itr_MIBCD, SEE_MIBCD] = run_msee_mibcd(params);
save(filename, 'Itr_MIBCD', 'SEE_MIBCD');
fprintf('MSEE-MI complete. Final MSEE = %.4f Mbits/Joule\n', ...
    params.BW * SEE_MIBCD(end) / params.P_lim / 1e6);

%% ======================================================================
%  Scheme 2: MSEE-Seq (Algorithm 2 — proposed, sequential BCD)
% =======================================================================
fprintf('\n====== MSEE-Seq (Algorithm 2) ======\n');
[Itr_SeqBCD, SEE_SeqBCD] = run_msee_seqbcd(params);
save(filename, 'Itr_SeqBCD', 'SEE_SeqBCD', '-append');
fprintf('MSEE-Seq complete. Final MSEE = %.4f Mbits/Joule\n', ...
    params.BW * SEE_SeqBCD(end) / params.P_lim / 1e6);

%% ======================================================================
%  Scheme 3: SEE-FTrj (fixed trajectory benchmark)
% =======================================================================
fprintf('\n====== SEE-FTrj (fixed-trajectory benchmark) ======\n');
[Itr_FixedTrj, SEE_FixedTrj] = run_msee_fixed_traj(params);
save(filename, 'Itr_FixedTrj', 'SEE_FixedTrj', '-append');

%% ======================================================================
%  Scheme 4: SEE-FPow (fixed power benchmark)
% =======================================================================
fprintf('\n====== SEE-FPow (fixed-power benchmark) ======\n');
[Itr_FixedPow, SEE_FixedPow] = run_msee_fixed_pow(params);
save(filename, 'Itr_FixedPow', 'SEE_FixedPow', '-append');

%% ======================================================================
%  Scheme 5: ASR-Seq (ASR maximization benchmark)
% =======================================================================
fprintf('\n====== ASR-Seq (ASR benchmark) ======\n');
[Itr_ASRoptim, SEE_ASRoptim] = run_asr_seq(params);
save(filename, 'Itr_ASRoptim', 'SEE_ASRoptim', '-append');

%% ======================================================================
%  Summary
% =======================================================================
fprintf('\n===== RESULTS SUMMARY (N=%d) =====\n', N);
schemes = {'MSEE-MI', 'MSEE-Seq', 'SEE-FTrj', 'SEE-FPow', 'ASR-Seq'};
SEEs    = {SEE_MIBCD, SEE_SeqBCD, SEE_FixedTrj, SEE_FixedPow, SEE_ASRoptim};
for i = 1:5
    s = SEEs{i};
    fprintf('  %s: %.4f Mbits/Joule\n', schemes{i}, ...
        params.BW * s(end) / params.P_lim / 1e6);
end
fprintf('\nAll results saved to: %s\n', filename);

%% ======================================================================
%  Generate paper figures
% =======================================================================
fprintf('\nGenerating paper figures...\n');
run_figure_2;   % Convergence: MSEE vs. iteration index
run_figure_3;   % UAV trajectories comparison
run_figure_4;   % Velocity and instantaneous flight power
% run_figure_5; % Transmission power profiles (requires manual inspection)
fprintf('Figures 2-4 generated. Run run_figure_6 through run_figure_9\n');
fprintf('separately after running the parameter sweep scripts.\n');
