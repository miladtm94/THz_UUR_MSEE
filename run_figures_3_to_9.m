% RUN_FIGURE_3  Reproduce Figure 3: Optimized UAV trajectories.
%
%   Plots the final converged UAV trajectory for each of the five schemes
%   overlaid on the network topology (annulus, BS, UEs).
%   Corresponds to Fig. 3 in the paper (Section IV).

function run_figure_3()
addpath(genpath(fileparts(mfilename('fullpath'))));
params = system_params();
fc = flight_constants();
params.nu0_ = fc.nu0;

filename = sprintf('results/myResults_N%d_default.mat', params.N);
if ~exist(filename, 'file')
    error('run_figure_3: result file not found. Run main.m first.');
end
load(filename, 'Itr_MIBCD', 'Itr_SeqBCD', 'Itr_FixedTrj', 'Itr_FixedPow', 'Itr_ASRoptim');

Itr_arr = {Itr_MIBCD, Itr_SeqBCD, Itr_FixedTrj, Itr_FixedPow, Itr_ASRoptim};

[~, ax_traj, lgd_traj] = plot_network_topology(params);

for i = 1:length(Itr_arr)
    Itr  = Itr_arr{i};
    Qu   = Itr(end).Trj';   % [3xN]
    plot_uav_trajectory(Qu(1,:), Qu(2,:), ax_traj, lgd_traj, i, '');
end
grid(ax_traj, 'on');

if ~exist('figures', 'dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig3_trajectories.fig');
saveas(gcf, 'figures/fig3_trajectories.png');
fprintf('Figure 3 saved.\n');
end

% =========================================================================

% RUN_FIGURE_4  Reproduce Figure 4: Velocity and instantaneous flight power.
%
%   Two-subplot figure: top panel shows UAV speed [m/s] vs. time slot;
%   bottom panel shows instantaneous flight power consumption [W] vs. time.
%   Corresponds to Fig. 4 in the paper (Section IV).

function run_figure_4()
addpath(genpath(fileparts(mfilename('fullpath'))));
params = system_params();

filename = sprintf('results/myResults_N%d_default.mat', params.N);
if ~exist(filename, 'file')
    error('run_figure_4: result file not found. Run main.m first.');
end
load(filename, 'Itr_MIBCD', 'Itr_SeqBCD', 'Itr_FixedTrj', 'Itr_FixedPow', 'Itr_ASRoptim');

Itr_arr   = {Itr_MIBCD, Itr_SeqBCD, Itr_FixedTrj, Itr_FixedPow, Itr_ASRoptim};
lgd_names = {'SEE-MI', 'SEE-Seq', 'SEE-FTrj', 'SEE-FPow', 'ASR-Seq'};
colours   = {'#0072BD', '#D95319', '#7E2F8E', '#A2142F', '#77AC30'};
markers   = {'s', '^', 'p', 'd', 'o'};
linestyle = {'-', '--', ':', '-.', '-'};

T  = params.N * params.dt;
tt = linspace(0, T, params.N);

figure;
for i = 1:length(Itr_arr)
    Itr = Itr_arr{i};
    subplot(2,1,1); hold on;
    pp = plot(tt, vecnorm(Itr(end).Vel', 2));
    pp.LineStyle     = linestyle{i};  pp.Marker = markers{i};
    pp.MarkerIndices = i : 10-i : params.N;
    pp.LineWidth     = 2;

    subplot(2,1,2); hold on;
    qq = plot(tt, Itr(end).prop);
    qq.LineStyle     = linestyle{i};  qq.Marker = markers{i};
    qq.MarkerIndices = i : 10-i : params.N;
    qq.LineWidth     = 2;
end

subplot(2,1,1);
ylabel("UAV speed [m/s]", 'FontName', 'Times New Roman');
legend(lgd_names, 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

subplot(2,1,2);
ylabel("Instantaneous flight power [W]", 'FontName', 'Times New Roman');
xlabel('Time [s]', 'FontName', 'Times New Roman');
legend(lgd_names, 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
xlim([0, T]);
grid on;

if ~exist('figures', 'dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig4_velocity_flightpower.fig');
saveas(gcf, 'figures/fig4_velocity_flightpower.png');
fprintf('Figure 4 saved.\n');
end

% =========================================================================

% RUN_FIGURE_6  Reproduce Figure 6: MSEE vs. molecular absorption coefficient.
%
%   Sweeps af in [0.005, 0.025] (5 values) for all schemes.
%   Each sweep point is saved as results/myResults_N100_af<val>.mat.
%   Corresponds to Fig. 6 in the paper (Section IV).
%
%   ESTIMATED RUNTIME: ~5x longer than a single default run.

function run_figure_6()
addpath(genpath(fileparts(mfilename('fullpath'))));
if ~exist('results', 'dir'), mkdir('results'); end

af_vec = linspace(0.005, 0.025, 5);
N      = 100;
myRes  = zeros(length(af_vec), 5);

for nnn = 1:length(af_vec)
    params      = system_params();
    params.N    = N;
    params.af   = af_vec(nnn);
    params.Consts = {params.P_lim, params.Qa, params.Qb, ...
                     params.af, params.beta0, params.N0, params.BW};
    fc = flight_constants();  params.nu0_ = fc.nu0;

    fname = sprintf('results/myResults_N%d_af%.3f.mat', N, af_vec(nnn));

    [~, SEE_MI]  = run_msee_mibcd(params);     SEE_MIBCD   = SEE_MI;
    [~, SEE_Seq] = run_msee_seqbcd(params);    SEE_SeqBCD  = SEE_Seq;
    [~, SEE_FT]  = run_msee_fixed_traj(params);SEE_FixedTrj= SEE_FT;
    [~, SEE_FP]  = run_msee_fixed_pow(params); SEE_FixedPow= SEE_FP;
    [~, SEE_ASR] = run_asr_seq(params);        SEE_ASRoptim= SEE_ASR;

    save(fname, 'SEE_MIBCD', 'SEE_SeqBCD', 'SEE_FixedTrj', 'SEE_FixedPow', 'SEE_ASRoptim');

    SEEs = {SEE_MI, SEE_Seq, SEE_FT, SEE_FP, SEE_ASR};
    for i = 1:5
        myRes(nnn, i) = params.BW * SEEs{i}(end) / params.P_lim / 1e6;
    end
    fprintf('af = %.4f done.\n', af_vec(nnn));
end

% --- Plot ---
colours   = {'k', '#0072BD', '#D95319', '#7E2F8E', '#77AC30'};
markers   = {'none', '^', '*', '+', 'x'};
linestyle = {'-', '--', ':', '-.'};
lgd_names = {'MSEE-MI', 'MSEE-Seq', 'SEE-FTrj', 'SEE-FPow', 'ASR-Seq'};

figure; hold on;
for i = 1:5
    p = plot(af_vec, myRes(:, i), 'LineWidth', 2);
    p.LineStyle = linestyle{mod(i,4)+1};
    p.Marker    = markers{i};
    p.Color     = colours{i};
end
xlabel('Molecular absorption coefficient [m^{-1}]', 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Min. secrecy energy efficiency [Mbits/Joule]', 'FontName', 'Times New Roman', 'FontSize', 11);
lgd = legend(lgd_names);
lgd.Location = 'best';  lgd.FontSize = 10;  lgd.FontName = 'Times New Roman';
grid on;

if ~exist('figures','dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig6_vs_af.fig');
saveas(gcf, 'figures/fig6_vs_af.png');
fprintf('Figure 6 saved.\n');
end

% =========================================================================

% RUN_FIGURE_7  Reproduce Figure 7: MSEE vs. flight mission time T.
%
%   Sweeps N in [50, 70, 80, 90, 100, 120] (T = N * dt seconds).
%   Corresponds to Fig. 7 in the paper (Section IV).

function run_figure_7()
addpath(genpath(fileparts(mfilename('fullpath'))));
if ~exist('results', 'dir'), mkdir('results'); end

Nvec   = [50, 70, 80, 90, 100, 120];
dt     = system_params().dt;
myRes  = zeros(length(Nvec), 5);
myInit = zeros(length(Nvec), 5);

for nnn = 1:length(Nvec)
    params   = system_params();
    params.N = Nvec(nnn);
    fc = flight_constants();  params.nu0_ = fc.nu0;

    fname = sprintf('results/myResults_N%d.mat', Nvec(nnn));

    [~, SEE_MI]  = run_msee_mibcd(params);
    [~, SEE_Seq] = run_msee_seqbcd(params);
    [~, SEE_FT]  = run_msee_fixed_traj(params);
    [~, SEE_FP]  = run_msee_fixed_pow(params);
    [~, SEE_ASR] = run_asr_seq(params);

    SEE_MIBCD = SEE_MI; SEE_SeqBCD = SEE_Seq;
    SEE_FixedTrj = SEE_FT; SEE_FixedPow = SEE_FP; SEE_ASRoptim = SEE_ASR;
    save(fname, 'SEE_MIBCD', 'SEE_SeqBCD', 'SEE_FixedTrj', 'SEE_FixedPow', 'SEE_ASRoptim');

    SEEs = {SEE_MI, SEE_Seq, SEE_FT, SEE_FP, SEE_ASR};
    for i = 1:5
        myRes(nnn,i)  = params.BW * SEEs{i}(end) / params.P_lim / 1e6;
        myInit(nnn,i) = params.BW * SEEs{i}(1)   / params.P_lim / 1e6;
    end
    fprintf('N = %d done.\n', Nvec(nnn));
end

T_vec = Nvec * dt;
colours   = {'k', '#0072BD', '#D95319', '#7E2F8E', '#77AC30'};
markers   = {'p', '^', '*', '+', 'x'};
linestyle = {'-', '--', ':', '-.'};
lgd_names = {'MSEE-MI', 'MSEE-Seq', 'SEE-FTrj', 'SEE-FPow', 'ASR-Seq', 'Initial'};

figure; hold on;
for i = 1:5
    p = plot(T_vec, myRes(:,i), 'LineWidth', 2);
    p.LineStyle = linestyle{mod(i,4)+1};
    p.Marker    = markers{i};
    p.Color     = colours{i};
end
p_init = plot(T_vec, myInit(:,1), ':', 'Color', 'm', 'LineWidth', 2);

xlabel('Flight mission time T [s]', 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Min. secrecy energy efficiency [Mbits/Joule]', 'FontName', 'Times New Roman', 'FontSize', 11);
legend(lgd_names, 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

if ~exist('figures','dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig7_vs_T.fig');
saveas(gcf, 'figures/fig7_vs_T.png');
fprintf('Figure 7 saved.\n');
end

% =========================================================================

% RUN_FIGURE_8  Reproduce Figure 8: MSEE vs. average network transmit power.
%
%   Sweeps P_ave in {1, 2, 4, 8} W.
%   Corresponds to Fig. 8 in the paper (Section IV).

function run_figure_8()
addpath(genpath(fileparts(mfilename('fullpath'))));
if ~exist('results', 'dir'), mkdir('results'); end

Pave_vec = [1, 2, 4, 8];
N        = 100;
myRes    = zeros(length(Pave_vec), 5);

for nnn = 1:length(Pave_vec)
    params        = system_params();
    params.N      = N;
    params.P_ave  = Pave_vec(nnn);
    % Update all power budgets derived from P_ave
    params.Pa_ave  = 0.1 * params.P_ave;
    params.Pa_peak = 4   * params.Pa_ave;
    params.Pb_ave  = 0.5 * params.P_ave;
    params.Pb_peak = 4   * params.Pb_ave;
    params.Pu_ave  = 0.4 * params.P_ave;
    params.Pu_peak = 4   * params.Pu_ave;
    fc = flight_constants();  params.nu0_ = fc.nu0;

    fname = sprintf('results/myResults_N%d_Pave%d.mat', N, Pave_vec(nnn));

    [~, SEE_MI]  = run_msee_mibcd(params);
    [~, SEE_Seq] = run_msee_seqbcd(params);
    [~, SEE_FT]  = run_msee_fixed_traj(params);
    [~, SEE_FP]  = run_msee_fixed_pow(params);
    [~, SEE_ASR] = run_asr_seq(params);

    SEE_MIBCD = SEE_MI; SEE_SeqBCD = SEE_Seq;
    SEE_FixedTrj = SEE_FT; SEE_FixedPow = SEE_FP; SEE_ASRoptim = SEE_ASR;
    save(fname, 'SEE_MIBCD', 'SEE_SeqBCD', 'SEE_FixedTrj', 'SEE_FixedPow', 'SEE_ASRoptim');

    SEEs = {SEE_MI, SEE_Seq, SEE_FT, SEE_FP, SEE_ASR};
    for i = 1:5
        myRes(nnn,i) = params.BW * SEEs{i}(end) / params.P_lim / 1e6;
    end
    fprintf('P_ave = %d W done.\n', Pave_vec(nnn));
end

colours = {'k','#0072BD','#D95319','#7E2F8E','#77AC30'};
markers = {'p','^','*','+','x'};
lgd_names = {'MSEE-MI','MSEE-Seq','SEE-FTrj','SEE-FPow','ASR-Seq'};
linestyle = {'-','--',':','-.'};

figure; hold on;
for i = 1:5
    p = plot(Pave_vec, myRes(:,i), 'LineWidth', 2);
    p.LineStyle = linestyle{mod(i,4)+1};
    p.Marker    = markers{i};
    p.Color     = colours{i};
end
xlabel('Average network TX power [W]', 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Min. secrecy energy efficiency [Mbits/Joule]', 'FontName', 'Times New Roman', 'FontSize', 11);
legend(lgd_names, 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

if ~exist('figures','dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig8_vs_Pave.fig');
saveas(gcf, 'figures/fig8_vs_Pave.png');
fprintf('Figure 8 saved.\n');
end

% =========================================================================

% RUN_FIGURE_9  Reproduce Figure 9: MSEE vs. UAV average flight power limit.
%
%   Sweeps P_lim in [180, 190, 200, 210, 220, 230, 250] W.
%   Corresponds to Fig. 9 in the paper (Section IV).

function run_figure_9()
addpath(genpath(fileparts(mfilename('fullpath'))));
if ~exist('results', 'dir'), mkdir('results'); end

Plim_vec = [180, 190, 200, 210, 220, 230, 250];
N        = 100;
myRes    = zeros(length(Plim_vec), 5);

for nnn = 1:length(Plim_vec)
    params       = system_params();
    params.N     = N;
    params.P_lim = Plim_vec(nnn);
    params.Consts = {params.P_lim, params.Qa, params.Qb, ...
                     params.af, params.beta0, params.N0, params.BW};
    fc = flight_constants();  params.nu0_ = fc.nu0;

    fname = sprintf('results/myResults_N%d_Plim%d.mat', N, Plim_vec(nnn));

    [~, SEE_MI]  = run_msee_mibcd(params);
    [~, SEE_Seq] = run_msee_seqbcd(params);
    [~, SEE_FT]  = run_msee_fixed_traj(params);
    [~, SEE_FP]  = run_msee_fixed_pow(params);
    [~, SEE_ASR] = run_asr_seq(params);

    SEE_MIBCD = SEE_MI; SEE_SeqBCD = SEE_Seq;
    SEE_FixedTrj = SEE_FT; SEE_FixedPow = SEE_FP; SEE_ASRoptim = SEE_ASR;
    save(fname, 'SEE_MIBCD', 'SEE_SeqBCD', 'SEE_FixedTrj', 'SEE_FixedPow', 'SEE_ASRoptim');

    SEEs = {SEE_MI, SEE_Seq, SEE_FT, SEE_FP, SEE_ASR};
    for i = 1:5
        myRes(nnn,i) = params.BW * SEEs{i}(end) / Plim_vec(nnn) / 1e6;
    end
    fprintf('P_lim = %d W done.\n', Plim_vec(nnn));
end

% Convert to kW for x-axis (N * P_lim / 1000)
x_vals = Plim_vec * N / 1e3;

colours = {'k','#0072BD','#D95319','#7E2F8E','#77AC30'};
markers = {'p','^','*','+','x'};
lgd_names = {'MSEE-MI','MSEE-Seq','SEE-FTrj','SEE-FPow','ASR-Seq'};
linestyle = {'-','--',':','-.'};

figure; hold on;
for i = 1:5
    p = plot(x_vals, myRes(:,i), 'LineWidth', 2);
    p.LineStyle = linestyle{mod(i,4)+1};
    p.Marker    = markers{i};
    p.Color     = colours{i};
end
xlabel('Flight power limit P_{lim} [kW total]', 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Min. secrecy energy efficiency [Mbits/Joule]', 'FontName', 'Times New Roman', 'FontSize', 11);
legend(lgd_names, 'Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

if ~exist('figures','dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig9_vs_Plim.fig');
saveas(gcf, 'figures/fig9_vs_Plim.png');
fprintf('Figure 9 saved.\n');
end
