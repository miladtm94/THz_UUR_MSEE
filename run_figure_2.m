% RUN_FIGURE_2  Reproduce Figure 2: MSEE convergence vs. iteration index.
%
%   Generates the convergence plot showing normalized MSEE [Mbits/Joule]
%   vs. outer BCD iteration index for all five schemes.
%   Corresponds to Fig. 2 in the paper (Section IV).
%
%   PREREQUISITE: Run main.m first (or at least the default N=100 case).

function run_figure_2()
addpath(genpath(fileparts(mfilename('fullpath'))));
params = system_params();

filename = sprintf('results/myResults_N%d_default.mat', params.N);
if ~exist(filename, 'file')
    error('run_figure_2: result file not found. Run main.m first.');
end
load(filename, 'SEE_MIBCD', 'SEE_SeqBCD', 'SEE_FixedTrj', 'SEE_FixedPow');

colours   = {'k',       '#0072BD', '#D95319', '#7E2F8E'};
markers   = {'none',    '^',       '*',       '+'};
linestyle = {'-',       '--',      ':',       '-.'};
lgd_names = {'MSEE-MI', 'MSEE-Seq', 'SEE-FTrj', 'SEE-FPow'};

SEE_arr = {SEE_MIBCD, SEE_SeqBCD, SEE_FixedTrj, SEE_FixedPow};
lens    = cellfun(@numel, SEE_arr);
l_max   = max(lens);

figure; hold on;
for i = 1:4
    SEE = SEE_arr{i};
    % Pad shorter curves with their final value for alignment
    SEE_padded = [SEE, repmat(SEE(end), 1, l_max - length(SEE))];
    p = plot(0:l_max-1, params.BW * SEE_padded / params.P_lim / 1e6, ...
        'LineWidth', 2);
    p.LineStyle     = linestyle{i};
    p.Marker        = markers{i};
    p.MarkerIndices = 1:3:l_max;
    p.Color         = colours{i};
end

xlabel('Iteration index', 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Min. secrecy energy efficiency [Mbits/Joule]', ...
    'FontName', 'Times New Roman', 'FontSize', 11);
lgd = legend(lgd_names);
lgd.Location = 'southeast';
lgd.FontSize = 10;
lgd.FontName = 'Times New Roman';
grid on;
xlim([0, l_max-1]);

if ~exist('figures', 'dir'), mkdir('figures'); end
saveas(gcf, 'figures/fig2_convergence.fig');
saveas(gcf, 'figures/fig2_convergence.png');
fprintf('Figure 2 saved to figures/fig2_convergence.*\n');
end
