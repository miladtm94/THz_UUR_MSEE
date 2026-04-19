function [fig_traj, ax_traj, lgd_traj] = plot_network_topology(params)
% PLOT_NETWORK_TOPOLOGY  Visualize the 2D network layout and UAV region.
%
%   [fig_traj, ax_traj, lgd_traj] = plot_network_topology(params)
%
%   Draws the circular annulus region [R1, R2], base station (BS),
%   ground UE positions, and the UAV's initial/final depot location.
%   Returns handles for subsequent trajectory overlays via plot_uav_trajectory.
%
%   INPUTS:
%       params : System parameter struct from system_params()
%
%   OUTPUTS:
%       fig_traj  : Figure handle
%       ax_traj   : Axes handle
%       lgd_traj  : Legend handle
%
% See also: plot_uav_trajectory, system_params.

N    = params.N;
R1   = params.R1;
R2   = params.R2;
Ru   = params.Ru;
Qa   = params.Qa;
Qb   = params.Qb;
K    = params.K;

theta = linspace(0, 2*pi, 200);

fig_traj      = figure;
fig_traj.Name = 'UAV Trajectory';
ax_traj       = axes(fig_traj);

% Inner and outer boundary circles
plot(ax_traj, R1*cos(theta), R1*sin(theta), '--r', 'LineWidth', 1.2);
hold(ax_traj, 'on');
plot(ax_traj, R2*cos(theta), R2*sin(theta), '--b', 'LineWidth', 1.2);

% Ground UE positions
xa = Qa(1,:);  ya = Qa(2,:);
h_ue = scatter(ax_traj, xa, ya, 60, 'm', 'filled', '^');
for k = 1:K
    text(xa(k)-1, ya(k)-3, sprintf('UE%d', k), 'FontSize', 8);
end

% Base station at origin
h_bs = scatter(ax_traj, Qb(1), Qb(2), 100, 'g', 'filled', 's');
text(Qb(1)-1, Qb(2)-3, 'BS', 'FontSize', 8);

% UAV initial/final depot
h_depot = scatter(ax_traj, Ru, 0, 120, 'r', 'filled', 'h');

axis(ax_traj, [-(1.2*R2), (1.2*R2), -(1.2*R2), (1.2*R2)]);
xlabel(ax_traj, 'x [m]');
ylabel(ax_traj, 'y [m]');
lgd_traj = legend(ax_traj, h_depot, {'Init./Final Loc.'}, ...
    'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
grid(ax_traj, 'on');

end

% =========================================================================

function plot_uav_trajectory(xu, yu, ax_traj, lgd_traj, scheme_idx, iter_label)
% PLOT_UAV_TRAJECTORY  Overlay one UAV trajectory on the network plot.
%
%   plot_uav_trajectory(xu, yu, ax_traj, lgd_traj, scheme_idx, iter_label)
%
%   Adds a directed trajectory curve (with arrows showing flight direction)
%   to an existing axes from plot_network_topology.
%
%   INPUTS:
%       xu, yu      : UAV horizontal positions [1xN]
%       ax_traj     : Axes handle from plot_network_topology
%       lgd_traj    : Legend handle to update
%       scheme_idx  : Integer 1-5 indexing the scheme color/marker
%       iter_label  : String label shown in the legend
%
% See also: plot_network_topology, plot_trajectory_dir.

colours   = {'k', '#0072BD', '#D95319', '#7E2F8E', '#77AC30'};
markers   = {'o', 'p',       'd',       's',       'h'};
lgd_names = {'SEE-MI', 'SEE-Seq', 'SEE-FTrj', 'SEE-FPow', 'ASR-Seq'};

hold(ax_traj, 'on');
[h_line, ~] = plot_trajectory_dir(ax_traj, xu, yu);

h_line.LineStyle     = '-';
h_line.Marker        = markers{scheme_idx};
h_line.MarkerSize    = 6;
h_line.MarkerIndices = 1 : round(length(xu)/10) : length(xu);
h_line.Color         = colours{scheme_idx};
h_line.LineWidth     = 1.5;

if scheme_idx <= length(lgd_names)
    lgd_traj.String{end+1} = lgd_names{scheme_idx};
end

end

% =========================================================================

function [h_line, h_arrows] = plot_trajectory_dir(ax, vX, vY)
% PLOT_TRAJECTORY_DIR  Plot an (x,y) curve with direction-indicating arrows.
%
%   [h_line, h_arrows] = plot_trajectory_dir(ax, vX, vY)
%
%   Draws the trajectory as a line and adds quiver arrows between
%   consecutive points to indicate the direction of motion.
%
%   INPUTS:
%       ax    : Target axes handle
%       vX, vY: Position vectors [1xN]
%
%   OUTPUTS:
%       h_line   : Line plot handle
%       h_arrows : Quiver handle (excluded from legend)
%
% See also: plot_uav_trajectory.

n   = length(vX);
idx = 1:(n-1);

dX = vX(idx+1) - vX(idx);
dY = vY(idx+1) - vY(idx);

h_line   = plot(ax, vX, vY);  hold(ax, 'on');
h_arrows = quiver(ax, vX(idx), vY(idx), dX, dY, 0);
h_arrows.Annotation.LegendInformation.IconDisplayStyle = 'off';
grid(ax, 'on');
axis(ax, 'equal');

end
