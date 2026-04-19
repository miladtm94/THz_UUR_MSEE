function [Q_uav, V_uav, xu, yu, vx, vy, P_k, P_u, P_b, zeta, T] = ...
    feasible_init(params)
% FEASIBLE_INIT  Compute a feasible initial point for the MSEE optimization.
%
%   [Q_uav, V_uav, xu, yu, vx, vy, P_k, P_u, P_b, zeta, T] = ...
%       feasible_init(params)
%
%   Constructs a feasible initial trajectory and power allocation that
%   satisfies all constraints C1-C14 from Eq. (3)-(13). The trajectory is
%   initialized as either a circular or piriform (cycloidal) path, depending
%   on whether N is sufficient for a full circle at Vmax.
%
%   TRAJECTORY INITIALIZATION LOGIC:
%   1. Determine the minimum N for a circular path (Nmin_circular) such that
%      the speed and acceleration limits are not violated.
%   2. If N >= Nmin_circular: use a circular trajectory of radius Ru.
%   3. If Nmin_cyclic <= N < Nmin_circular: use a piriform (lemniscate-like)
%      cyclic path that fits within the region while satisfying kinematics.
%   4. Otherwise: display a warning and request a larger N.
%
%   POWER INITIALIZATION:
%   - P_b(n) = Pb_ave - eps   (just below average limit)
%   - P_u(n) = Pu_ave - eps   (just below average limit)
%   - zeta(k,n): time-division scheduling, each UE assigned c = ceil(N/K) slots
%   - P_k(k,n) = Pa_ave * zeta(k,n)
%
%   INPUTS:
%       params  : System parameter struct from system_params()
%
%   OUTPUTS:
%       Q_uav   : Initial UAV 3D trajectory [3xN]
%       V_uav   : Initial UAV velocity [2xN]
%       xu, yu  : Horizontal UAV positions [1xN]
%       vx, vy  : UAV velocity components [1xN]
%       P_k     : Initial UE power schedule [KxN]
%       P_u     : Initial UAV relay power [1xN]
%       P_b     : Initial BS jamming power [1xN]
%       zeta    : Initial user scheduling [KxN]
%       T       : Total flight mission duration [s]
%
%   NOTE:
%       The initialization does NOT guarantee optimality, only feasibility.
%       The BCD algorithm (Algorithms 2 and 3) refines from this point.
%
% See also: system_params, flight_constants, see_compute.

N       = params.N;
K       = params.K;
Vmax    = params.Vmax;
Amax    = params.Amax;
Ru      = params.Ru;
H       = params.H;
Qb      = params.Qb;
dt      = params.dt;
Pa_ave  = params.Pa_ave;
Pb_ave  = params.Pb_ave;
Pu_ave  = params.Pu_ave;
QI      = params.QI;

epsilon = 1e-3;   % Small offset to keep power strictly below average limit

%% -----------------------------------------------------------------
%  Step 1: Determine minimum N for a feasible circular trajectory
% ------------------------------------------------------------------
j = 3;
while true
    theta_circ = linspace(0, 2*pi, j);
    xu_c = Ru * cos(theta_circ);
    yu_c = Ru * sin(theta_circ);
    yu_c(end) = 0;
    vx_c = [diff(xu_c), 0] / dt;
    vy_c = [diff(yu_c), 0] / dt;
    if max(norms([vx_c; vy_c])) <= Vmax && ...
       max(norms([diff(vx_c); diff(vy_c)])) <= Amax
        break;
    end
    j = j + 1;
end
Nmin_circular = j;

%% -----------------------------------------------------------------
%  Step 2: Determine minimum N for a feasible cyclic (piriform) path
% ------------------------------------------------------------------
NN = 3;
while NN < Nmin_circular
    t_cyc = linspace(pi/2, 5*pi/2, NN);
    x_cyc = Ru * (sin(t_cyc) + 1) / 2;
    y_cyc = Ru * (1 - sin(t_cyc)) .* cos(t_cyc);
    vx_c2 = [diff(x_cyc), 0] / dt;
    vy_c2 = [diff(y_cyc), 0] / dt;
    if max(norms([vx_c2; vy_c2])) <= Vmax && ...
       max(norms([diff(vx_c2); diff(vy_c2)])) <= Amax
        break;
    end
    NN = NN + 1;
end
Nmin_cyclic = NN;

fprintf('Min N for circular trajectory: %d\n', Nmin_circular);
fprintf('Min N for cyclic (piriform) trajectory: %d\n', Nmin_cyclic);

%% -----------------------------------------------------------------
%  Step 3: Select and construct the trajectory
% ------------------------------------------------------------------
if N >= Nmin_circular
    % Circular trajectory of radius Ru (constraint C14 satisfied by construction)
    theta  = linspace(0, 2*pi, N);
    xu     = Ru * cos(theta);
    yu     = Ru * sin(theta);
    yu(end) = 0;    % Enforce return to initial x-axis point (C10)
    vx     = [diff(xu), 0] / dt;
    vy     = [diff(yu), 0] / dt;
    fprintf('Circular trajectory initialized with radius %.1f m\n', Ru);

elseif N >= Nmin_cyclic
    % Piriform-based cyclic trajectory (scaled to satisfy kinematics)
    a = Ru;
    while true
        t_cyc = linspace(pi/2, 5*pi/2, N);
        xu    = Ru * (sin(t_cyc) + 1) / 2;
        yu    = a  * (1 - sin(t_cyc)) .* cos(t_cyc);
        vx    = [diff(xu), 0] / dt;
        vy    = [diff(yu), 0] / dt;
        if max(norms([vx; vy])) <= Vmax && ...
           max(norms([diff(vx); diff(vy)])) <= Amax
            break;
        end
        a = a - 1;
    end
    fprintf('Cyclic piriform trajectory initialized (amplitude a = %.1f m)\n', a);

else
    error('feasible_init:insufficientN', ...
        'N = %d is too small for a cyclic flight. Increase N to at least %d.', ...
        N, Nmin_cyclic);
end

fprintf('Max UAV speed at initialization: %.3f m/s\n', ...
        max(norms([diff(xu); diff(yu)] / dt)));

%% -----------------------------------------------------------------
%  Step 4: Assemble trajectory matrices
% ------------------------------------------------------------------
Q_uav = [xu; yu; repmat(H, 1, N)];   % [3xN]
V_uav = [vx; vy];                     % [2xN]

%% -----------------------------------------------------------------
%  Step 5: Initialize transmit powers
% ------------------------------------------------------------------
P_b = repmat(Pb_ave - epsilon, 1, N);   % [1xN]
P_u = repmat(Pu_ave - epsilon, 1, N);   % [1xN]

%% -----------------------------------------------------------------
%  Step 6: Initialize user scheduling (round-robin TDMA)
% ------------------------------------------------------------------
c    = ceil(N / K);
base = [ones(1, c), zeros(1, N - c)];
zeta = zeros(K, N);
for k = 1:K
    zeta(k, :) = base;
    base = circshift(base, [1, c]);
end

% Relax binary to continuous: push toward 0/1 but avoid strict boundary
zeta(zeta == 0) = epsilon;
zeta(zeta == 1) = 1 - K * epsilon;

% Initial UE transmit power proportional to scheduling variable
P_k = repmat(Pa_ave - epsilon, K, N) .* zeta;

%% -----------------------------------------------------------------
%  Step 7: Compute mission time
% ------------------------------------------------------------------
T = N * dt;   % Total flight mission duration [s]

end
