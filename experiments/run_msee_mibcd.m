function [Itr, SEE] = run_msee_mibcd(params)
% RUN_MSEE_MIBCD  Algorithm 3: Maximum-Improvement BCD for MSEE maximization.
%
%   [Itr, SEE] = run_msee_mibcd(params)
%
%   Implements Algorithm 3 (MSEE-MI) from Section III-E of the paper:
%   a greedy Block Coordinate Descent (BCD) algorithm that at each
%   iteration solves all four subproblems (P1)-(P4) in parallel, then
%   commits only the update that yields the greatest MSEE improvement.
%
%   ALGORITHM 3 SUMMARY (Section III-E):
%   1. At iteration l, solve (P1), (P2), (P3), (P4) with current values
%      of all OTHER blocks held fixed (parallel evaluation).
%   2. Evaluate the MSEE improvement that each subproblem's solution gives.
%   3. Accept the update with the largest improvement; keep others unchanged.
%   4. Repeat until fractional MSEE increase < eps_bcd.
%
%   CONVERGENCE:
%   The monotone non-decrease of MSEE follows from the max-improvement
%   selection rule; the sequence is bounded above by a finite MSEE value
%   (proof sketch in Section III-E of the paper). This algorithm converges
%   faster than Algorithm 2 (MSEE-Seq) but requires more computation per
%   iteration as all four subproblems are solved every step.
%
%   INPUTS:
%       params  : System parameter struct from system_params()
%
%   OUTPUTS:
%       Itr     : Struct array of length (total iterations), each entry:
%                   .Trj     - UAV trajectory [NxK]
%                   .Vel     - UAV velocity   [Nx2]
%                   .usrPow  - UE power       [NxK]
%                   .uavPow  - UAV relay power [Nx1]
%                   .bsPow   - BS jammer power [Nx1]
%                   .usrSch  - Scheduling      [NxK]
%                   .prop    - Propulsion power [Nx1]
%       SEE     : MSEE value at each iteration [1 x total_iters]
%
% See also: run_msee_seqbcd, run_msee_fixed_traj, run_msee_fixed_pow.

%% ---- Setup ----
N    = params.N;
K    = params.K;

% Get feasible initial point
[Q_uav, V_uav, xu, yu, vx, vy, P_k, P_u, P_b, zeta, ~] = feasible_init(params);

l = 1;
[SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
    see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);

Itr(l) = pack_iteration(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

% Dinkelbach tracking arrays (pre-allocated for P4 inner loop)
MAX_DINKL = 20;
lambda_frac     = zeros(MAX_DINKL, MAX_DINKL);
ASR_approx      = zeros(MAX_DINKL, MAX_DINKL);
Pflight_approx  = zeros(MAX_DINKL, MAX_DINKL);
psi_gap         = ones(MAX_DINKL, MAX_DINKL);
SEE_approx      = zeros(MAX_DINKL, MAX_DINKL);

fprintf('MSEE-MI (Algorithm 3) optimization started.\n');

%% ---- Main BCD Loop ----
while true
    l = l + 1;
    fprintf('\n--- Outer iteration %d ---\n', l-1);

    SEE_candidate = zeros(1, 4);

    %% -- (P1): Joint UE scheduling + transmit power --
    lambda_n = params.P_lim / sum(P_flight) / (2*log(2));
    B_kn = g_ku ./ (P_b .* g_ub + 1);
    C_n  = P_u  .* g_ub;
    D_kn = (g_ub .* (P_u + P_b) + 1) ./ g_ku;

    tStart = tic;
    [P_k_, zeta_, ~] = usr_sched_pow_optim( ...
        P_k, zeta, params.Pa_peak, params.Pa_ave, lambda_n, B_kn, C_n, D_kn);
    fprintf('  (P1) solved in %.2f s\n', toc(tStart));
    [SEE_candidate(1), ~, ~, ~, ~, ~] = ...
        see_compute(P_k_, zeta_, P_b, P_u, Q_uav, V_uav, params);

    %% -- (P2): UAV relay transmit power --
    lambda_kn = params.P_lim * zeta ./ sum(P_flight) / (2*log(2));
    E_kn = P_k .* g_ku;
    F_kn = (P_k .* g_ku + P_b .* g_ub + 1) ./ g_ub;
    G_k  = sum(params.P_lim * zeta .* log2(1 + snr_uav), 2) ./ sum(P_flight) / 2;

    tStart = tic;
    [P_u_, ~] = uav_relay_pow_optim(P_u, params.Pu_peak, params.Pu_ave, ...
        lambda_kn, E_kn, F_kn, G_k, zeta);
    fprintf('  (P2) solved in %.2f s\n', toc(tStart));
    [SEE_candidate(2), ~, ~, ~, ~, ~] = ...
        see_compute(P_k, zeta, P_b, P_u_, Q_uav, V_uav, params);

    %% -- (P3): BS cooperative jamming power --
    lambda_kn = params.P_lim * zeta ./ sum(P_flight) / (2*log(2));
    H_kn = g_ku .* P_k .* P_u;
    I_kn = (P_u .* g_ub + P_k .* g_ku + 1) ./ g_ub;
    J_kn = P_k .* g_ku ./ g_ub;
    K_n  = 1 ./ g_ub;

    tStart = tic;
    [P_b_, ~] = bs_jamming_pow_optim(P_b, params.Pb_peak, params.Pb_ave, ...
        lambda_kn, H_kn, I_kn, J_kn, K_n, zeta);
    fprintf('  (P3) solved in %.2f s\n', toc(tStart));
    [SEE_candidate(3), ~, ~, ~, ~, ~] = ...
        see_compute(P_k, zeta, P_b_, P_u, Q_uav, V_uav, params);

    %% -- (P4): Joint trajectory and velocity (Dinkelbach, Algorithm 1) --
    epsilon = 1e-7;
    for k = 1:K
        un(k,:) = norms([Q_uav - params.Qa(:,k)]) - epsilon;
    end
    k0 = (P_u + P_b) ./ (P_k .* P_u);
    k1 = 1 ./ P_u;
    k2 = P_k;
    k3 = P_b;
    sn = 1 ./ g_ku + epsilon;
    wn = 1 ./ g_ub + epsilon;

    mu = sqrt(sqrt(1 + norms(V_uav).^4/(4*params.nu0_^4)) ...
              - norms(V_uav).^2 / (2*params.nu0_^2)) + 1e-2;

    [SEE_orig_0, ~, ~, ~, ~, Pf0] = see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    m = 1;
    lambda_frac(l-1, m) = SEE_orig_0;
    Pflight_approx(l-1, m) = mean(Pf0 / params.P_lim);
    ASR_approx(l-1, m)    = lambda_frac(l-1,m) * Pflight_approx(l-1,m);

    Q_state = [xu', yu', vx', vy', mu', un', sn', wn', zeros(N, 3)];

    tStart = tic;
    while true
        Q_new = trajectory_optim(Q_state, params.H, params.Vmax, params.P_lim, ...
            params.N0, params.beta0, params.af, params.Qa, params.Qb, ...
            params.QI(1:2), params.dt, 0, params.R2, k0, k1, k2, k3, zeta, ...
            lambda_frac(l-1, m));
        m = m + 1;
        xu_new = Q_new(:,1)';  yu_new = Q_new(:,2)';
        vx_new = Q_new(:,3)';  vy_new = Q_new(:,4)';
        Quav_new = [xu_new; yu_new; repmat(params.H,1,N)];
        Vuav_new = [vx_new; vy_new];

        mu = sqrt(sqrt(1 + norms(Vuav_new).^4/(4*params.nu0_^4)) ...
                  - norms(Vuav_new).^2 / (2*params.nu0_^2)) + 1e-2;
        for k = 1:K
            un(k,:) = norms([Quav_new - params.Qa(:,k)]) - epsilon;
        end
        [SEE_new, g_ku_t, g_ub_t, ~, ~, Pf_new] = ...
            see_compute(P_k, zeta, P_b, P_u, Quav_new, Vuav_new, params);
        sn = 1 ./ g_ku_t + epsilon;
        wn = 1 ./ g_ub_t + epsilon;

        Pflight_approx(l-1,m)  = mean(Pf_new / params.P_lim);
        ASR_approx(l-1,m)      = Pflight_approx(l-1,m) * SEE_new;
        lambda_frac(l-1,m)     = SEE_new;
        psi_gap(l-1, m-1) = ASR_approx(l-1,m-1) - lambda_frac(l-1,m) * Pflight_approx(l-1,m-1);

        if abs(psi_gap(l-1,m-1)) <= params.eps_dinkelbch || psi_gap(l-1,m-1) >= 0
            fprintf('  (P4) Dinkelbach converged at inner iter %d (gap=%.2e)\n', ...
                m-1, psi_gap(l-1,m-1));
            Q_uav_ = Quav_new;  V_uav_ = Vuav_new;
            xu_ = xu_new;       yu_ = yu_new;
            vx_ = vx_new;       vy_ = vy_new;
            break;
        else
            Q_state = [xu_new', yu_new', vx_new', vy_new', mu', un', sn', wn', zeros(N,3)];
        end
    end
    fprintf('  (P4) solved in %.2f s\n', toc(tStart));
    [SEE_candidate(4), ~, ~, ~, ~, ~] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav_, V_uav_, params);

    %% -- Greedy selection: accept the block with maximum improvement --
    [max_improvement, best_idx] = max(SEE_candidate);

    if max_improvement <= 0
        warning('run_msee_mibcd: No improving update found. Terminating.');
        break;
    end

    switch best_idx
        case 1
            P_k   = P_k_;    zeta  = zeta_;
            fprintf('  -> Selected: (P1) joint UE scheduling + power\n');
        case 2
            P_u   = P_u_;
            fprintf('  -> Selected: (P2) UAV relay power\n');
        case 3
            P_b   = P_b_;
            fprintf('  -> Selected: (P3) BS jamming power\n');
        case 4
            Q_uav = Q_uav_;  V_uav = V_uav_;
            xu = xu_;  yu = yu_;  vx = vx_;  vy = vy_;
            fprintf('  -> Selected: (P4) trajectory & velocity\n');
    end

    %% -- Re-evaluate MSEE with the accepted update --
    [SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_iteration(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    F_algo = (SEE(l) - SEE(l-1)) / SEE(l-1);
    fprintf('  Fractional MSEE increase: %.6f\n', F_algo);

    if F_algo <= params.eps_bcd
        fprintf('MSEE-MI converged after %d outer iterations.\n', l-1);
        break;
    end
end

end

%% ---- Helper: pack one iteration's state ----
function s = pack_iteration(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight)
s.usrPow = P_k';
s.bsPow  = P_b';
s.uavPow = P_u';
s.usrSch = zeta';
s.Trj    = Q_uav';
s.Vel    = V_uav';
s.prop   = P_flight';
end
