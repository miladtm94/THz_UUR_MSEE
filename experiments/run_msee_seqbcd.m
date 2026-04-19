function [Itr, SEE] = run_msee_seqbcd(params)
% RUN_MSEE_SEQBCD  Algorithm 2: Sequential BCD for MSEE maximization.
%
%   [Itr, SEE] = run_msee_seqbcd(params)
%
%   Implements Algorithm 2 (MSEE-Seq) from Section III-E of the paper:
%   a sequential Block Coordinate Descent algorithm that solves subproblems
%   (P1) -> (P2) -> (P3) -> (P4) in order, immediately using each
%   subproblem's solution to warm-start the next.
%
%   ALGORITHM 2 SUMMARY (Section III-E):
%   1. Solve (P1): update (P_k, zeta) given (P_u, P_b, Q, V).
%   2. Solve (P2): update P_u given updated (P_k, zeta) and (P_b, Q, V).
%   3. Solve (P3): update P_b given updated (P_k, zeta, P_u) and (Q, V).
%   4. Solve (P4): run Dinkelbach inner loop to update (Q, V) given
%      the latest (P_k, zeta, P_u, P_b).
%   5. Compute MSEE; repeat until fractional increase < eps_bcd.
%
%   CONVERGENCE (Eq. 52, Section III-E):
%   Each sequential update is non-decreasing in MSEE (proved via
%   equality of objective at the expansion point and tightness of the
%   SCA bounds), so the overall objective is monotonically non-decreasing.
%   Combined with the compactness of the feasible set, Algorithm 2
%   converges to at least a stationary point.
%
%   INPUTS:
%       params : System parameter struct from system_params()
%
%   OUTPUTS:
%       Itr    : Struct array (see run_msee_mibcd for field definitions)
%       SEE    : MSEE trajectory over sub-iterations [1 x total_sub_iters]
%
% See also: run_msee_mibcd, run_msee_fixed_traj, run_msee_fixed_pow.

%% ---- Setup ----
N  = params.N;
K  = params.K;
fc = flight_constants();

[Q_uav, V_uav, xu, yu, vx, vy, P_k, P_u, P_b, zeta, ~] = feasible_init(params);

l  = 1;
ll = 1;  % Outer iteration counter (one outer iter = all 4 subproblems)
[SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
    see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

% Dinkelbach bookkeeping
MAX_D = 20;
lambda_frac    = zeros(MAX_D, MAX_D);
ASR_approx     = zeros(MAX_D, MAX_D);
Pflight_approx = zeros(MAX_D, MAX_D);
psi_gap        = ones(MAX_D, MAX_D);

fprintf('MSEE-Seq (Algorithm 2) optimization started.\n');

%% ---- Main BCD Loop ----
while true
    fprintf('\n--- Outer iteration %d ---\n', ll);

    %% -- (P1): Joint UE scheduling + transmit power --
    l = l + 1;
    lambda_n = params.P_lim / sum(P_flight) / (2*log(2));
    B_kn = g_ku ./ (P_b .* g_ub + 1);
    C_n  = P_u  .* g_ub;
    D_kn = (g_ub .* (P_u + P_b) + 1) ./ g_ku;

    tStart = tic;
    [P_k, zeta, ~] = usr_sched_pow_optim( ...
        P_k, zeta, params.Pa_peak, params.Pa_ave, lambda_n, B_kn, C_n, D_kn);
    fprintf('  (P1) converged in %.2f s\n', toc(tStart));

    [SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    %% -- (P2): UAV relay transmit power --
    l = l + 1;
    lambda_kn = params.P_lim * zeta ./ sum(P_flight) / (2*log(2));
    E_kn = P_k .* g_ku;
    F_kn = (P_k .* g_ku + P_b .* g_ub + 1) ./ g_ub;
    G_k  = sum(params.P_lim * zeta .* log2(1 + snr_uav), 2) ./ sum(P_flight) / 2;

    tStart = tic;
    [P_u, ~] = uav_relay_pow_optim(P_u, params.Pu_peak, params.Pu_ave, ...
        lambda_kn, E_kn, F_kn, G_k, zeta);
    fprintf('  (P2) converged in %.2f s\n', toc(tStart));

    [SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    %% -- (P3): BS cooperative jamming power --
    l = l + 1;
    lambda_kn = params.P_lim * zeta ./ sum(P_flight) / (2*log(2));
    H_kn = g_ku .* P_k .* P_u;
    I_kn = (P_u .* g_ub + P_k .* g_ku + 1) ./ g_ub;
    J_kn = P_k .* g_ku ./ g_ub;
    K_n  = 1 ./ g_ub;

    tStart = tic;
    [P_b, ~] = bs_jamming_pow_optim(P_b, params.Pb_peak, params.Pb_ave, ...
        lambda_kn, H_kn, I_kn, J_kn, K_n, zeta);
    fprintf('  (P3) converged in %.2f s\n', toc(tStart));

    [SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    %% -- (P4): Joint trajectory and velocity (Dinkelbach inner loop) --
    l = l + 1;

    epsilon = 1e-7;
    un = zeros(K, N);
    for k = 1:K
        un(k,:) = norms([Q_uav - params.Qa(:,k)]) - epsilon;
    end
    k0 = (P_u + P_b) ./ (P_k .* P_u);
    k1 = 1 ./ P_u;
    k2 = P_k;
    k3 = P_b;
    sn = 1 ./ g_ku + epsilon;
    wn = 1 ./ g_ub + epsilon;
    mu = sqrt(sqrt(1 + norms(V_uav).^4/(4*fc.nu0^4)) ...
              - norms(V_uav).^2/(2*fc.nu0^2)) + 1e-2;

    [SEE0, ~, ~, ~, ~, Pf0] = see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    m = 1;
    lambda_frac(ll, m)    = SEE0;
    Pflight_approx(ll, m) = mean(Pf0 / params.P_lim);
    ASR_approx(ll, m)     = lambda_frac(ll,m) * Pflight_approx(ll,m);
    Q_state = [xu', yu', vx', vy', mu', un', sn', wn', zeros(N, 3)];

    tStart = tic;
    while true
        Q_new = trajectory_optim(Q_state, params.H, params.Vmax, params.P_lim, ...
            params.N0, params.beta0, params.af, params.Qa, params.Qb, ...
            params.QI(1:2), params.dt, 0, params.R2, k0, k1, k2, k3, zeta, ...
            lambda_frac(ll, m));
        m = m + 1;
        xu = Q_new(:,1)';  yu = Q_new(:,2)';
        vx = Q_new(:,3)';  vy = Q_new(:,4)';
        mu = Q_new(:,5)';
        un = Q_new(:, 6:5+K)';
        sn = Q_new(:, 6+K:5+2*K)';
        wn = Q_new(:, 6+2*K)';
        Quav_new = [xu; yu; repmat(params.H,1,N)];
        Vuav_new = [vx; vy];

        [SEE_m, g_ku, g_ub, ~, ~, Pf_m] = ...
            see_compute(P_k, zeta, P_b, P_u, Quav_new, Vuav_new, params);

        Pflight_approx(ll,m) = mean(Pf_m / params.P_lim);
        ASR_approx(ll,m)     = Pflight_approx(ll,m) * SEE_m;
        lambda_frac(ll,m)    = SEE_m;
        psi_gap(ll,m-1) = ASR_approx(ll,m-1) - lambda_frac(ll,m) * Pflight_approx(ll,m-1);

        SEE(l) = SEE_m;   % Track within Dinkelbach

        if abs(psi_gap(ll,m-1)) <= params.eps_dinkelbch || psi_gap(ll,m-1) >= 0
            fprintf('  (P4) Dinkelbach converged at inner iter %d (gap=%.2e)\n', ...
                m-1, psi_gap(ll,m-1));
            Q_uav = Quav_new;  V_uav = Vuav_new;
            break;
        else
            Q_state = [xu', yu', vx', vy', mu', un', sn', wn', zeros(N,3)];
        end
    end
    fprintf('  (P4) solved in %.2f s\n', toc(tStart));

    [SEE(l), g_ku, g_ub, snr_bs, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    %% -- Check outer convergence --
    ll = ll + 1;
    F_outer = (SEE(ll) - SEE(ll-1)) / SEE(ll-1);
    fprintf('  Outer fractional MSEE increase: %.6f\n', F_outer);
    if F_outer <= params.eps_bcd
        fprintf('MSEE-Seq converged after %d outer iterations.\n', ll-1);
        break;
    end
end

end

function s = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight)
s.usrPow = P_k';  s.bsPow = P_b';  s.uavPow = P_u';
s.usrSch = zeta'; s.Trj = Q_uav';  s.Vel = V_uav';  s.prop = P_flight';
end
