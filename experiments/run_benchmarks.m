function [Itr, SEE] = run_msee_fixed_traj(params)
% RUN_MSEE_FIXED_TRAJ  Benchmark: MSEE maximization with fixed trajectory.
%
%   [Itr, SEE] = run_msee_fixed_traj(params)
%
%   Implements the "SEE-FTrj" benchmark from Section IV of the paper:
%   optimizes UE power (P_k, zeta), UAV relay power (P_u), and BS jamming
%   power (P_b) via the sequential BCD approach, but keeps the UAV
%   trajectory FIXED at the initial circular/piriform path throughout.
%
%   This benchmark isolates the contribution of trajectory optimization
%   to MSEE performance. Its MSEE is lower than MSEE-Seq and MSEE-MI
%   because the trajectory is not optimized.
%
%   INPUTS:
%       params : System parameter struct from system_params()
%
%   OUTPUTS:
%       Itr    : Struct array (see run_msee_mibcd for field definitions)
%       SEE    : MSEE at each sub-iteration
%
% See also: run_msee_mibcd, run_msee_seqbcd, run_msee_fixed_pow.

N = params.N;
K = params.K;

[Q_uav, V_uav, ~, ~, ~, ~, P_k, P_u, P_b, zeta, ~] = feasible_init(params);

l = 1;
[SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
    see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

fprintf('SEE-FTrj (fixed-trajectory benchmark) started.\n');

while true
    l = l + 1;

    % (P1)
    lambda_n = params.P_lim / sum(P_flight) / (2*log(2));
    B_kn = g_ku ./ (P_b .* g_ub + 1);
    C_n  = P_u  .* g_ub;
    D_kn = (g_ub .* (P_u + P_b) + 1) ./ g_ku;
    [P_k, zeta, ~] = usr_sched_pow_optim( ...
        P_k, zeta, params.Pa_peak, params.Pa_ave, lambda_n, B_kn, C_n, D_kn);
    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);
    l = l + 1;

    % (P2)
    lambda_kn = params.P_lim * zeta ./ sum(P_flight) / (2*log(2));
    E_kn = P_k .* g_ku;
    F_kn = (P_k .* g_ku + P_b .* g_ub + 1) ./ g_ub;
    G_k  = sum(params.P_lim * zeta .* log2(1 + snr_uav), 2) ./ sum(P_flight) / 2;
    [P_u, ~] = uav_relay_pow_optim(P_u, params.Pu_peak, params.Pu_ave, ...
        lambda_kn, E_kn, F_kn, G_k, zeta);
    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);
    l = l + 1;

    % (P3)
    lambda_kn = params.P_lim * zeta ./ sum(P_flight) / (2*log(2));
    H_kn = g_ku .* P_k .* P_u;
    I_kn = (P_u .* g_ub + P_k .* g_ku + 1) ./ g_ub;
    J_kn = P_k .* g_ku ./ g_ub;
    K_n  = 1 ./ g_ub;
    [P_b, ~] = bs_jamming_pow_optim(P_b, params.Pb_peak, params.Pb_ave, ...
        lambda_kn, H_kn, I_kn, J_kn, K_n, zeta);
    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    F_algo = (SEE(l) - SEE(l-3)) / SEE(l-3);
    fprintf('  Outer fractional MSEE increase: %.6f\n', F_algo);
    if F_algo <= params.eps_bcd
        fprintf('SEE-FTrj converged.\n');
        break;
    end
end

end

% =========================================================================

function [Itr, SEE] = run_msee_fixed_pow(params)
% RUN_MSEE_FIXED_POW  Benchmark: MSEE with fixed powers, trajectory only.
%
%   [Itr, SEE] = run_msee_fixed_pow(params)
%
%   Implements the "SEE-FPow" benchmark from Section IV: optimizes only
%   the UAV trajectory (Q, V) via Dinkelbach + SCA, while keeping all
%   transmit powers fixed at their initial values.
%
%   This isolates the contribution of power allocation to MSEE.
%
%   INPUTS:
%       params : System parameter struct from system_params()
%
%   OUTPUTS:
%       Itr, SEE : see run_msee_mibcd
%
% See also: run_msee_fixed_traj.

N  = params.N;
K  = params.K;
fc = flight_constants();

[Q_uav, V_uav, xu, yu, vx, vy, P_k, P_u, P_b, zeta, ~] = feasible_init(params);

l = 1;
[SEE(l), g_ku, g_ub, ~, ~, P_flight] = ...
    see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

MAX_D = 20;
lambda_frac    = zeros(MAX_D, MAX_D);
Pflight_approx = zeros(MAX_D, MAX_D);
ASR_approx     = zeros(MAX_D, MAX_D);
psi_gap        = ones(MAX_D, MAX_D);

fprintf('SEE-FPow (fixed-power benchmark) started.\n');

outer = 1;
while true
    outer = outer + 1;
    l = l + 1;

    epsilon = 1e-7;
    un = zeros(K, N);
    for k = 1:K
        un(k,:) = norms([Q_uav - params.Qa(:,k)]) - epsilon;
    end
    k0 = (P_u + P_b) ./ (P_k .* P_u);
    k1 = 1 ./ P_u;  k2 = P_k;  k3 = P_b;
    sn = 1 ./ g_ku + epsilon;
    wn = 1 ./ g_ub + epsilon;
    mu = sqrt(sqrt(1 + norms(V_uav).^4/(4*fc.nu0^4)) ...
              - norms(V_uav).^2/(2*fc.nu0^2)) + 1e-2;

    [SEE0, ~, ~, ~, ~, Pf0] = see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    m = 1;
    lambda_frac(outer,m)    = SEE0;
    Pflight_approx(outer,m) = mean(Pf0 / params.P_lim);
    ASR_approx(outer,m)     = lambda_frac(outer,m) * Pflight_approx(outer,m);
    Q_state = [xu', yu', vx', vy', mu', un', sn', wn', zeros(N,3)];

    tStart = tic;
    while true
        Q_new = trajectory_optim(Q_state, params.H, params.Vmax, params.P_lim, ...
            params.N0, params.beta0, params.af, params.Qa, params.Qb, ...
            params.QI(1:2), params.dt, 0, params.R2, k0, k1, k2, k3, zeta, ...
            lambda_frac(outer,m));
        m = m + 1;
        xu = Q_new(:,1)';  yu = Q_new(:,2)';
        vx = Q_new(:,3)';  vy = Q_new(:,4)';
        mu = Q_new(:,5)';
        un = Q_new(:,6:5+K)';
        sn = Q_new(:,6+K:5+2*K)';
        wn = Q_new(:,6+2*K)';
        Quav_new = [xu; yu; repmat(params.H,1,N)];
        Vuav_new = [vx; vy];

        [SEE_m, g_ku, g_ub, ~, ~, Pf_m] = ...
            see_compute(P_k, zeta, P_b, P_u, Quav_new, Vuav_new, params);
        Pflight_approx(outer,m) = mean(Pf_m / params.P_lim);
        ASR_approx(outer,m)     = Pflight_approx(outer,m) * SEE_m;
        lambda_frac(outer,m)    = SEE_m;
        psi_gap(outer,m-1) = ASR_approx(outer,m-1) - ...
            lambda_frac(outer,m) * Pflight_approx(outer,m-1);
        SEE(l) = SEE_m;

        if abs(psi_gap(outer,m-1)) <= params.eps_dinkelbch || psi_gap(outer,m-1) >= 0
            fprintf('  (P4) Dinkelbach converged at inner iter %d\n', m-1);
            Q_uav = Quav_new;  V_uav = Vuav_new;
            break;
        else
            Q_state = [xu', yu', vx', vy', mu', un', sn', wn', zeros(N,3)];
        end
    end
    fprintf('  (P4) solved in %.2f s\n', toc(tStart));

    [SEE(l), g_ku, g_ub, ~, ~, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    F_algo = (SEE(l) - SEE(l-1)) / SEE(l-1);
    fprintf('  Fractional MSEE increase: %.6f\n', F_algo);
    if F_algo <= params.eps_bcd
        fprintf('SEE-FPow converged.\n');
        break;
    end
end

end

% =========================================================================

function [Itr, SEE] = run_asr_seq(params)
% RUN_ASR_SEQ  Benchmark: ASR maximization without flight power in objective.
%
%   [Itr, SEE] = run_asr_seq(params)
%
%   Implements the "ASR-Seq" benchmark from Section IV: a sequential BCD
%   algorithm maximizing the minimum average secrecy rate (ASR) directly,
%   without the propulsion power in the denominator.
%
%   This benchmark uses the same power optimizers (P1)-(P3) but with
%   lambda_n = 1/(2*ln2) (no flight power normalization), and the ASR
%   trajectory optimizer (trajectory_optim_asr) for the trajectory step.
%
% See also: run_msee_seqbcd.

N = params.N;
K = params.K;

[Q_uav, V_uav, xu, yu, vx, vy, P_k, P_u, P_b, zeta, ~] = feasible_init(params);

l = 1;
[SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
    see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

fprintf('ASR-Seq (no flight power) benchmark started.\n');

epsilon = 1e-7;
un = zeros(K, N);
for k = 1:K
    un(k,:) = norms([Q_uav - params.Qa(:,k)]) - epsilon;
end
k0 = (P_u + P_b) ./ (P_k .* P_u);
k1 = 1 ./ P_u;  k2 = P_k;  k3 = P_b;
sn = 1 ./ g_ku + epsilon;
wn = 1 ./ g_ub + epsilon;
Q_asr_state = [xu', yu', vx', vy', un', sn', wn', zeros(N, 1)];

while true
    l = l + 1;

    % (P1): lambda_n = 1/(2*ln2), no P_lim normalization
    lambda_n = 1 / (2*log(2));
    B_kn = g_ku ./ (P_b .* g_ub + 1);
    C_n  = P_u  .* g_ub;
    D_kn = (g_ub .* (P_u + P_b) + 1) ./ g_ku;
    [P_k, zeta, ~] = usr_sched_pow_optim( ...
        P_k, zeta, params.Pa_peak, params.Pa_ave, lambda_n, B_kn, C_n, D_kn);
    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);
    l = l + 1;

    % (P2)
    lambda_kn = zeta ./ (2*log(2));
    E_kn = P_k .* g_ku;
    F_kn = (P_k .* g_ku + P_b .* g_ub + 1) ./ g_ub;
    G_k  = sum(zeta .* log2(1 + snr_uav), 2) / 2;
    [P_u, ~] = uav_relay_pow_optim(P_u, params.Pu_peak, params.Pu_ave, ...
        lambda_kn, E_kn, F_kn, G_k, zeta);
    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);
    l = l + 1;

    % (P3)
    lambda_kn = zeta ./ (2*log(2));
    H_kn = g_ku .* P_k .* P_u;
    I_kn = (P_u .* g_ub + P_k .* g_ku + 1) ./ g_ub;
    J_kn = P_k .* g_ku ./ g_ub;
    K_n  = 1 ./ g_ub;
    [P_b, ~] = bs_jamming_pow_optim(P_b, params.Pb_peak, params.Pb_ave, ...
        lambda_kn, H_kn, I_kn, J_kn, K_n, zeta);
    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);
    l = l + 1;

    % (P4): ASR trajectory optimizer (no flight power in objective)
    for k = 1:K
        un(k,:) = norms([Q_uav - params.Qa(:,k)]) - epsilon;
    end
    k0 = (P_u + P_b) ./ (P_k .* P_u);
    k1 = 1 ./ P_u;  k2 = P_k;  k3 = P_b;
    sn = 1 ./ g_ku + epsilon;
    wn = 1 ./ g_ub + epsilon;
    Q_asr_state = [xu', yu', vx', vy', un', sn', wn', zeros(N,1)];

    Q_new = trajectory_optim_asr(Q_asr_state, params.H, params.Vmax, params.N0, ...
        params.beta0, params.af, params.Qa, params.Qb, params.QI(1:2), ...
        params.dt, params.R2, k0, k1, k2, k3, zeta);
    xu = Q_new(:,1)';  yu = Q_new(:,2)';
    vx = Q_new(:,3)';  vy = Q_new(:,4)';
    Q_uav = [xu; yu; repmat(params.H,1,N)];
    V_uav = [vx; vy];

    Q_asr_state = [xu', yu', vx', vy', un', sn', wn', zeros(N,1)];

    [SEE(l), g_ku, g_ub, ~, snr_uav, P_flight] = ...
        see_compute(P_k, zeta, P_b, P_u, Q_uav, V_uav, params);
    Itr(l) = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight);

    % ASR convergence check (use ASR = SEE * avg_normalized_Pflight)
    ASR_l     = SEE(l) * mean(P_flight / params.P_lim);
    ASR_l_prev = SEE(l-4) * mean(Itr(l-4).prop / params.P_lim);
    F_algo = (ASR_l - ASR_l_prev) / ASR_l_prev;
    fprintf('  Fractional ASR increase: %.6f\n', F_algo);
    if F_algo <= params.eps_bcd
        fprintf('ASR-Seq converged.\n');
        break;
    end
end

end

% =========================================================================
function s = pack_itr(P_k, P_b, P_u, zeta, Q_uav, V_uav, P_flight)
s.usrPow = P_k';  s.bsPow = P_b';  s.uavPow = P_u';
s.usrSch = zeta'; s.Trj = Q_uav';  s.Vel = V_uav';  s.prop = P_flight';
end
