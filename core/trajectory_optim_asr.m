function Q_new = trajectory_optim_asr(Q_old, H, Vmax, N0, ...
    beta0, af, Qa, Qb, QI, dt, R2, k0, k1, k2, k3, zeta)
% TRAJECTORY_OPTIM_ASR  Trajectory optimization for the ASR benchmark.
%
%   Q_new = trajectory_optim_asr(Q_old, H, Vmax, N0, beta0, af, Qa,
%       Qb, QI, dt, R2, k0, k1, k2, k3, zeta)
%
%   Solves the trajectory optimization subproblem for the ASR-Seq benchmark
%   (Section IV, "ASR-Seq") which maximizes the minimum average secrecy
%   rate WITHOUT accounting for the propulsion power in the denominator
%   (i.e., no average flight power constraint C9).
%
%   This benchmark corresponds to a standard max-min ASR design that
%   ignores UAV energy efficiency. It is used in the paper to demonstrate
%   the gain of the proposed MSEE formulation over pure secrecy-rate
%   optimization (Figures 2 and 4 of the paper).
%
%   The formulation is identical to trajectory_optim except:
%   - The Dinkelbach lambda parameter is absent (pure ASR, not SEE).
%   - Constraint C9 (average flight power <= P_lim) is removed.
%   - The mu_opt variable and its associated constraint (Eq. 42) are
%     removed since flight power does not appear in the objective.
%
%   INPUTS / OUTPUTS: see trajectory_optim for full documentation.
%   The difference is that this function maximizes only psi (the ASR
%   lower bound), not psi - lambda * P_bar_f.
%
% See also: trajectory_optim, run_asr_seq.

[K, N] = size(zeta);

xu  = Q_old(:, 1)';
yu  = Q_old(:, 2)';
vx  = Q_old(:, 3)';
vy  = Q_old(:, 4)';
un  = Q_old(:, 5 : 4+K)';
snu = Q_old(:, 5+K : 4+2*K)';
wn  = Q_old(:, 5+2*K)';

fc = flight_constants();

cvx_begin quiet
    variables xu_opt(1,N) yu_opt(1,N)
    variables vx_opt(1,N) vy_opt(1,N)
    variable  u_opt(K, N)  nonnegative
    variable  tt_opt(1, N) nonnegative
    variable  uu_opt(1, N) nonnegative
    variable  ttt_opt(K, N) nonnegative
    variable  uuu_opt(K, N) nonnegative
    variable  w_opt(1, N)   nonnegative
    variable  sl_opt(K, N)  nonnegative
    variable  su_opt(K, N)  nonnegative
    variable  y_opt(K, N)   nonnegative
    variable  y1_opt(K, N)  nonnegative
    variable  y2_opt(1, N)  nonnegative
    variable  psi

    maximize(psi)   % Maximize ASR lower bound (no flight-power term)

    subject to
        % Mobility constraints (C10-C14)
        xu_opt(1)   == QI(1);
        yu_opt(1)   == QI(2);
        xu_opt(N)   == QI(1);
        yu_opt(N)   == QI(2);
        xu_opt(2:N) - xu_opt(1:N-1) == vx_opt(1:N-1) * dt;
        yu_opt(2:N) - yu_opt(1:N-1) == vy_opt(1:N-1) * dt;
        norms([vx_opt; vy_opt]) <= Vmax;
        norms([vx_opt(2:N)-vx_opt(1:N-1); vy_opt(2:N)-vy_opt(1:N-1)]) <= fc.Amax;
        norms([xu_opt(2:N)-Qb(1); yu_opt(2:N)-Qb(2)]) <= R2;

        % BS inverse channel upper bound (Eq. 43f)
        tt_opt' >= square_pos(norms([repmat(H,1,N); xu_opt-Qb(1); yu_opt-Qb(2)])');
        uu_opt' + (beta0/N0) * rel_entr((N0/beta0)*tt_opt', w_opt') <= 0;
        uu_opt' >= af * pow_p(tt_opt', 3/2);

        % Secrecy rate SCA lower bound coefficients
        lambda_kn = zeta / (2 * log(2) * N);
        A0 = log(1 + 1./(k0.*snu + k1.*wn)) + log(1 + k3./wn) ...
             + (k0.*snu + k1.*wn) ./ ((k0.*snu + k1.*wn) .* (k0.*snu + k1.*wn + 1)) ...
             + k3 ./ (wn .* (k3 + wn));
        A1 = -k0 ./ ((k0.*snu + k1.*wn) .* (k0.*snu + k1.*wn + 1));
        A2 = -k1 ./ ((k0.*snu + k1.*wn) .* (k0.*snu + k1.*wn + 1)) ...
             - k3 ./ (wn .* (k3 + wn));

        sl_opt ./ k2 >= exp(-y1_opt);
        w_opt  ./ k3 >= exp(-y2_opt);

        for k = 1:K
            sum(lambda_kn(k,:)' .* (A0(k,:)' + A1(k,:)' .* su_opt(k,:)' ...
                + A2(k,:)' .* w_opt' - y_opt(k,:)')) >= psi;
            y_opt(k,:)' >= log_sum_exp([zeros(N,1), y1_opt(k,:)', y2_opt'], 2);

            ttt_opt(k,:)' >= square_pos( ...
                norms([xu_opt-Qa(1,k); yu_opt-Qa(2,k); repmat(H,1,N)])');;
            uuu_opt(k,:)' + (beta0/N0) * ...
                rel_entr((N0/beta0)*ttt_opt(k,:)', su_opt(k,:)') <= 0;
            uuu_opt(k,:)' >= af * pow_p(ttt_opt(k,:)', 3/2);

            (N0/beta0) * (un(k,:)'.^2 .* exp(af*un(k,:)') ...
                + un(k,:)' .* exp(af*un(k,:)') .* (af*un(k,:)'+2) ...
                .* (u_opt(k,:) - un(k,:))') >= sl_opt(k,:)';

            -(xu'.^2 + yu'.^2) + 2*(xu_opt)'.*(xu - Qa(1,k))' ...
                + 2*(yu_opt)'.*(yu - Qa(2,k))' ...
                + repmat(norms([H; Qa(1,k); Qa(2,k)])^2, N, 1) ...
                >= square_pos(u_opt(k,:)');
        end

cvx_end

if strcmp(cvx_status, 'Failed') || strcmp(cvx_status, 'Infeasible')
    warning('trajectory_optim_asr: CVX failed. Returning previous iterate.');
    Q_new = [Q_old, zeros(N, 1)];
    return;
end

Q_new = [full(xu_opt)', full(yu_opt)', full(vx_opt)', full(vy_opt)', ...
         full(u_opt)', full(su_opt)', full(w_opt)', repmat(psi, N, 1)];

end
