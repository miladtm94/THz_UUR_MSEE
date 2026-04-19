function Q_new = trajectory_optim(Q_old, H, Vmax, P_lim, N0, ...
    beta0, af, Qa, Qb, QI, dt, R1, R2, k0, k1, k2, k3, zeta, lambda_dinkl)
% TRAJECTORY_OPTIM  Solve one Dinkelbach iteration of subproblem (P4).
%
%   Q_new = trajectory_optim(Q_old, H, Vmax, P_lim, N0, beta0, af, Qa,
%       Qb, QI, dt, R1, R2, k0, k1, k2, k3, zeta, lambda_dinkl)
%
%   Solves one convex approximation (P4.5) from Section III-D for the
%   joint UAV trajectory and velocity optimization within the Dinkelbach
%   fractional-programming framework (Algorithm 1).
%
%   DINKELBACH REFORMULATION (Eq. 37-38):
%   The fractional objective min_k R_bar_k / P_bar_f is converted to
%   an equivalent parametric subtractive form (Eq. 38):
%       maximize  min_k R_bar_k(q,v) - lambda^(m) * P_bar_f(v)
%   where lambda^(m) = min_k R_bar_k / P_bar_f at the m-th iteration.
%   Convergence to lambda* is guaranteed (Dinkelbach, [47]).
%
%   CONVEX APPROXIMATIONS APPLIED (Lemma 3, Eqs. 44-47):
%   1. Flight power upper bound P_f^ub (Eq. 41): replaces the non-convex
%      induced power term using a slack variable mu[n] (Eq. 39).
%   2. Constraint on mu[n] (Eq. 40e linearized via SCA, Eq. 42):
%      Taylor expansion of the convex LHS at (mu^m, v^m).
%   3. Secrecy rate lower bound: slack variables (r, w, s, u) introduced
%      to decouple UAV position from the logarithmic expressions.
%      Constraints (43d)-(43f) enforce the THz channel variable bounds.
%   4. SCA for the non-convex distance constraint (Eq. 49c):
%      Taylor expansion of ||q-qk||^2 at q^m.
%   5. CVX-incompatible constraints handled via Remark 4 reformulations
%      using log-sum-exp, relative entropy, and exponential cone proxies.
%
%   INPUTS:
%       Q_old       : Current iterate state vector [N x (5+3K+1)]
%                     Columns: [xu, yu, vx, vy, mu, un(K), sn(K), wn, ...]
%       H           : UAV flight altitude [m]
%       Vmax        : Maximum UAV speed [m/s] (C12)
%       P_lim       : Average flight power budget [W] (C9)
%       N0          : Noise power [W]
%       beta0       : Reference channel gain at 1m
%       af          : Molecular absorption coefficient [1/m]
%       Qa          : UE 3D positions [3xK]
%       Qb          : BS 3D position [3x1]
%       QI          : UAV initial/final position [x; y] [2x1]
%       dt          : Time slot duration [s]
%       R1          : Inner boundary radius (unused in current formulation)
%       R2          : Outer boundary radius (C14)
%       k0, k1, k2, k3 : Secrecy rate expression parameters [KxN] / [1xN]
%                      k0 = (P_u+P_b)/(P_k*P_u),  k1 = 1/P_u
%                      k2 = P_k,  k3 = P_b
%       zeta        : User scheduling [KxN]
%       lambda_dinkl: Current Dinkelbach parameter lambda^(m)
%
%   OUTPUT:
%       Q_new  : Updated state [N x (5+3K+1+3)] containing:
%                [xu, yu, vx, vy, mu, u_opt(K), su_opt(K), w_opt, psi, P_bar_f, obj]
%
%   NOTE:
%       The (1/v0^2) factor in the Taylor expansion of constraint (40e)
%       corresponds to the coefficient 1/(2*v0^2) in the paper's notation
%       (Eq. 42); this arises because the SCA takes the gradient of
%       ||v||^2/(2*v0^2). The code uses 1/v0^2 to be consistent with
%       the explicit Taylor derivation in Appendix B (Lemma 3).
%
% See also: uav_relay_pow_optim, bs_jamming_pow_optim, see_compute.

[K, N] = size(zeta);

% Unpack current iterate from Q_old
xu  = Q_old(:, 1)';   % [1xN] x-position
yu  = Q_old(:, 2)';   % [1xN] y-position
vx  = Q_old(:, 3)';   % [1xN] x-velocity
vy  = Q_old(:, 4)';   % [1xN] y-velocity
mu  = Q_old(:, 5)';   % [1xN] induced-power auxiliary variable (Eq. 39)
un  = Q_old(:, 6 : 5+K)';           % [KxN] distance lower-bound slack
snu = Q_old(:, 6+K : 5+2*K)';       % [KxN] channel gain upper-bound slack (s_k)
wn  = Q_old(:, 6+2*K)';             % [1xN] BS channel gain upper-bound slack (w)

fc = flight_constants();

cvx_begin quiet
    % ---- Decision variables ----
    variables xu_opt(1,N) yu_opt(1,N)       % UAV horizontal positions
    variables vx_opt(1,N) vy_opt(1,N)       % UAV velocities
    variable  mu_opt(1,N)  nonnegative      % Induced power slack (Eq. 39)
    variable  u_opt(K, N)  nonnegative      % Distance lower-bound [Kx N]
    variable  tt_opt(1, N) nonnegative      % BS distance^2 auxiliary
    variable  uu_opt(1, N) nonnegative      % BS channel exp-cone aux
    variable  ttt_opt(K, N) nonnegative     % UE distance^2 auxiliary
    variable  uuu_opt(K, N) nonnegative     % UE channel exp-cone aux
    variable  w_opt(1, N)  nonnegative      % Inverse BS channel gain (w[n])
    variable  sl_opt(K, N) nonnegative      % Lower bound on 1/g_ku (s_lb)
    variable  su_opt(K, N) nonnegative      % Upper bound on 1/g_ku (s_ub)
    variable  y_opt(K, N)  nonnegative      % log-sum-exp auxiliary
    variable  y1_opt(K, N) nonnegative      % log(k2/s) auxiliary
    variable  y2_opt(1, N) nonnegative      % log(k3/w) auxiliary
    variable  psi                           % MSEE lower bound (inner min)

    % ---- Objective: Dinkelbach subtractive form (Eq. 38) ----
    expression P_bar_f_ub   % Average normalized flight power upper bound
    P_bar_f_ub = (1/N) * sum( fc.P0 * (1 + 2 * square_pos(norms([vx_opt; vy_opt])) ...
                                         / (fc.Omega^2 * fc.R_rotor^2)) ...
                              + 0.5 * fc.d0 * fc.rho * fc.s * fc.A ...
                                * pow_pos(norms([vx_opt; vy_opt]), 3) ...
                              + fc.Pi * mu_opt ) / P_lim;

    maximize(psi - lambda_dinkl * P_bar_f_ub)

    subject to
        % ---- UAV mobility constraints (C10-C14, Eq. 13) ----
        xu_opt(1)   == QI(1);                      % Initial x position (C10)
        yu_opt(1)   == QI(2);                      % Initial y position (C10)
        xu_opt(N)   == QI(1);                      % Final x position   (C10)
        yu_opt(N)   == QI(2);                      % Final y position   (C10)
        xu_opt(2:N) - xu_opt(1:N-1) == vx_opt(1:N-1) * dt;  % Kinematics (C11)
        yu_opt(2:N) - yu_opt(1:N-1) == vy_opt(1:N-1) * dt;
        norms([vx_opt; vy_opt]) <= Vmax;           % Speed limit (C12)
        norms([vx_opt(2:N)-vx_opt(1:N-1); vy_opt(2:N)-vy_opt(1:N-1)]) <= fc.Amax; % (C13)
        norms([xu_opt(2:N)-Qb(1); yu_opt(2:N)-Qb(2)]) <= R2;  % Region (C14)

        % ---- Average flight power constraint (C9, Eq. 12) ----
        P_bar_f_ub <= 1;

        % ---- Induced power auxiliary constraint (Eq. 42, SCA of Eq. 40e) ----
        % Taylor expansion of LHS of (40e) at current iterate (mu^m, v^m):
        %   -mu^2 + 2*mu^m*mu + (1/v0^2)*(-||v^m||^2 + 2*(v^m)^T*v) >= 1/mu^2
        mu_opt >= 0;
        -mu'.^2 + 2*mu'.*mu_opt' + (1/fc.nu0^2) ...
            .* (-(vx'.^2 + vy'.^2) + 2*(vx'.*vx_opt' + vy'.*vy_opt')) ...
            >= pow_p(mu_opt', -2);

        % ---- BS-to-UUR inverse channel bound (Eq. 43f, Remark 4 handling) ----
        % Enforce: N0/beta0 * d_ub^2 * exp(af*d_ub) <= w_opt
        % Via Remark 4 (Eq. 50): introduce t1 = d_ub^2, t2 >= af*t1^(3/2)
        tt_opt' >= square_pos(norms([repmat(H,1,N); xu_opt-Qb(1); yu_opt-Qb(2)])');
        uu_opt' + (beta0/N0) * rel_entr((N0/beta0)*tt_opt', w_opt') <= 0;
        uu_opt' >= af * pow_p(tt_opt', 3/2);

        % ---- Secrecy rate lower bound (Eq. 43c, 48c via Lemma 3) ----
        % Linearized secrecy rate expression; SCA at (snu, wn) (Eq. 44-47)
        lambda_kn = zeta / (2 * log(2) * N);

        % Coefficients of the SCA linear lower bound at current iterate (snu, wn)
        A0 = log(1 + 1./(k0.*snu + k1.*wn)) + log(1 + k3./wn) ...
             + (k0.*snu + k1.*wn) ./ ((k0.*snu + k1.*wn) .* (k0.*snu + k1.*wn + 1)) ...
             + k3 ./ (wn .* (k3 + wn));
        A1 = -k0 ./ ((k0.*snu + k1.*wn) .* (k0.*snu + k1.*wn + 1));
        A2 = -k1 ./ ((k0.*snu + k1.*wn) .* (k0.*snu + k1.*wn + 1)) ...
             - k3 ./ (wn .* (k3 + wn));

        % Reformulate the term -log(1+k2/s+k3/w) via log-sum-exp (Remark 4, Eq. 51)
        %   sl_opt / k2 >= exp(-y1_opt)   <=>  -y1_opt >= log(k2/sl_opt)
        %   w_opt  / k3 >= exp(-y2_opt)   <=>  -y2_opt >= log(k3/w_opt)
        sl_opt ./ k2 >= exp(-y1_opt);
        w_opt  ./ k3 >= exp(-y2_opt);

        for k = 1:K
            % Secrecy rate lower bound constraint (Eq. 48c)
            sum(lambda_kn(k,:)' .* (A0(k,:)' + A1(k,:)' .* su_opt(k,:)' ...
                + A2(k,:)' .* w_opt' - y_opt(k,:)')) >= psi;

            % log-sum-exp constraint (Remark 4, Eq. 51)
            y_opt(k,:)' >= log_sum_exp([zeros(N,1), y1_opt(k,:)', y2_opt'], 2);

            % ---- UE-to-UUR inverse channel upper bound (Eq. 43e) ----
            % N0/beta0 * d_ku^2 * exp(af*d_ku) <= su_opt  (convex, kept as-is)
            ttt_opt(k,:)' >= square_pos( ...
                norms([xu_opt-Qa(1,k); yu_opt-Qa(2,k); repmat(H,1,N)])');;
            uuu_opt(k,:)' + (beta0/N0) * ...
                rel_entr((N0/beta0)*ttt_opt(k,:)', su_opt(k,:)') <= 0;
            uuu_opt(k,:)' >= af * pow_p(ttt_opt(k,:)', 3/2);

            % ---- UE-to-UUR inverse channel lower bound (Eq. 48d via Eq. 46) ----
            % N0/beta0 * d_ku^2 * exp(af*d_ku) >= sl_opt
            % SCA lower bound using Lemma 3 Eq. (46): f43^lb at un^m
            (N0/beta0) * (un(k,:)'.^2 .* exp(af*un(k,:)') ...
                + un(k,:)' .* exp(af*un(k,:)') .* (af*un(k,:)'+2) ...
                .* (u_opt(k,:) - un(k,:))') >= sl_opt(k,:)';

            % ---- SCA for distance constraint (Eq. 49c) ----
            % Taylor expansion of ||q-qk||^2 at q^m:
            %  -||q^m||^2 + 2*(q^m-qk)^T*q + ||qk||^2 + H^2 >= u_opt^2
            -( xu'.^2 + yu'.^2) + 2*(xu_opt)'.*(xu - Qa(1,k))' ...
                + 2*(yu_opt)'.*(yu - Qa(2,k))' ...
                + repmat(norms([H; Qa(1,k); Qa(2,k)])^2, N, 1) ...
                >= square_pos(u_opt(k,:)');
        end

cvx_end

if strcmp(cvx_status, 'Failed') || strcmp(cvx_status, 'Infeasible')
    warning('trajectory_optim: CVX failed. Returning previous iterate.');
    Q_new = [Q_old, zeros(N, 3)];
    return;
end

% Pack outputs: [position, velocity, mu, u, su, w, psi, P_bar_f_ub, obj_val]
Q_new = [full(xu_opt)', full(yu_opt)', full(vx_opt)', full(vy_opt)', ...
         full(mu_opt)', full(u_opt)', full(su_opt)', full(w_opt)', ...
         repmat(psi, N, 1), repmat(P_bar_f_ub, N, 1), ...
         repmat(psi - lambda_dinkl * P_bar_f_ub, N, 1)];

end
