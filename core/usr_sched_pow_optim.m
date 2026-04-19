function [P_k_new, zeta_new, obj_val] = usr_sched_pow_optim( ...
    P_k, zeta, Pa_peak, Pa_ave, lambda_n, B_kn, C_n, D_kn)
% USR_SCHED_POW_OPTIM  Solve subproblem (P1): joint UE scheduling & power.
%
%   [P_k_new, zeta_new, obj_val] = usr_sched_pow_optim(
%       P_k, zeta, Pa_peak, Pa_ave, lambda_n, B_kn, C_n, D_kn)
%
%   Solves the convex approximation of the joint user scheduling and UE
%   transmit power subproblem (P1.2) from Section III-A of the paper.
%
%   The binary scheduling variables zeta_k[n] are first relaxed to
%   continuous values in [0,1] (time-sharing relaxation), and then the
%   Penalty-SCA (PSCA) technique is used to force them back toward 0/1
%   via a penalty on the bilinear term zeta*(1-zeta) (Eq. 25).
%
%   CHANGE OF VARIABLES (Section III-A):
%       p_tilde_k[n] = p_k[n] * zeta_k[n]
%   transforms the bilinear coupling of P_k and zeta into a product of
%   two jointly optimized variables, enabling the use of CVX.
%
%   SCA APPROXIMATION (Lemma 1, Eq. 23):
%       Term II in (22a) is the concave function Z2(zeta, p_tilde; B_kn).
%       Its linear upper bound f1^ub is used to make (22a) convex.
%
%   CVX REFORMULATION (Remark 3, Eq. 26-28):
%       The concave function Z1 is expressed via relative entropy functions
%       Erel(x,y) = x*log(x/y) which are natively supported by CVX's
%       disciplined convex programming (DCP) rules.
%
%   PENALTY-SCA (Eq. 25):
%       Objective: maximize  psi - mu * eta
%       Constraint: sum(sum((1-2*zeta^(l))*zeta + (zeta^(l))^2)) <= eta
%       The penalty weight mu is doubled at each outer iteration until
%       eta < eta_min (binary constraint approximately satisfied).
%
%   INPUTS:
%       P_k      : Current UE power [KxN], from previous BCD iteration
%       zeta     : Current scheduling [KxN], continuous in [0,1]
%       Pa_peak  : UE peak power constraint [W]   (C4)
%       Pa_ave   : UE average power constraint [W] (C3)
%       lambda_n : Scaling factor = P_lim / sum(P_flight) / (2*ln2)  [1xN]
%       B_kn     : Eavesdropping SNR term [KxN]
%                  B_kn(k,n) = g_ku(k,n) / (P_b(n)*g_ub(n) + 1)
%       C_n      : Relay-link numerator term [1xN]
%                  C_n(n) = P_u(n) * g_ub(n)
%       D_kn     : Relay-link denominator term [KxN]
%                  D_kn(k,n) = (g_ub(n)*(P_u(n)+P_b(n))+1) / g_ku(k,n)
%
%   OUTPUTS:
%       P_k_new  : Optimized UE power schedule [KxN]
%       zeta_new : Optimized scheduling variables [KxN] in [0,1]
%       obj_val  : Achieved value of the MSEE lower bound (psi)
%
% See also: uav_relay_pow_optim, bs_jamming_pow_optim, trajectory_optim.

[K, N] = size(zeta);

% Penalty-SCA parameters
mu       = 0.001;    % Initial penalty weight (small to allow initial infeasibility)
eta_min  = 1e-2;     % Convergence threshold for binary constraint violation
max_iter = 50;       % Safety cap on penalty iterations

% Current product variable: p_tilde = p_k * zeta
P_tilde = P_k .* zeta;

for iter = 1:max_iter
    cvx_begin quiet
        variable P_tilde_opt(K, N) nonnegative
        variable zeta_opt(K, N) nonnegative
        variables psi eta

        maximize (psi - mu * eta)

        subject to
            % --- Secrecy rate lower bound constraint (P1.2, Eq. 25a) ---
            % Uses CVX relative-entropy reformulation of Z1 (Remark 3, Eq. 26)
            % and SCA linear upper bound of Z2 (Lemma 1, Eq. 23)
            for k = 1:K
                a_n = C_n;
                b_n = D_kn(k, :);

                % Z1 reformulated via rel_entr (Eq. 26):
                %   Z1 = -(1+a)/(a*b) * Erel(p+b*z, (a+1)*p+b*z)
                %        - 1/(a*b)     * Erel((a+1)*p+b*z, p+b*z)
                x_n = P_tilde_opt(k,:) + b_n .* zeta_opt(k,:);    % p_tilde + b*zeta
                y_n = (a_n + 1) .* P_tilde_opt(k,:) + b_n .* zeta_opt(k,:);

                % Linear upper bound of Z2 at current iterate (Lemma 1, Eq. 23):
                %   f1^ub(zeta, p_tilde; zeta^l, p_tilde^l, B_kn)
                f0_z2  = lambda_n .* zeta(k,:) .* log(1 + B_kn(k,:) .* P_k(k,:));
                df_dz  = lambda_n .* (log(1 + B_kn(k,:) .* P_k(k,:)) - ...
                          B_kn(k,:) .* P_k(k,:) ./ (1 + B_kn(k,:) .* P_k(k,:)));
                df_dp  = lambda_n .* B_kn(k,:) ./ (1 + B_kn(k,:) .* P_k(k,:));

                sum(lambda_n .* ( -(a_n + 1)./(a_n .* b_n) .* rel_entr(x_n, y_n) ...
                                  - 1./(a_n .* b_n) .* rel_entr(y_n, x_n) )) ...
                    - ( f0_z2 * ones(N,1) ...
                        + df_dz * (zeta_opt(k,:) - zeta(k,:))' ...
                        + df_dp * (P_tilde_opt(k,:) - P_tilde(k,:))' ) >= psi;
            end

            % --- Linearized binary constraint (P1.2, Eq. 25b) ---
            % First-order Taylor of zeta^2 at current iterate zeta^(l)
            sum(sum((1 - 2 * zeta) .* zeta_opt + zeta.^2)) <= eta;

            % --- Power constraints (C3, C4 in terms of p_tilde) ---
            sum(P_tilde_opt, 2) / N <= Pa_ave;          % Avg power (C3)
            0 <= P_tilde_opt <= zeta_opt * Pa_peak;     % Peak power (C4)

            % --- Scheduling constraints (C1, C2 relaxed) ---
            0 <= zeta_opt <= 1;                         % Relaxed binary (C1)
            (ones(1, K) * zeta_opt)' <= 1;             % At most one UE per slot (C2)

    cvx_end

    if strcmp(cvx_status, 'Failed') || strcmp(cvx_status, 'Infeasible')
        warning('usr_sched_pow_optim: CVX failed at penalty iter %d.', iter);
        break;
    end

    zeta   = max(0, full(zeta_opt));
    P_tilde = max(0, full(P_tilde_opt));

    mu = mu * 2;    % Double penalty weight to enforce binary constraint
    if eta <= eta_min
        break;
    end
end

% Recover P_k from p_tilde = P_k * zeta (guard against division by zero)
safe_zeta = max(zeta, 1e-10);
P_k_new  = P_tilde ./ safe_zeta;
zeta_new = zeta;
obj_val  = psi;

end
