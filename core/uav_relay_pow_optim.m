function [P_u_new, obj_val] = uav_relay_pow_optim( ...
    P_u, Pu_peak, Pu_ave, lambda_kn, E_kn, F_kn, G_k, zeta)
% UAV_RELAY_POW_OPTIM  Solve subproblem (P2): UUR relay transmit power.
%
%   [P_u_new, obj_val] = uav_relay_pow_optim(
%       P_u, Pu_peak, Pu_ave, lambda_kn, E_kn, F_kn, G_k, zeta)
%
%   Solves the convex subproblem (P2) from Section III-B of the paper,
%   which optimizes the UAV relay transmit power P_u[n] for all time slots.
%
%   SUBPROBLEM (P2) (Eq. 29):
%       maximize  psi
%       s.t.  sum_n lambda_kn(k,n) * ln(1 + E_kn*P_u/(P_u+F_kn)) - G_k >= psi, forall k
%             (1/N) * sum_n P_u[n] <= Pu_ave
%             0 <= P_u[n] <= Pu_peak
%
%   CONVEXITY (Lemma 2, Eq. 30):
%       The function f2(P_u; E, F) = ln(1 + E*P_u/(P_u+F)) is concave when
%       E >= 1 (i.e., E_kn = P_k*g_ku > 0 and F_kn = (P_k*g_ku+P_b*g_ub+1)/g_ub > 0).
%       This follows from the concavity of ln(1 + q*x) with q >= 0 composed
%       with the concave (when ad >= bc) function (ax+b)/(cx+d) (Lemma 2).
%       Hence subproblem (P2) is convex and solvable directly via CVX.
%
%   CVX REFORMULATION:
%       The ln(1 + (ax+b)/(cx+d)) form is rewritten to comply with CVX DCP
%       rules using auxiliary variables and inequality constraints as:
%           (Pu+bn)*(1+an)/(an*bn) >= x_temp(k,n)  (epigraph of 1/P_u)
%           exp(y_temp) <= 1 - u_temp               (log-exp cone)
%           u_temp >= 1 / x_temp                     (perspective of 1/x)
%
%   INPUTS:
%       P_u      : Current UAV relay power [1xN] (warm-start)
%       Pu_peak  : Peak relay power constraint [W]  (C6)
%       Pu_ave   : Average relay power constraint [W] (C5)
%       lambda_kn: Scaling factor [KxN]
%                  lambda_kn(k,n) = P_lim * zeta(k,n) / sum(P_flight) / (2*ln2)
%       E_kn     : Relay SNR numerator term [KxN]
%                  E_kn(k,n) = P_k(k,n) * g_ku(k,n)
%       F_kn     : Relay SNR denominator term [KxN]
%                  F_kn(k,n) = (E_kn + P_b*g_ub + 1) / g_ub(n)
%       G_k      : Eavesdropping rate offset [Kx1]
%                  G_k(k) = sum_n[ P_lim*zeta(k,n)*log2(1+snr_uav(k,n)) ]
%                           / sum(P_flight) / 2
%       zeta     : Current scheduling variables [KxN]
%
%   OUTPUTS:
%       P_u_new  : Optimized UAV relay power schedule [1xN]
%       obj_val  : Achieved MSEE lower bound value
%
% See also: usr_sched_pow_optim, bs_jamming_pow_optim.

[K, N] = size(zeta);

cvx_begin quiet
    variable P_u_opt(1, N) nonnegative
    variable x_temp(K, N)
    variable y_temp(K, N)
    variable u_temp(K, N)
    variable psi

    maximize(psi)

    subject to
        for k = 1:K
            a_n = E_kn(k, :);    % P_k * g_ku [1xN]
            b_n = F_kn(k, :);    % (E_kn + P_b*g_ub + 1) / g_ub [1xN]

            % Reformulate ln(1 + a*P_u/(P_u+b)) via auxiliary variables
            % to comply with CVX DCP rules (see Lemma 2 reformulation):
            %   ln(1 + a*x/(x+b)) = ln((1+a)*(x+b/(1+a)) / (x+b))
            %                     = ln(1+a) + ln(1 - 1/((x+b)*(1+a)/b))
            % which becomes concave in x and tractable via log1p epigraph
            (P_u_opt + b_n) .* (1 + a_n) ./ (a_n .* b_n) >= x_temp(k, :);
            sum(lambda_kn(k,:) .* (log(1 + a_n) + y_temp(k,:))) - G_k(k) >= psi;
            exp(y_temp(k,:))  <= 1 - u_temp(k,:);
            u_temp(k,:)       >= inv_pos(x_temp(k,:));
            x_temp(k,:)       >= 1;
        end

        % Power constraints (C5, C6)
        0 <= P_u_opt <= Pu_peak;
        sum(P_u_opt) / N <= Pu_ave;

cvx_end

if strcmp(cvx_status, 'Failed') || strcmp(cvx_status, 'Infeasible')
    warning('uav_relay_pow_optim: CVX solver failed. Returning previous iterate.');
    P_u_new = P_u;
    obj_val = -Inf;
    return;
end

P_u_new = full(P_u_opt);
obj_val = psi;

end
