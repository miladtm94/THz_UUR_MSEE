function [P_b_new, obj_val] = bs_jamming_pow_optim( ...
    P_b, Pb_peak, Pb_ave, lambda_kn, H_kn, I_kn, J_kn, K_n, zeta)
% BS_JAMMING_POW_OPTIM  Solve subproblem (P3): BS cooperative jamming power.
%
%   [P_b_new, obj_val] = bs_jamming_pow_optim(
%       P_b, Pb_peak, Pb_ave, lambda_kn, H_kn, I_kn, J_kn, K_n, zeta)
%
%   Solves the SCA-approximated subproblem (P3.1) from Section III-C,
%   which optimizes the BS cooperative jamming power schedule P_b[n].
%
%   The BS jamming signal (sent in Phase 1 of DACJ) simultaneously degrades
%   the UUR's wiretap ability while contributing to the noise at UUR.
%   Increasing P_b improves physical-layer security but reduces the
%   relay's received SNR (both effects captured in the secrecy rate).
%
%   SUBPROBLEM (P3) CONSTRAINT (Eq. 32b):
%       sum_n lambda_kn(k,n) * [ ln(1 + H_kn/(P_b+I_kn))
%                               - ln(1 + J_kn/(P_b+K_n)) ] >= psi,  forall k
%
%   This constraint is of the form "convex minus convex" in P_b.
%   The first log term is concave (Lemma 2), so the SCA linearizes it at
%   the current iterate P_b^(l) to obtain the convex lower bound:
%       f3(P_b; P_b^(l), H_kn, I_kn)  [Eq. 34]
%
%   The second log term ln(1 + J_kn/(P_b+K_n)) is convex in P_b and is
%   kept as-is (subtracted), yielding a concave overall approximation.
%
%   INPUTS:
%       P_b      : Current BS jamming power [1xN] (warm-start)
%       Pb_peak  : Peak jamming power [W]   (C8)
%       Pb_ave   : Average jamming power [W] (C7)
%       lambda_kn: Scaling factor [KxN]
%                  lambda_kn(k,n) = P_lim * zeta(k,n) / sum(P_flight) / (2*ln2)
%       H_kn     : End-to-end SNR numerator [KxN]
%                  H_kn(k,n) = P_k(k,n) * g_ku(k,n) * P_u(n)
%       I_kn     : End-to-end SNR denominator offset [KxN]
%                  I_kn(k,n) = (P_u*g_ub + P_k*g_ku + 1) / g_ub(n)
%       J_kn     : Wiretap SNR numerator [KxN]
%                  J_kn(k,n) = P_k(k,n) * g_ku(k,n) / g_ub(n)
%       K_n      : Wiretap SNR constant [1xN]
%                  K_n(n) = 1 / g_ub(n)
%       zeta     : Current scheduling variables [KxN]
%
%   OUTPUTS:
%       P_b_new  : Optimized BS jamming power schedule [1xN]
%       obj_val  : Achieved MSEE lower bound value
%
% See also: usr_sched_pow_optim, uav_relay_pow_optim.

[K, N] = size(zeta);

cvx_begin quiet
    variable P_b_opt(1, N) nonnegative
    variable psi
    variables t_aux(K, N) u_aux(K, N)

    maximize(psi)

    subject to
        for k = 1:K
            % SCA lower bound of the concave first term (Eq. 34):
            %   f3(P_b; P_b^l, H_kn, I_kn)
            %   = ln(1 + H_kn/(P_b^l + I_kn))
            %     - H_kn / ((P_b^l+I_kn)*(P_b^l+H_kn+I_kn)) * (P_b - P_b^l)
            f3_lb = log(1 + H_kn(k,:) ./ (P_b + I_kn(k,:))) ...
                  - H_kn(k,:) ./ ((P_b + I_kn(k,:)) .* (P_b + H_kn(k,:) + I_kn(k,:))) ...
                    .* (P_b_opt - P_b);

            % y(k,n) = K_n/J_kn + P_b_opt/J_kn  (arg of ln(1+J/P_b+K))
            % Reformulate convex term -ln(1+J/y) as:
            %   -ln(1+J/(P_b+K)) = -(t_aux) where exp(-t_aux) <= 1 - u_aux
            %                                and u_aux >= 1/(y+1)
            y_n = K_n ./ J_kn(k,:) + P_b_opt ./ J_kn(k,:);
            sum(lambda_kn(k,:) .* (f3_lb - t_aux(k,:))) >= psi;
            exp(-t_aux(k,:)) <= 1 - u_aux(k,:);
            u_aux(k,:) >= inv_pos(y_n + 1);
        end

        % Power constraints (C7, C8)
        0 <= P_b_opt <= Pb_peak;
        sum(P_b_opt) / N <= Pb_ave;

cvx_end

if strcmp(cvx_status, 'Failed') || strcmp(cvx_status, 'Infeasible')
    warning('bs_jamming_pow_optim: CVX solver failed. Returning previous iterate.');
    P_b_new = P_b;
    obj_val = -Inf;
    return;
end

P_b_new = full(P_b_opt);
obj_val  = psi;

end
