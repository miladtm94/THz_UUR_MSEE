function [SEE, g_ku, g_ub, snr_bs, snr_uav, P_flight] = see_compute( ...
    P_k, zeta, P_b, P_u, Q_uav, V_uav, params)
% SEE_COMPUTE  Evaluate the minimum SEE and all intermediate metrics.
%
%   [SEE, g_ku, g_ub, snr_bs, snr_uav, P_flight] = see_compute(
%       P_k, zeta, P_b, P_u, Q_uav, V_uav, params)
%
%   Computes the minimum Secrecy Energy Efficiency (MSEE) over all K UEs,
%   as defined in Eq. (20) of the paper, together with the channel gains,
%   per-UE instantaneous SNRs, and propulsion power.
%
%   INPUTS:
%       P_k     : UE transmit power schedule [KxN], P_k(k,n) is the power
%                 of UE k at time slot n [W]
%       zeta    : User scheduling variables [KxN], zeta(k,n) in [0,1]
%                 (relaxed binary; 1 = UE k scheduled at slot n)
%       P_b     : BS jamming power schedule [1xN] [W]
%       P_u     : UAV relay transmit power [1xN] [W]
%       Q_uav   : UAV 3D position [3xN], Q_uav = [xu; yu; H*ones(1,N)]
%       V_uav   : UAV velocity [2xN], V_uav = [vx; vy]
%       params  : System parameter struct from system_params()
%
%   OUTPUTS:
%       SEE      : Minimum SEE over all K UEs [bits/Joule] normalized by
%                  (B, P_lim) per Remark 2 in the paper
%       g_ku     : THz channel gains UE-k to UUR [KxN] (Eq. 1)
%       g_ub     : THz channel gain UUR to BS [1xN]
%       snr_bs   : Effective end-to-end SNR at BS (Eq. 17 argument) [KxN]
%       snr_uav  : Wiretap SNR at UUR (Eq. 18 argument) [KxN]
%       P_flight : Instantaneous propulsion power [1xN] [W]
%
%   SEE DEFINITION (Eq. 20, normalized by B and P_lim per Remark 2):
%       SEEk = (1/N) * sum_n [ (1/2) * zeta_k[n] *
%              max(0, log2(1+snr_bs_k[n]) - log2(1+snr_uav_k[n])) ]
%              / mean(P_f[n] / P_lim)
%       MSEE = min_k SEEk
%
%   CHANNEL MODEL (Eq. 1-2):
%       g_ku[n] = beta0 * exp(-af * d_ku[n]) / (d_ku[n]^2 * N0)
%       where d_ku[n] = sqrt(||q[n]-qk||^2 + H^2)
%
%   NOTE:
%       The (.)^+ operator in the ASR (Eq. 19) is applied via max(0,...).
%       The factor 1/2 accounts for two-phase transmission (Eq. 19).
%
% See also: system_params, flight_power, usr_sched_pow_optim.

P_lim = params.P_lim;
Qa    = params.Qa;
Qb    = params.Qb;
af    = params.af;
beta0 = params.beta0;
N0    = params.N0;
K     = size(zeta, 1);

% ---- THz channel gain: UUR-to-BS link (Eq. 1 applied to UUR-BS) ----
d_ub = vecnorm(Qb - Q_uav);                              % [1xN]
g_ub = beta0 * exp(-af * d_ub) ./ (d_ub.^2 * N0);       % [1xN]

% ---- Per-UE THz channel gains and SNRs ----
g_ku   = zeros(K, size(zeta, 2));
snr_bs = zeros(K, size(zeta, 2));
snr_uav = zeros(K, size(zeta, 2));
R_sec  = zeros(size(zeta));

for k = 1:K
    % UE-k to UUR distance and channel gain (Eq. 1-2)
    d_ku      = vecnorm(Q_uav - Qa(:, k));               % [1xN]
    g_ku(k,:) = beta0 * exp(-af * d_ku) ./ (d_ku.^2 * N0);

    % End-to-end SNR at BS for the AF relay link (Eq. 17, argument of log2)
    % gamma_b^k[n] = P_k * g_ku * P_u * g_ub
    %                / ((P_u + P_b) * g_ub + P_k * g_ku + 1)
    snr_bs(k,:) = (P_k(k,:) .* P_u .* g_ku(k,:) .* g_ub) ./ ...
                  ((P_u + P_b) .* g_ub + P_k(k,:) .* g_ku(k,:) + 1);

    % Wiretap SNR at UUR (Eq. 18, argument of log2)
    % gamma_u^k[n] = P_k * g_ku / (P_b * g_ub + 1)
    snr_uav(k,:) = (P_k(k,:) .* g_ku(k,:)) ./ (P_b .* g_ub + 1);

    % Instantaneous secrecy rate (Eq. 19): (1/2) * zeta_k * max(0, R_b - R_u)
    R_sec(k,:) = 0.5 * zeta(k,:) .* ...
                 max(0, log2(1 + snr_bs(k,:)) - log2(1 + snr_uav(k,:)));
end

% ---- Propulsion power (Eq. 11) ----
P_flight = flight_power(V_uav);  % [1xN]

% ---- MSEE computation (Eq. 20, normalized by B and P_lim per Remark 2) ----
% ASR of UE k = (1/N) * sum_n R_sec(k,n) [bits/s/Hz]
ASR_per_ue = mean(R_sec, 2);             % [Kx1]
MSEE_unnorm = min(ASR_per_ue);           % min over k
% Denominator: average normalized flight power
avg_norm_flight_pow = mean(P_flight / P_lim);

SEE = MSEE_unnorm / avg_norm_flight_pow;

end
