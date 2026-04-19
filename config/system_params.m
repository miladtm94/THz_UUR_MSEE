function params = system_params()
% SYSTEM_PARAMS  Return all simulation parameters for the THz-UUR-MSEE system.
%
%   params = system_params() returns a struct containing all system-level
%   parameters used throughout the simulation, as defined in Section II of:
%
%   M. T. Mamaghani and Y. Hong, "Terahertz Meets Untrusted UAV-Relaying:
%   Minimum Secrecy Energy Efficiency Maximization via Trajectory and
%   Communication Co-design," IEEE Transactions on Communications, 2022.
%
%   OUTPUT:
%       params  - struct with the following fields:
%
%   --- THz Channel Parameters (Section II-A) ---
%       fc       : Carrier frequency [Hz]          (default: 0.8 THz)
%       BW       : Communication bandwidth [Hz]    (default: 10 GHz)
%       N0       : Noise power [W]                 (PSD x BW)
%       beta0    : Reference channel gain at 1 m   (Eq. 1)
%       af       : Molecular absorption coefficient [1/m] (default: 0.005)
%
%   --- Network Topology ---
%       K        : Number of ground UEs
%       H        : UAV flight altitude [m]
%       R1       : Inner boundary radius [m] (minimum UE-BS distance)
%       R2       : Outer boundary / max UAV range [m]
%       Qb       : BS horizontal coordinates [1x2], placed at origin
%       Qa       : UE 3D coordinates [3xK]
%       QI       : UAV initial/final 3D position [3x1]
%
%   --- Power Budget Constraints (Eqs. 5-10) ---
%       P_ave    : Overall average power reference [W]
%       Pa_ave   : UE average transmit power [W]   (0.1 * P_ave)
%       Pa_peak  : UE peak transmit power [W]      (4 * Pa_ave)
%       Pu_ave   : UAV relay average transmit power [W]  (0.4 * P_ave)
%       Pu_peak  : UAV relay peak transmit power [W]     (4 * Pu_ave)
%       Pb_ave   : BS jammer average transmit power [W]  (0.5 * P_ave)
%       Pb_peak  : BS jammer peak transmit power [W]     (4 * Pb_ave)
%       P_lim    : UAV average flight power budget [W]   (Eq. 12, C9)
%
%   --- UAV Flight Parameters (Eqs. 11-13, C10-C14) ---
%       Vmax     : Maximum UAV speed [m/s]
%       Amax     : Maximum UAV acceleration [m/s^2]
%       Ru       : Nominal circular trajectory radius [m]
%       N        : Number of time slots
%       dt       : Duration of each time slot [s] (T/N)
%
%   --- Algorithm Convergence Thresholds ---
%       eps_bcd      : BCD outer-loop fractional termination threshold
%       eps_dinkelbch: Dinkelbach inner-loop termination threshold
%
% NOTES:
%   - All power values are in Watts (linear scale).
%   - User positions are drawn uniformly in the annulus [R1, R2] using a
%     fixed random seed (rng default) for reproducibility (Section IV).
%   - The Consts cell array packs {P_lim, Qa, Qb, af, beta0, N0, BW}
%     for convenient passing to SEE evaluation functions.
%
% See also: flight_constants, ue_random_placement, feasible_init.

%% =====================================================================
%  THz Channel and Noise Parameters
% ======================================================================
C           = 3e8;               % Speed of light [m/s]
params.fc   = 0.8e12;            % Carrier frequency: 0.8 THz [Hz]
PSD         = 1e-3 * db2pow(-196); % Noise PSD: -174 dBm/Hz -> W/Hz
params.BW   = 10e9;              % THz bandwidth [Hz] (10 GHz sub-band)
params.N0   = PSD * params.BW;   % Thermal noise power [W]

% Reference path gain at 1 m: beta0 = (C / 4*pi*fc)^2  (Eq. 1)
params.beta0 = (C / (4 * pi * params.fc))^2;

% Molecular absorption coefficient [m^-1] (default case; paper sweeps 0.005-0.025)
params.af = 0.005;

%% =====================================================================
%  Transmit Power Budgets (Constraints C3-C8, Eqs. 5-10)
% ======================================================================
params.P_ave    = 1;                     % Baseline average network power [W]
params.Pa_ave   = 0.1 * params.P_ave;   % UE average TX power [W]
params.Pa_peak  = 4   * params.Pa_ave;  % UE peak TX power [W]
params.Pb_ave   = 0.5 * params.P_ave;   % BS jammer average TX power [W]
params.Pb_peak  = 4   * params.Pb_ave;  % BS jammer peak TX power [W]
params.Pu_ave   = 0.4 * params.P_ave;   % UAV relay average TX power [W]
params.Pu_peak  = 4   * params.Pu_ave;  % UAV relay peak TX power [W]
params.P_lim    = 200;                  % UAV average propulsion power limit [W] (C9)

%% =====================================================================
%  Network Topology (Section II, 3D Cartesian system)
% ======================================================================
params.H  = 10;    % Fixed UAV altitude [m] (Footnote 2: fixed-altitude operation)
params.R1 = 20;    % Annulus inner radius [m] (min reliable UE-to-BS distance)
params.R2 = 30;    % Annulus outer radius [m] (UAV permitted flying boundary)

% Base station at origin (Section II)
params.Qb = [0; 0; 0];

% UE placement: K users uniformly distributed in annulus [R1, R2]
params.K = 5;     % Number of ground UEs (K)
rng default;      % Fixed seed for reproducibility (Section IV)
[xa, ya] = ue_random_placement(0, 0, params.R1, params.R2, params.K);
params.Qa = [xa; ya; zeros(1, params.K)];  % UE 3D coordinates [3xK]

% UAV initial/final (depot) position: on the annulus midline (C10)
params.Ru = mean([params.R1, params.R2]);  % Default circular path radius
params.QI = [params.Ru; 0; params.H];      % [x; y; z] of depot

%% =====================================================================
%  UAV Mobility Parameters (Constraints C10-C14, Eq. 13)
% ======================================================================
% Flight constants (blade profile power, induced power, etc.) loaded below
fc_params = flight_constants();
params.Vmax = fc_params.Vmax;  % Maximum speed [m/s]
params.Amax = fc_params.Amax;  % Maximum acceleration [m/s^2]

% Time discretization: N slots of duration dt = T/N (Section II)
% N is set externally (default 100), but dt is derived from velocity condition
% dt = 2/Vmax ensures that the maximum displacement per slot << H
params.N  = 100;
params.dt = 2 / params.Vmax;   % Time slot duration [s]  (d_max_t << H condition)

%% =====================================================================
%  Algorithm Convergence Thresholds (Sections III-D, III-E)
% ======================================================================
params.eps_bcd       = 1e-3;  % BCD outer-loop fractional increase threshold (eps_1)
params.eps_dinkelbch = 1e-4;  % Dinkelbach inner-loop absolute gap threshold (eps_2, Alg. 1)

%% =====================================================================
%  Convenience Cell Array for SEE Evaluation
% ======================================================================
params.Consts = {params.P_lim, params.Qa, params.Qb, ...
                 params.af, params.beta0, params.N0, params.BW};

end
