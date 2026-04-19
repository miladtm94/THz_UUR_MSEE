function fc = flight_constants()
% FLIGHT_CONSTANTS  Return rotary-wing UAV propulsion model constants.
%
%   fc = flight_constants() returns a struct containing the mechanical
%   constants for the fixed-altitude rotary-wing UAV propulsion power model
%   defined in Eq. (11) of:
%
%   M. T. Mamaghani and Y. Hong, "Terahertz Meets Untrusted UAV-Relaying:
%   Minimum Secrecy Energy Efficiency Maximization via Trajectory and
%   Communication Co-design," IEEE Transactions on Communications, 2022.
%
%   The instantaneous propulsion power at time slot n is given by (Eq. 11):
%
%   P_f[n] = P0 * (1 + 3*||v[n]||^2 / (Omega^2 * R_u^2))   (blade profile)
%           + 0.5 * d0 * rho * s * A * ||v[n]||^3           (parasite)
%           + Pi * sqrt( sqrt(1 + ||v[n]||^4/(4*nu0^4))
%                        - ||v[n]||^2 / (2*nu0^2) )          (induced)
%
%   where the hovering blade-profile power is:
%       P0 = delta * rho * s * A * Omega^3 * R_u^3 / 8
%   and the hovering induced power is:
%       Pi = (1 + k_corr) * W^(3/2) / sqrt(2 * rho * A)
%
%   OUTPUT:
%       fc  - struct with fields:
%       Omega  : Blade angular velocity [rad/s]          (omega)
%       R_rotor: Rotor radius [m]                        (r)
%       rho    : Air density [kg/m^3]
%       s      : Rotor solidity (dimensionless)
%       A      : Rotor disk area [m^2]
%       nu0    : Mean rotor induced velocity in hover [m/s]
%       d0     : Fuselage drag ratio (dimensionless)
%       delta  : Profile drag coefficient (dimensionless)
%       k_corr : Induced power correction factor
%       W      : UAV gross weight [N]
%       P0     : Blade profile power in hover [W]
%       Pi     : Induced power in hover [W]
%       Vmax   : Maximum UAV speed [m/s]   (from [27] Table I)
%       Amax   : Maximum UAV acceleration [m/s^2]
%
% REFERENCE:
%   Y. Zeng, J. Xu, and R. Zhang, "Energy minimization for wireless
%   communication with rotary-wing UAV," IEEE Trans. Wireless Commun.,
%   vol. 18, no. 4, pp. 2329-2345, Apr. 2019. (Ref. [27] in the paper)

% Rotor mechanics
fc.Omega   = 300;    % Blade angular velocity [rad/s]
fc.R_rotor = 0.4;    % Rotor radius [m]
fc.rho     = 1.225;  % Air density [kg/m^3] (at sea level, 15°C)
fc.s       = 0.05;   % Rotor solidity (blade area / disk area)
fc.A       = 0.503;  % Rotor disk area [m^2]
fc.nu0     = 4.03;   % Mean rotor induced velocity in hover [m/s]
fc.d0      = 0.6;    % Fuselage drag ratio

% Power model coefficients
fc.delta   = 0.012;  % Profile drag coefficient
fc.k_corr  = 0.1;    % Induced power correction factor
fc.W       = 20;     % UAV gross weight [N] (approx. 2 kg UAV)

% Derived hover powers (from Eq. 11 definitions)
fc.P0 = fc.delta * fc.rho * fc.s * fc.A * fc.Omega^3 * fc.R_rotor^3 / 8;
% Nominal value: P0 ≈ 79.86 W
fc.Pi = (1 + fc.k_corr) * fc.W^(3/2) / sqrt(2 * fc.rho * fc.A);
% Nominal value: Pi ≈ 88.63 W

% Mobility limits (Constraint C12, C13 from Eq. 13)
fc.Vmax = 20;  % Maximum UAV speed [m/s]
fc.Amax = 5;   % Maximum UAV acceleration [m/s^2]

end
