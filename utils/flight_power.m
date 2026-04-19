function P_flight = flight_power(V_uav)
% FLIGHT_POWER  Compute the instantaneous rotary-wing UAV propulsion power.
%
%   P_flight = flight_power(V_uav) evaluates the mechanical propulsion
%   power model from Eq. (11) for all time slots simultaneously.
%
%   INPUT:
%       V_uav    : UAV velocity matrix [2xN], where each column is the
%                  horizontal velocity vector [vx; vy] at a time slot.
%
%   OUTPUT:
%       P_flight : Instantaneous propulsion power [1xN] in Watts.
%
%   MODEL (Eq. 11):
%   P_f[n] = P0*(1 + 3*||v||^2/(Omega^2*R^2))   <- blade profile
%           + 0.5*d0*rho*s*A*||v||^3             <- parasite
%           + Pi*sqrt(sqrt(1+||v||^4/(4*nu0^4)) - ||v||^2/(2*nu0^2)) <- induced
%
%   The three terms correspond to blade profile, parasite (fuselage drag),
%   and induced (lift-related) power components, as detailed in:
%   Y. Zeng et al., IEEE Trans. Wireless Commun., 2019 (Ref. [27]).
%
%   ASSUMPTIONS:
%   - Fixed-altitude flight (no vertical velocity component).
%   - Quasi-static aerodynamic model valid for moderate speeds.
%
% See also: flight_constants, feasible_init.

fc = flight_constants();

speed_sq = norms(V_uav).^2;   % ||v[n]||^2 for all n, [1xN]
speed    = norms(V_uav);       % ||v[n]||   for all n, [1xN]

% Induced power auxiliary variable: mu[n] = sqrt(sqrt(1+v^4/4/nu0^4) - v^2/2/nu0^2)
% See Eq. (39) and Eq. (41) in the paper
mu_induced = sqrt( sqrt(1 + speed_sq.^2 / (4 * fc.nu0^4)) ...
                   - speed_sq / (2 * fc.nu0^2) );

% Blade profile power component
P_blade = fc.P0 * (1 + 3 * speed_sq / (fc.Omega^2 * fc.R_rotor^2));

% Parasite (fuselage drag) power component
P_parasite = 0.5 * fc.d0 * fc.rho * fc.s * fc.A * speed.^3;

% Induced (lift) power component
P_induced = fc.Pi * mu_induced;

P_flight = P_blade + P_parasite + P_induced;

end
