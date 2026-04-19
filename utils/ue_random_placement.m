function [x, y] = ue_random_placement(xc, yc, r_min, r_max, K)
% UE_RANDOM_PLACEMENT  Generate K UE positions uniformly in an annulus.
%
%   [x, y] = ue_random_placement(xc, yc, r_min, r_max, K) places K ground
%   UEs uniformly at random in an annular region centered at (xc, yc) with
%   inner radius r_min and outer radius r_max.
%
%   The angular positions are evenly spaced (2*pi/K apart) and the radial
%   distances are drawn to ensure area-uniform distribution via the square
%   root transform:
%       r = sqrt(r_min^2 + (r_max^2 - r_min^2) * U)
%   where U ~ Uniform(0,1), following the inverse CDF of the uniform
%   annulus radial distribution.
%
%   INPUTS:
%       xc    : x-coordinate of annulus center [m]
%       yc    : y-coordinate of annulus center [m]
%       r_min : Inner radius of annulus [m]  (= R1 in the paper, Section II)
%       r_max : Outer radius of annulus [m]  (= R2 in the paper)
%       K     : Number of UEs to place
%
%   OUTPUTS:
%       x  : [1xK] UE x-coordinates [m]
%       y  : [1xK] UE y-coordinates [m]
%
%   NOTE:
%       The caller should set rng default (or a specific seed) before
%       calling this function to ensure reproducibility (Section IV).
%
% See also: system_params.

% Area-uniform radial distances (inverse CDF method)
r = sqrt(r_min^2 + (r_max^2 - r_min^2) * rand(1, K));

% Uniformly-spaced angular positions [rad]
theta = (2 * pi / K) * (1:K);

% Cartesian coordinates
x = xc + r .* cos(theta);
y = yc + r .* sin(theta);

end
