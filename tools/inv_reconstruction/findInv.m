%% Reconstructs inversion based on GB2001
% Assumes that between three mass levels, there resides an ambiguous layer containing the inversion layer at M + 1, with layer M being the layer immediately below the inversion.
%  ------- M + 3/2
% |       | <-- This point obtained from extrapolation of profile aloft
% |  M+1  | <-- Inversion somewhere in here
% |       | <-- This point = mean of thermodynamic values below
%  ------- M + 1/2

% In terms of the MYNN grid, this is:
%  ------- k + 1 <-- fluxes, w
% |       |
% |  (k)  | <-- theta, q
% |       |
%  ------- k <-- fluxes, w

% Outputs
% pi: inversion height pressure
% dThetavl: inversion jump in dThetavl
% pzoom: high resolution of pressure within ambiguous layer
% thetavlzoom: high resolution of thetavl within ambiguous layer

% Inputs
% tempVar: Should be a temperature variable that is ideally conserved
% within the STBL, i.e. theta_e or theta_l, theta_vl, etc.
% p: pressure
% p2: pressure between grid levels, at edges
% ki: index of layer directly below ambiguous layer

function [pi, mu, dTempVar, pzoom, tempZoom] = findInv(tempVar, p, p2, ki)
%% Set threshold

dThresh = 0.1; % 0.1 K is the threshold for temperature jump across inversion

%% GB2001 algorithm

% Initialize slope vector s
s = zeros(size(p) + 1); % Flux grid

s(ki+2) = (tempVar(ki+2) - tempVar(ki+1)) ./ (p(ki+2) - p(ki+1)); % Slope of theta directly above
s(ki+3) = (tempVar(ki+3) - tempVar(ki+2)) ./ (p(ki+3) - p(ki+2)); % Slope of theta one point above

splus = max(s(ki+2), s(ki+3)); % Choose slope closer to 0

tempVar2(ki+2) = tempVar(ki+2) +  splus .* (p2(ki+2) - p(ki+2)); % Interpolated temperature between grid points

% Solve quadratic equation -- the smaller root is always the physically meaningful solution with 0 < mu <= 1. Mu =1 corresponds to level k, i.e. the bottom boundary of the ambiguous layer.

% mu = (p2(ki+1) - pi) ./ (p2(ki+1) - p2(ki)); % pi = pressure at inversion
a = - (splus ./ 2) .* (p2(ki+2) - p2(ki+1));
b = (tempVar2(ki+2) - tempVar(ki));
c = - (tempVar(ki+1) - tempVar(ki));

mu = (-b + sqrt(b.^2 - 4*a*c)) ./ (2.*a);

dTempVar = tempVar2(ki+2) - splus .* mu .* (p2(ki+2) - p2(ki+1)) - tempVar(ki); % Inversion jump
pi = -(p2(ki+2) - p2(ki+1)).*mu + p2(ki+2); % Inversion pressure level

%% Zoom into inversion and plot

pzoom = p(ki):-1:p(ki+2);

for n = 1:length(pzoom)
	if pzoom(n) > pi
		tempZoom(n) = tempVar(ki);
	elseif pzoom(n) < pi
		tempZoom(n) = tempVar2(ki+2) + splus .* (pzoom(n) - p2(ki+2));
	end
end

end