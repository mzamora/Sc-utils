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
% dX: inversion jump in X
% pzoom: high resolution of pressure within ambiguous layer
% Xzoom: high resolution of X within ambiguous layer

% Inputs
% X: Variable for which you wish to determine inversion jump
% p: pressure
% p2: pressure between grid levels, at edges
% ki: index of layer directly below ambiguous layer
% mu: above-inversion mass fraction--can be interpreted as inversion level
%     within grid point
% pi: inversion pressure level

function [dX, pzoom, Xzoom] = findJump(X, p, p2, ki, mu, pi)
%% GB2001 algorithm

s32 = (X(ki+2) - X(ki+1)) ./ (p(ki+2) - p(ki+1)); % Slope of theta directly above
s52 = (X(ki+3) - X(ki+2)) ./ (p(ki+3) - p(ki+2)); % Slope of theta one point above

splus = max(s32, s52); % Choose slope closer to 0
% splus = s32;
X32 = X(ki+2) + splus .* (p2(ki+2) - p(ki+2));

dX = 2 .* ( (X(ki+1) - (1 - mu).*X(ki)) ./ mu ) - X(ki) - X32;

%%
pzoom = p(ki):-1:p(ki+2);

Xzoom = zeros(size(pzoom)); % Init

for n = 1:length(pzoom)
	if pzoom(n) > pi
		Xzoom(n) = X(ki);
	elseif pzoom(n) < pi
		Xzoom(n) = X32 + splus .* (pzoom(n) - p2(ki+2));
	end
end

end