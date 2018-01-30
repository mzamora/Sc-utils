% Generates new figure at location defined by id because I'm too lazy to move things around
% Author: Handa Yang

function out = newFig(id, w, h, fs)

if nargin < 1, id = 1; end
if (nargin < 2 || nargin < 3), w = 1680; h = 1050; end % Default resolution is 1680 x 1050 (EBU2 304 desktop)
if nargin < 4, fs = 16;

% Arrange figures in 2x3 grid (i.e. 2 rows of 3), starting from left monitor to right monitor
for i = 1:3
	location{i} = [ (10 + ((i-1)/3)*w), (h/2), ( (w-60)/3 ), ( (h-200)/2 ) ];
	embed{i} = [ 0.05 + (i-1)*.325, 0.58, 0.28, 0.4 ];
end

for i = 4:6
	location{i} = location{i-3} - [0, (h/2 - 5), 0, 0];
	embed{i} = embed{i-3} - [0, 0.5, 0, 0];
end

for i = 7:12
	location{i} = location{i-6} + [w, 0, 0, 0];
end

if id < 13
	out = figure('Position', location{id}, 'Color', [1 1 1]);
else
	out(1) = figure('Position', [40 25 1600 900], 'Color', [1 1 1]);
	for j = 1:6
		out(j+1) = axes('Position', embed{j}, 'fontsize', fs);
	end
end

end