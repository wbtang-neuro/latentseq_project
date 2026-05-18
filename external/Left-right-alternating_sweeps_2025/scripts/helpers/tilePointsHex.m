function [xt, yt, tileShifts] = tilePointsHex(x, y, gridAxes)
% Generate all seven possible tilings of a set of points for hexagonal 
% grid axes

x = x(:);
y = y(:);
n = numel(x);

% p = [x, y]; 

xt = zeros(n, 7);
yt = zeros(n, 7);

xt(:, 1) = x;
yt(:, 1) = y;
c = 1;

tileShifts = zeros(7, 2);

for a = 1:3 % iterate through axes
    for pol = [-1, 1] % both polarities
        c = c+1;
        gax = pol*gridAxes(a, :);
        xt(:, c) = x + gax(1);
        yt(:, c) = y + gax(2);
        tileShifts(c, :) = gax;
    end
end

end