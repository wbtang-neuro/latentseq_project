function [axes, angles] = idealGridAxes()
% Create idealized hexagonal grid axes
angles = [0, pi/3, 2*pi/3];
axes = zeros(3, 2);
for a = 1:3
    axes(a, :) = rotate2d([1, 0], angles(a));
end
end