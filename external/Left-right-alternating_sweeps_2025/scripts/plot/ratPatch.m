function hRat = ratPatch(p)
%RATPATCH Summary of this function goes here
%   Detailed explanation goes here
arguments
    p.color = [1,1,1];
    p.edgeColor = [0,0,0];
    p.sizeMeters = .035;
    p.orientation = pi/2;
    p.position = [0,0];
    p.hRat = [];
    p.ax = gca;
end
ydata = ([0,.5,1,.5,0]-.5).*p.sizeMeters;
xdata = ([0,.3,0,1,0]-.5).*p.sizeMeters;
[xdata, ydata] = rotate2d(xdata, ydata, p.orientation);
xdata = xdata + p.position(:, 1); ydata = ydata + p.position(:, 2);
if ~isempty(p.hRat)
    hRat = p.hRat;
    hRat.XData = xdata;
    hRat.YData = ydata;
else
    hRat = patch(p.ax, xdata, ydata, p.color, 'edgeColor', p.edgeColor);
end
end

