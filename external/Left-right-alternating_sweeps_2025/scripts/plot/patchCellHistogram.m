function hpatch=patchCellHistogram(p)
%PATCHCELLHISTOGRAM Summary of this function goes here
%   Detailed explanation goes here
arguments
    p.units
    p.idx
    p.edges
    p.xpos
    p.color
end
S = SweepsSettings;
centers = edg2cen(p.edges)';
bw = abs(centers(2)-centers(1));

occ = histcounts([p.units.shankPos], p.edges); 
h = histcounts([p.units(p.idx).shankPos], p.edges);
h = 100*h./occ;
h = 2*h;
h(isnan(h))=0;
x = ([h', h']);
x = vec(x');
y = ([centers, centers]+.5*bw*[-1,1]);
y = vec(y');
x = [0;x;0;0];
x = x +p.xpos;
y = [y(1);y;y(end);y(1)];

hpatch = patch(x, y, p.color, 'EdgeColor', 'none', 'FaceAlpha', .5);

end

