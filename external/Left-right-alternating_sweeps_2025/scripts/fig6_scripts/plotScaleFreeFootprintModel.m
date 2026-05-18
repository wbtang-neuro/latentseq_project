function plotScaleFreeFootprintModel()

ngrid = 101;

% Create 2D coordinate grid and calculate distances and directions from
% each point to the origin
gv = linspace(-1, 1, ngrid);
[xx, yy] = meshgrid(gv);
dd = hypot(xx,yy);
aa = atan2(yy,xx);

% Calculate a Von-Mises function
kappa = 5;
zVonMises = circ_vmpdf(aa, pi/2, kappa);

% Calculate inverse-square distance
zInvSqDist = 1./(dd.^2);
 
zFootprint = zVonMises .* zInvSqDist;

zdata = {zVonMises, zInvSqDist, zFootprint};
names = ["Angular intensity function", "Distance intensity function", "Product"];

figure(colormap=bone());
for i = 1:3
    nexttile()
    z = zdata{i};
    z = sqrt(z); % some tone mapping to aid visualization
    imagesc(gv, gv, z);
    clim(prctile(z(:), [0, 99]));
    xline(0, 'w');
    yline(0, 'w');
    axis xy image off
    title(names(i));
    plotAgent(0, 0, height=0.6);
end

end