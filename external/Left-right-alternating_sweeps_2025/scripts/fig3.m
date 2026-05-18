%% Run CCG analysis
sweepsSetup
computeXcorrs(S.recs_of_mec(:), S)
%%
fld = fullfile(S.dataRoot, "results","xcorrs");
extractConnections(sourceFolder=fld)
%% Load cross-correlation analyses
sweepsSetup
fld = fullfile(S.dataRoot, "results","xcorrs");
p.basedir = fullfile(fld, "filtered");
%% Load
p.recNames = S.recs_of_mec(:);
p.loadHeavy = 1; % Set to 1 for examples
p.loadTuning = 1;
allrecs = loadConnRes(p);
%% Filter
p.min_spikes = 2000;
p.min_count = 1000; %
p.min_peakZ = 5;
p.min_peakExcess = 0;
p.only_good = 0;
p.only_functional = 1;
p.lump_dir = 0; 
p.lump_grid = 1;
recs = filterConnRes(allrecs(:), p); % Fewer detected conns bc only running periods are incl in CCG calc and not immobility/sleep periods (ref Schwindel et al 2014)

%% Plot example connections
clf
plotConnExamples(recs)

%% Plot tuning relationships
clf
tiledlayout(1,3);
[res, recs] = plotDirTuningConns(recs, "id_id", 1);
[res, recs] = plotDirTuningConns(recs, "id_conjunctive", 1);
[res, recs] = plotPhaseOffsetDirConns(recs, "conjunctive_grid", 1);

%% Schematic: Plot a conjunctive layer projecting onto a grid cell layer to generate a sweep
tl = tiledlayout(2,1);
boxcorners = [-1,1,1,-1,-1;-1,-1,1,1,-1]'.*2;
nexttile;
[corners, centers] = getGridTiles();
sigma = .01;
gv = linspace(-2,2,200);
Y = getGridFrame(0,0,sigma, gv);
% clf
alpha = reshapeSquare(Y);
img = [];
img(1,1,:)= S.col_conj;
img = repmat(img, [size(alpha), 1]);
himg = imshow(img, 'XData', gv, 'YData', gv);
himg.AlphaData = alpha./max(alpha(:));
plot(corners(:, 1), corners(:, 2), 'color', [1,1,1]*.5, 'LineWidth',.1)
ax = gca; ax.XLim = minmax(gv); ax.YLim = minmax(gv);
view(10, 15)
ax.ZLim = [-2,0];
plot3([0,0], [0,0], [0,-2], 'k');
plot3([0,.5], [0,0], [0,-2], 'color', S.col_conj);

plot(boxcorners(:, 1), boxcorners(:, 2), 'color', [1,1,1]*.1, 'LineWidth',1)
camproj('perspective')
% Make a sweep
nexttile;
frame = alpha;
nsteps = 10;
for s = 1:nsteps
    color =[];
    color(1, 1, :) = S.col_grid+[.5,1,.5]*nsteps*.1-[1,1,.4]*s*.1;
    color = repmat(color, [size(frame), 1]);
    hdec(s) = imshow(color, 'XData',gv, 'YData',gv);
    frame = getGridFrame(.05*s,0,sigma, gv);
    frame = reshapeSquare(frame);
    hdec(s).AlphaData = gather(frame./max(frame(:)));
end
plot(corners(:, 1), corners(:, 2), 'color', [1,1,1]*.5, 'LineWidth',.1)
plot(boxcorners(:, 1), boxcorners(:, 2), 'color', [1,1,1]*.1, 'LineWidth',1)
ax = gca; ax.XLim = minmax(gv); ax.YLim = minmax(gv);
view(10, 15)
ax.ZLim = [-2,0];
camproj('perspective')

%% Schematic: Layer of ID cells projecting to conjunctive cells
figure
clf;
tl = tiledlayout(1,2);
nexttile
nlayers = 6;
pts = 1:nlayers;
cols = mapcolors(pts', [0,nlayers], 'hsv');
cols = S.col_id + [0,.5,0].*linspace(+.2, .8, nlayers)';
scatter(pts'*0, pts',1550,cols, 'MarkerEdgeColor','k', 'MarkerFaceAlpha',.5)
plot(vec(pts*0+[0.5;2;nan]), vec(pts-.02+[0;0;nan]),'Color','k', 'LineWidth',2.2)
x = pts*0+[0;0;nan];
y = pts+[0;0;nan];
angles = linspace(0, 2*pi, nlayers+1);
angles = fliplr(angles)
[dx, dy] = pol2cart(angles(2:end), .5);
x(2, :) = .5*dx; x(1, :) = -.5*dx; 
y(2, :)= y(2, :)+.5*dy; y(1, :)= y(1, :)-.5*dy;
plot(vec(x)-0, vec(y),'Color','k', 'LineWidth',4)

for i = 1:nlayers
    ratPatch(color=[0,0,0], position=[x(2, i), y(2,i)], orientation=angles(i+1), sizeMeters=.3);
end
ax = gca;
ax.YLim = [0, 2*nlayers];
axis off equal

nexttile
%
%# create stacked images (I am simply repeating the same image 5 times)
fig = gcf;
fig.Renderer = "opengl";
[corners, centers] = getGridTiles();
sigma = .01;
sigma = .03;
gv = linspace(-1.7,1.7,200);
Y = getGridFrame(0,0,sigma, gv);
% clf
alpha = reshapeSquare(Y);
alpha = alpha./max(alpha(:));
clear img;
img = repmat(alpha, 1, 1, 3);
img(:, :, [2,3])=0;
img(:, :, [1])=1;
whiteimg = img;
whiteimg(:) = 1;
% img.X = 255*alpha./max(alpha(:));
% I = repmat(img,[1 1 5]);
% cmap = img.map;

%# coordinates
Z = ones(size(img,1),size(img,2));

%# plot each slice as a texture-mapped surface (stacked along the Z-dimension)
boxcorners = [-1,1,1,-1,-1;-1,-1,1,1,-1]'.*1.7;

% shadow
gvShade = linspace(-1.7,1,200);
alphaShade = alpha*0;
alphaShade(end-180:end, 1:180)=.2;
alphaShade = imgaussfilt(alphaShade, 10);
colShade = whiteimg*0;

% nlayers = 5;
for k=1:nlayers
    surface('XData',gv, 'YData',gv, 'ZData',Z.*k*1, ...
        'CData',whiteimg, 'CDataMapping','direct', ...
        'EdgeColor','none')
    color =[];
    color(1,1,:) = mapcolors(k, [0,nlayers], 'hsv');
    color(1,1,:) = S.col_conj+[1,1,.4]*nlayers*.1-[1,1,.4]*k*.1;
    color = repmat(color, [size(alpha), 1]);
    surface('XData',gv, 'YData',gv, 'ZData',Z.*k*1, ...
        'CData',color, 'CDataMapping','direct', ...
        'EdgeColor','none', 'FaceAlpha', 'flat','AlphaData', alpha)

    % shade
    if k<nlayers
        surface('XData',gv, 'YData',gv, 'ZData',Z.*k*1+.03, ...
        'CData',colShade, 'CDataMapping','direct', ...
        'EdgeColor','none', 'FaceAlpha', 'flat','AlphaData', alphaShade)
    end
    plot3(corners(:, 1), corners(:, 2),1*k+0*corners(:, 2), 'color', [1,1,1]*.7, 'LineWidth',.1)
    plot3(boxcorners(:, 1), boxcorners(:, 2), 1*k+0*boxcorners(:, 2), 'color', [1,1,1]*.1, 'LineWidth',1)
end


ax = gca; ax.XLim = minmax(gv); ax.YLim = minmax(gv); ax.ZLim = [0, 2*k];
camproj('perspective')
view(10, 15)
axis off
%% 
function frame = getGridFrame(x0, y0, sigma, gv)
    [corners, centers] = getGridTiles();
    [xi, yi] = meshgrid(gv);
    X = [xi(:), yi(:)];
    Y = zeros(length(X), 1);
    for i = 1:length(centers)
        Y = Y+mvnpdf(X, centers(i, :)+[x0, y0], [1,1]*sigma);
    end
    frame = Y;
end