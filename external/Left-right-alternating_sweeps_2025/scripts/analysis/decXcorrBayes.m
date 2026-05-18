function [peakPos, peakProb, torusPhase, info, Xcorrs, P] = decXcorrBayes(rateMaps, spikeCounts, sigma, gvTile, peakDetectionMethod)
% Decode single-module grid activity using LMT
%
% INPUTS
% mPos  : LMT position model
% Us    : array of Unit structs
% t     : bin times for decoding (must be contiguous)
% sigma : standard deviation of gaussian smoothing window
% gvTile: vector of gridded points on phase tile
%
% OUTPUTS
% peakPos   : position of peak in decoded probability distribution [time, xy]
% peakProb  : probability values at the positions in "peakPos"
% torusPhase: circular torus phase values of "peakPos"
% info      : various other useful stuff

if nargin < 6 || isempty(peakDetectionMethod), peakDetectionMethod = "regionalMax"; end

if isempty(gvTile)
    gvTile = (-1.5:0.01:1.5)';
end
%% CHUNK_SIZE = 1e4;
CHUNK_SIZE = 5e3;
% CHUNK_SIZE = 1e4;
nt = size(spikeCounts, 1);

S = SweepSettings;
% Loate the tuning of the supplied units in the LMT F matrix.
% cases where we *
F = rateMaps;
Fsize = sqrt(size(F, 1))*[1, 1];

acorr = calcAcorr(F);
[xxTile, yyTile] = meshgrid(gvTile);
tileGridSize = size(xxTile);
tilePoints = cat(3, xxTile, yyTile);

% Normalize
spikeCounts = gsmooth(spikeCounts, sigma);
Y = gather(spikeCounts);
tuning = gather(F./mean(F+eps, 'omitnan'));
[tform, interpPointsTf, idealAxes, idealAngles, pointsInTile, acorrRotation, gridAxes, spacing] = fitAcorrToIdealGrid(acorr, tilePoints);

c = 1;
nchunks = 0;

peakPos = zeros(nt, 2);
peakProb = zeros(nt, 1);
torusPhase = zeros(nt, 2);

if nargout > 4
   Xcorrs = zeros([tileGridSize, nt], "single");
end

%% Decode
tic; disp('Decoding')
if sigma<8
    probChunk = simpleBayesian(Y, tuning, S.dt, 0.5);
else
    probChunk = simpleBayesian(Y, tuning, S.dt, 15);
end
toc
P = probChunk;

% Do xcorr stuff
while c<nt
    nchunks = nchunks + 1;
    inds = c + (1:CHUNK_SIZE) - 1;
    if inds(end) > nt
        inds = c:nt;
    end
    fprintf("Chunk %u, bins [%u - %u], total %u\n", nchunks, inds([1, end]), nt);
    Pchunk = gpuArray(P(inds, :));

    [ppos, peakProb(inds), xc] = decodeOneChunk(Pchunk, Fsize, acorr, gvTile,...
        interpPointsTf, pointsInTile, peakDetectionMethod);

    [peakPos(inds, :)] = wrapPosToGridTile(ppos, idealAxes);
    if nargout > 4
        disp('\n')
        Xcorrs(:, :, inds) = xc;
    end

    c = c + numel(inds);
end

info = struct( ...
    "tform", tform, ...
    "gvTile", gvTile, ...
    "interpPointsTf", interpPointsTf, ...
    "sigma", sigma, ...
    "acorr", acorr, ...
    "acorrRotation", acorrRotation, ...
    "idealAxes", idealAxes, ...
    "gridAxes", gridAxes, ...
    "gridSpacing", spacing, ...
    "idealAngles", idealAngles, ...
    "pointsInTile", pointsInTile);

end

function [nearestPeakPos, peakProbability, xcorrs] = decodeOneChunk( ...
    P, fsz, acorr, gvTile, interpPointsTf, ...
    pointsInTile, peakDetectionMethod)

nt = size(P, 1);

% Caculate xcorr of the decoded probability distribution vs. the
% grid autocorrelogram, at every time step.
P3 = reshape(P', [fsz, nt]);
xcorrs3 = normxcorr2_general_stack_helper(acorr, P3);
peakProbability = max(P, [], 2);

% Trim the xcorrs to the same size as the acorr
gxc = centeredGridVector(size(xcorrs3, 1))';
vxc = abs(gxc)<=((size(acorr,1)/2)-1);
xcorrs3 = xcorrs3(vxc, vxc, :);
xcorrs3(isnan(xcorrs3)) = 0;
gxc = gxc(vxc);


% % Transform the xcorrs and acorr onto the standard phase tile, by
% % re-interpolating the xcorr at the coordinates of the rotated grid of
% % 2-D tile points
%
% Cubic interpolation works best, because the output is smooth and the
% peaks may lie anywhere on the interpolation grid.
xcorrsTf3 = gpuArray.zeros([size(pointsInTile), nt], "single");
for i=1:nt
    xcorrsTf3(:, :, i) = interp2(gxc, gxc, xcorrs3(:, :, i), interpPointsTf(:, :, 1), interpPointsTf(:, :, 2), 'cubic');
end
xcorrsTf2 = reshape(xcorrsTf3, [], nt);

if nargout==3
    xcorrs = gather(xcorrsTf2);
    xcorrs = reshape(xcorrs, [size(pointsInTile), nt]);
end

if strcmpi(peakDetectionMethod, "regionalMax")
    % Detects all peaks in each xcorr, then find the peak nearest to the center.
    % xcorrsTf3 = gather(xcorrsTf3);
    nans = isnan(xcorrsTf3);
    if any(nans, "all")
        warning("NaNs in interpolated xcorr! Probably trying to extrapolate beyond xcorr limits.");
        xcorrsTf3(nans) = 0;
    end
    ismax = imregionalmaxStack(xcorrsTf3, 8);
    ismax = gather(ismax);
    % ismax = imregionalmax(xcorrsTf3, 8);
    inds = find(ismax);
    [i,j,k] = ind2sub(size(xcorrsTf3), inds);
    mxx = gvTile(j);
    myy = gvTile(i);
    mxr = hypot(mxx, myy);
    C = extractGroupedDataMex(k, mxr, mxx, myy);
    nearestPeakPos = nan(nt, 2);
    for n=1:min(nt, max(k))
        Cn = C(n, :);
        rad = Cn{1};
        [~, imn] = min(rad);
        x = Cn{2}(imn);
        y = Cn{3}(imn);
        if isempty(x)
            x = nan;
            y = nan;
        end
        nearestPeakPos(n, :) = [x, y];
    end
elseif strcmpi(peakDetectionMethod, "max")
    % Find the maximum value within the central tile region. Much faster
    % than the above method, but perhaps slightly less accurate (not sure,
    % haven't tested).
    xcorrsTf2(~pointsInTile(:), :) = nan;
    [~, imx] = max(xcorrsTf2);
    [i, j] = ind2sub(size(pointsInTile), imx(:));
    nearestPeakPos = gvTile([j, i]);
end


end

function [axes, angles] = idealGridAxes()
% Create idealized hexagonal grid axes
angles = [0, pi/3, 2*pi/3];
axes = zeros(3, 2);
for a = 1:3
    axes(a, :) = rotate2d([1, 0], angles(a));
end
end

function [acorr, gv] = calcAcorr(F)
clear acorr
for u = 1:size(F, 2)
    z = reshapeSquare(F(:, u));
    acorr(:, :, u) = autocorrelation(z);
end
acorr = median(acorr, 3);
acorr(isnan(acorr)) = 0;
acorr = gather(single(acorr));
acorrSz = size(acorr);
gv = centeredGridVector(acorrSz(1))';

end

function [tform, pointsTf, idealAxes, idealAngles, pointsInTile, acorrRotation, gridAxes, spacing] = fitAcorrToIdealGrid(acorr, interpPoints)
% Find the grid axes of the autocorrelogram
[~, gstats] = gridnessScore(acorr);
[x, y] = pol2cart(deg2rad(gstats.orientation), gstats.spacing);
gridAxes = [x, y];
gridAxes(3, :) = gridAxes(2, :) - gridAxes(1, :);
spacing = mean(gstats.spacing, 'omitnan');
% Fit a "standardization" transformation to map the grid axes
% onto a perfect template grid
[idealAxes, idealAngles] = idealGridAxes();
idealAxes = gridAxes./spacing;
tform = fitgeotform2d(gridAxes, idealAxes, 'affine');

interpGridSize = size(interpPoints, [1, 2]);
interpPoints = reshape(interpPoints, [], 2);
pointsTf = tform.transformPointsInverse(interpPoints);
pointsTf = reshape(pointsTf, [interpGridSize, 2]);

% Calculate which points in the tile grid are within the central tile
[tx, ty] = gridHexTileCoords(idealAxes);
pointsInTile = inpolygon(interpPoints(:, 1), interpPoints(:, 2), tx, ty);
pointsInTile = reshape(pointsInTile, interpGridSize);

% Calulate the rotation of the transformation
T = tform.T;
rot1 = atan2(T(2, 1), T(1, 1)); % from rot. matrix first column: [sin(alpha); cos(alpha)]
rot2 = atan2(-T(1, 2), T(2, 2)); % from second column: [-sin(alpha);
fprintf("Rotation of two grid axes: %.0f°, %.0f° ...", rad2deg(rot1), rad2deg(rot2));
acorrRotation = circ_mean([rot1; rot2]);

end



