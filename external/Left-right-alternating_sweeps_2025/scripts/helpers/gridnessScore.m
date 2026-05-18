% Calculate gridness score for an autocorrelogram
%
% Calculates a gridness score by expanding a circle around the centre field and
% calculating a correlation value of the expanded circle with it's rotated versions.
% The expansion is done up until the smallest side of the autocorrelogram.
% Can also calculate grid statistics.
%
%  USAGE
%   [score, <stats>] = analyses.gridnessScore(aCorr, <options>)
%   aCorr       A 2D autocorrelogram.
%   <options>   Optional list of property-value pairs (see table below)
%
%   ==============================================================================================
%    Properties    Values
%   ----------------------------------------------------------------------------------------------
%    'threshold'        Normalized threshold value used to search for peaks on the autocorrelogram.
%                       Ranges from 0 to 1, default value is 0.2.
%    'minOrientation'   Value of minimal difference of inner fields orientation (in degrees). If
%                       there are fields that differ in orientation for less than minOrientation,
%                       then only the closest to the centre field are left. Default value is 15.
%   ==============================================================================================
%   score       Gridness score. Ranges from -2 to 2. 2 is more a theoretical bound for a perfect grid.
%               More practical value is around 1.3.
%   stats       If this variable is requested, then it is a structure with the following statistics:
%       spacing         3-element vector with distances from the centre field to neighbour fields.
%       orientation     3-element vector with orientations between the centre field and neighbour fields.
%       ellipse         Ellipse fitted to the grid. Contains the centre, radii and orientation in
%                       radians, stored as [Cx, Cy, Rx, Ry, theta].
%       ellipseTheta    Radius of the ellipse in degrees wrapped in range [0..180].
%
function [gscore, varargout] = gridnessScore(aCorr, varargin)
    nout = max(nargout, 1) - 1;
    inp = inputParser;
    defaultFieldThreshold = 0.2;
    defaultMinOrientation = 15;
    gridStat.orientation = nan(3, 1);
    gridStat.spacing = nan(3, 1);
    gridStat.ellipse = nan(1, 5);
    gridStat.ellipseTheta = nan;

    % input argument check functions
    checkThreshold = @(x) helpers.isdscalar(x, '>0', '<1');
    checkDScalar = @(x) helpers.isdscalar(x, '>0');

    % fill input parser object
    addRequired(inp, 'aCorr');
    addParameter(inp, 'threshold', defaultFieldThreshold, checkThreshold);
    addParameter(inp, 'minOrientation', defaultMinOrientation, checkDScalar);
    addParameter(inp, 'minCenterFieldRadius', 2);

    parse(inp, aCorr, varargin{:});

    % get parsed arguments
    fieldThreshold = inp.Results.threshold;
    minOrientation = inp.Results.minOrientation;
    minCenterFieldRadius = inp.Results.minCenterFieldRadius;

    halfSize = ceil(size(aCorr)/2);
    half_height = halfSize(1);
    half_width = halfSize(2);
    aCorrRad = min(halfSize);
    aCorrSize = size(aCorr);
    
    if aCorrSize(1) == 1 || aCorrSize(2) == 1
        gscore = nan;
        if nout > 0
            varargout{1} = gridStat;
            varargout{2} = 0;
            varargout{3} = nan(6, 2);
            varargout{4} = 0;
            varargout{5} = [];
        end
        return;
    end

    % contourc is efficient if aCorr is normalized
%     maxValue = max(max(aCorr));
%     if maxValue ~= 1
%         aCorr = aCorr / maxValue;
%     end
    aCorr = aCorr / max(aCorr, [], "all");
    cFieldRadius = findCentreRadius(aCorr, fieldThreshold);
%     aCorr = single(aCorr);
    % if cFieldRadius == -1
    %     % let's try with increased threshold, this might help
    %     cFieldRadius = findCentreRadius(aCorr, fieldThreshold + 0.1, half_width, half_height);
    % end
    if cFieldRadius<minCenterFieldRadius || cFieldRadius>=aCorrRad % RJG 2022-11-09, making it work for coarse-binned acorrs
        gscore = nan;
        if nout > 0
            varargout{1} = gridStat;
            varargout{2} = 0;
            varargout{3} = nan(6, 2);
            varargout{4} = 0;
        end
        return;
    end

    % Meshgrid for expanding circle
    [rr, cc] = meshgrid(1:size(aCorr, 2), 1:size(aCorr, 1));

    % Define iteration radius step size for the gridness score
    radSteps = cFieldRadius:aCorrRad;
    radSteps(1) = [];
    numSteps = length(radSteps);

    GNS = zeros(numSteps, 2);
    GNS(:, 2) = radSteps;
    
    nrot = 5;
    rotAngles_deg = 30*(1:nrot);

    distFromCenter = sqrt((cc - half_height).^2 + (rr - half_width).^2);
    cFieldExclusionMask = distFromCenter > cFieldRadius;
    distFromCenter1 = distFromCenter(:);
    cFieldExclusionMask1 = cFieldExclusionMask(:);

    sz = [size(aCorr), nrot];
    rotACorr = zeros(sz);
%     aCorr = single(aCorr);
%     rotACorr = zeros(sz, "single");
    for i = 1:nrot    
        rotACorr(:, :, i) = images.internal.builtins.imrotate(aCorr, rotAngles_deg(i), size(aCorr), 'bilinear');
    end
    rotACorr = reshape(rotACorr, numel(aCorr), nrot);
    rotACorr = rotACorr(cFieldExclusionMask1, :);
    aCorr0 = aCorr;
    aCorr = aCorr(cFieldExclusionMask);
    distFromCenter1 = distFromCenter1(cFieldExclusionMask1);

    % Define expanding ring of autocorrellogram and do x30 correlations
    for i = 1:numSteps
        ind = distFromCenter1 < radSteps(i);
        %         ind = (cFieldExclusionMask1 & (distFromCenter1 < radSteps(i)));
        refCircle = aCorr(ind);
        rotCircle = rotACorr(ind, :);
        rotCorr = corrPearson(refCircle, rotCircle);
        GNS(i, 1) = min(rotCorr([2, 4])) - max(rotCorr([1, 3, 5]));
    end
    aCorr = aCorr0;

    % Find the biggest gridness score and radius
    numGridnessRadii = 3;
    numStep = numSteps - numGridnessRadii;
    if numStep < 1
        numStep = 1;
    end

    i0 = 1:numGridnessRadii;
    ii = (1:numStep) + i0' - 1;
    meanGridnessArray = mean(reshape(GNS(ii, 1), size(ii)), "omitnan");
    [gscore, gInd] = max(meanGridnessArray);
    gscoreLoc = gInd + (numGridnessRadii-1)/2;

    varargout{4} = radSteps(gscoreLoc);

    % Return if we do not need to calculate grid statistics
    if nout < 1
        return;
    end

    % Calculate gridness score statistics
    grad = radSteps(gscoreLoc);
    w = grad / 4;
    maskOuter = distFromCenter < grad+w;
%     maskInner = distFromCenter > grad-w; % RG 2023-01-16: add inner constraint to prevent weak maxima "winning"
    aCorrSm = imgaussfilt(aCorr, grad/(2*pi)/2, "padding", "replicate"); % RG 2023-01-16 eliminate spurious maxima
    varargout{5} = aCorrSm;

    bestCorr = maskOuter .* aCorrSm;

    regionalMaxMap = imregionalmax(bestCorr, 4);
    se = strel('square', 3);
    im2 = imdilate(regionalMaxMap, se); % dilate map to eliminate fragmentation
    cc = bwconncomp(im2, 8);
    stats = regionprops(cc, 'Centroid');

    if length(stats) < 5
        warning('BNT:numFields', 'Not enough inner fields has been found. Can''t calculate grid properties');

        varargout{1} = gridStat;
        varargout{2} = cFieldRadius;
        varargout{3} = nan(6, 2);
        varargout{4} = radSteps(gscoreLoc);
        return;
    end

    allCoords = [stats(:).Centroid];
    centresOfMass(:, 1) = allCoords(1:2:end);
    centresOfMass(:, 2) = allCoords(2:2:end);

    % Calculate orientation for each field relative to the centre field
    orientation = (atan2(centresOfMass(:, 2) - half_height, centresOfMass(:, 1) - half_width)); % atan2(Y, X)
    peaksToCentre = sqDistance(centresOfMass', [half_width half_height]');
    zeroInd = find(orientation == 0, 1);
    orientation(zeroInd) = []; % remove zero value, so that we do not have a side effect with minOrientation
    stats(zeroInd) = [];
    peaksToCentre(zeroInd) = [];
    centresOfMass(zeroInd, :) = [];

    % filter fields that have similar orientation
    orientDistSq = circ_dist2(orientation);
    closeFields = abs(orientDistSq) < deg2rad(minOrientation);
    [rows, cols] = size(closeFields);
    closeFields(1:(rows+1):rows*cols) = 0; % assign zero to diagonal elements
    closeFields(tril(true(rows))) = 0; % assign zero to lower triangular of a matrix. Matrix is
                                       % symmetric and we do not need these values.
    [rows, cols] = find(closeFields); % find non-empty elements, they correspond to indices of close fields
    if ~isempty(rows)
        indToDelete = zeros(1, length(rows));
        for i = 1:length(rows)
            % fieldPeaks = [fields([rows(i) cols(i)]).peakX; fields([rows(i) cols(i)]).peakY];
            % peaksToCentre = sqDistance(fieldPeaks, [half_width; half_height]);
            if peaksToCentre(rows(i)) > peaksToCentre(cols(i))
                indToDelete(i) = rows(i);
            else
                indToDelete(i) = cols(i);
            end
        end
        indToDelete = unique(indToDelete);
        stats(indToDelete) = [];
        peaksToCentre(indToDelete) = [];

        if length(stats) < 4
            warning('BNT:numFields', 'Not enough inner fields has been found. Can''t calculate grid properties');

            varargout{1} = gridStat;
            varargout{2} = cFieldRadius;
            varargout{3} = nan(6, 2);
            varargout{4} = radSteps(gscoreLoc);
            return;
        end

        allCoords = [stats(:).Centroid];
        clear centresOfMass;
        centresOfMass(:, 1) = allCoords(1:2:end);
        centresOfMass(:, 2) = allCoords(2:2:end);
    end

    % % get fields peak coordinates
    % fieldPeaks = zeros(length(stats), 2);
    % for i = 1:length(stats)
    %     [~, maxInd] = max(bestCorr(stats(i).PixelIdxList));
    %     fieldPeaks(i, :) = stats(i).PixelList(maxInd, :);
    % end
    % % fieldPeaks = [fields(:).peakX; fields(:).peakY]'; %
%     peaksToCentre = sqDistance(centresOfMass', [half_width half_height]');
    [~, sortInd] = sort(peaksToCentre);
    stats = stats(sortInd);
    centresOfMass = centresOfMass(sortInd, :);

    % leave only 6 closest neighbours (if available)
    if length(stats) > 5
%         stats = stats(1:6);
        centresOfMass = centresOfMass(1:6, :);
    else
%         stats = stats(1:end);
        centresOfMass = centresOfMass(1:end, :);
    end

    % centresOfMass = [fields(:).x; fields(:).y]';
    % Calculate orientation for each field relative to the centre field
    orientation = rad2deg(atan2(centresOfMass(:, 2) - half_height, centresOfMass(:, 1) - half_width)); % atan2(Y, X)

    % Calculate distances between centre of masses for each field and the centre field
    spacing = sqrt((centresOfMass(:, 1) - half_width).^2 + (centresOfMass(:, 2) - half_height).^2);

%     % Plot grid polygon points
%     figure, plot.colorMap(bestCorr), hold on;
%     plot(centresOfMass(:, 1), centresOfMass(:, 2), '+k', 'markersize', 8);

    ell = fitEllipse(centresOfMass(:, 1), centresOfMass(:, 2));
    ellipseTheta = rad2deg(wrapToPi(ell(end)) + pi);
%     drawEllipse(ell, 'linewidth', 2, 'color', [1 1 1]);

    % Determine axes orientation, spacing and deviation
    [~, bBC] = sort(abs(orientation));
    [~, bBC2] = sort(abs(orientation - orientation(bBC(1))));

    % leave only three values, because autocorrelogram is symmetric
    orientation = orientation(bBC2(1:3));
    spacing = spacing(bBC2(1:3));
    [orientation, orientSortInd] = sort(orientation);

    spacing = spacing(orientSortInd);

    gridStat.orientation = orientation;
    gridStat.spacing = spacing;
    gridStat.ellipse = ell;
    gridStat.ellipseTheta = ellipseTheta;

    varargout{1} = gridStat;
    varargout{2} = cFieldRadius;
    varargout{3} = centresOfMass;
    varargout{4} = radSteps(gscoreLoc);
end

function D = sqDistance(X, Y)
    D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);
end