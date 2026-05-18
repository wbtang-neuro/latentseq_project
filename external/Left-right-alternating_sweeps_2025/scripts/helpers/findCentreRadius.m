function cFieldRadius = findCentreRadius(aCorr, fieldThreshold)
% Adapted from BNT analysis.gridnessScore. The only changes are:
% 1) remove two acorr size arguments which supplied redundant information
% 2) Add a default value for argument "fieldThreshold"

% Added RJG
if nargin < 2 || isempty(fieldThreshold), fieldThreshold = 0.2; end

halfSize = ceil(size(aCorr)/2);
half_height = halfSize(1);
half_width = halfSize(2);
% /RJG

%     figure, contour(aCorr, [fieldThreshold fieldThreshold]), hold on;

% create a contour plot of fields
aCorr = double(aCorr);
cField = contourc(aCorr, [fieldThreshold fieldThreshold]);
[~, fLoc] = find(cField(1, :) == fieldThreshold);
if length(fLoc) > 1
    fLocMeans = zeros(length(fLoc), 2);
    allFields = cell(length(fLoc), 1);
    for i = 1:length(fLoc)-1 % the last contour is processed differently.
        allFields{i} = cField(:, fLoc(i)+1:fLoc(i+1)-1);
        fLocMeans(i, :) = [mean(allFields{i}(1, :)) mean(allFields{i}(2, :))];
        %             text(fLocMeans(i, 1), fLocMeans(i, 2), num2str(i));
    end
    i = i + 1;
    allFields{i} = cField(:, fLoc(i)+1:end);
    fLocMeans(i, :) = [mean(allFields{i}(1, :)) mean(allFields{i}(2, :))];

    % The 'min' approach is a bit faster (~ 0.02s for 1000 calculations), however sqDistance
    % appraoch is somewhat 'more correct'.
    % !!! Min approach gives an error on this data:
    % points = helpers.hexGrid([0 0 50 50], 15);
    % rmap = helpers.gauss2d(points, 10*ones(size(points, 1), 1), [50 50]);
    % aCorr = xcorr2(rmap - mean(rmap(:)));

    %         [~, fLocMin] = min(abs(mean([fLocMeans(:, 1) - half_width fLocMeans(:, 2) - half_height], 2)));
    %         [~, fLocMin] = min(sqDistance(fLocMeans', [half_width half_height]')); % point should be in format [x y]

    % get all distances and check two minimums of them
    allDistances = sqDistance(fLocMeans', [half_width half_height]'); % point should be in format [x y]
    [~, sortIndices] = sort(allDistances);
    twoMinIndices = sortIndices(1:2);

    if abs(allDistances(twoMinIndices(1)) - allDistances(twoMinIndices(2))) < 1
        % two fields with close middle points. Let's select one with minimum square
        areas = zeros(length(twoMinIndices), 1);
        for i = 1:length(twoMinIndices)
            testedField = allFields{twoMinIndices(i)};
            areas(i) = polyarea(testedField(1, :), testedField(2, :));

            % check that this polygon actually contains the middle point
            % if ~inpolygon(half_width, half_height, testedField(1, :), testedField(2, :))
            %     areas(i) = Inf;
            % end
        end
        % if all(isinf(areas))
        %     cFieldRadius = -1;
        %     return;
        % end
        [~, fLocMin] = min(areas);
        fLocMin = twoMinIndices(fLocMin);
    else
        fLocMin = twoMinIndices(1); % get the first minimum
    end

    centerField = allFields{fLocMin};
else
    centerField = cField(:, fLoc+1:end);
end
cFieldRadius = floor(sqrt(polyarea(centerField(1, :), centerField(2, :))/pi));
end

function D = sqDistance(X, Y)
    D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);
end