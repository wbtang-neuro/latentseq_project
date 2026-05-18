function [corners, centers] = getGridTiles()
%PLOTGRIDTILES Plot seven grid tiles defined by a set of gridaxes
% GridAxes: 3x2 vector specifying grid axes
if nargin==2
    col = [1,1,1]*.5;
end
%% get the coordinates corresponding to the corners of the central hexagon
% gridAxes = permute(gridAxes, [1,3,2]);
gridAxes = fliplr(idealGridAxes);
gridAxes = gridAxes*1.5;
originalAxes = gridAxes;
gridAxes = gridAxes./(2*cosd(30));
corners = [gridAxes; -gridAxes]; % flip axes to get other side of hex
corners(end+1, :) = corners(1, :); % repeat first point to close the hexagon;
corners(end+1, :) = nan; % repeat first point to close the hexagon;


% get shift vectors for the neighboring tiles
% There are six displacement vectors, all have length=spacing
spacing = (hypot(gridAxes(1,1), gridAxes(1,2))*sqrt(3));
% the vectors should be pointing in multiples of 60 deg
[dx, dy] = pol2cart(deg2rad(30:60:360), spacing);

%%% Alternative way of computing shift vectors
% Sort axes by ascending orientation
gridAxes = [gridAxes; -gridAxes]; % flip axes to get other side of hex
orientations = mod(atan2(corners(1:6, 2), corners(1:6, 1)), 2*pi);
[~, isort] = sort(orientations, 'ascend');
% The center of each neighboring tile is specified by adding two adjecent axes 
for i = 1:6
    if i<6
        dx(i) = gridAxes(isort(i), 1)+gridAxes(isort(i+1), 1);
        dy(i) = gridAxes(isort(i), 2)+gridAxes(isort(i+1), 2);
    else
        dx(i) = gridAxes(isort(i), 1)+gridAxes(isort(1), 1);
        dy(i) = gridAxes(isort(i), 2)+gridAxes(isort(1), 2);
    end
end

x = repmat(corners(:, 1), 1, 7)+[0, dx];
y = repmat(corners(:, 2), 1, 7)+[0, dy];

%%% Another alternative that makes it easy to generate more tesselations
%Currently only work for horizontal, pointy-top tilings
% First find the translations for the center row (2
spacing = hypot(originalAxes(1, 1), originalAxes(1, 2));
dx = (-5:5)*spacing;
dx = repmat(dx, 1, 1);
shifts = zeros(11, 1); shifts(1:2:end) = spacing./2; % Add 1 spacing shift to every other row
dx = dx+shifts;

dy = (-5:5)'*spacing./(2*cosd(30))*1.5;
dy = repmat(dy, 1, 11);
% [x, y] = tilePointsHex(zeros(7, 1), zeros(7, 1), gridAxes*2);
% x(end+1, :) = nan; y(end+1, :) = nan;
% x = x(:); y = y(:);

x = repmat(corners(:, 1), 1, 121)+[dx(:)'];
y = repmat(corners(:, 2), 1, 121)+[dy(:)'];

corners = [x(:), y(:)];
centers = [dx(:), dy(:)];
% plot
% hTiles = plot(ax, x(:), y(:), 'Color', col);
end
