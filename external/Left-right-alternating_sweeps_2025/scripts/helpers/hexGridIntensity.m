function z = hexGridIntensity(x, y, varargin)
% HEXGRIDINTENSITY construct hexagonal grid intensity function

inp = inputParser();
inp.KeepUnmatched = false;
inp.addParameter('orientation', 0,      @(x) validateattributes(x, {'double'}, {'scalar'}));
inp.addParameter('phase', [0 0],        @(x) validateattributes(x, {'double'}, {'size', [1 2]}));
inp.addParameter('fieldWidth', 0.2,     @(x) validateattributes(x, {'double'}, {'scalar', 'nonnegative'} ));
inp.addParameter('range', 2,            @(x) validateattributes(x, {'double'}, {'scalar', 'increasing'} ));
inp.addParameter('spacing', 1,          @(x) validateattributes(x, {'double'}, {'scalar', 'nonnegative'} ));
inp.addParameter('pdfType', 'gaussian', @(x) validateattributes(x, {'char', 'string'}, {'nrows', 1}));
inp.parse(varargin{:});
P = inp.Results;

[xx, yy] = meshgrid(x, y);
z = zeros(size(xx));

if strcmpi(P.pdfType, 'interference')
    % Create PDF as an interference pattern of cosine waves

    for n = 1:3
        ang = 2*pi * n/6;
        xxR = (xx-P.phase(1))*cos(ang) + (yy-P.phase(2))*sin(ang);
        z = z + (1 + cos(2*pi * xxR/P.spacing) );
    end
    z = z-min(z(:));
else
    % Create PDF by summing individual field PDFs
    
    % Wrap phase onto rhombus. The rhombus wrapping doesn't take grid
    % orientation into account, so must unrotate and rerotate phase
    
    phase = rotate2d(P.phase, -P.orientation);
    phase = wrapPhaseRhombus(phase / P.spacing) * P.spacing;
    phase = rotate2d(phase, P.orientation);
    
    % Get node coords
    maxrng = max( max(abs(x)), max(abs(y)) );
    nRings = ceil((1.5*maxrng + 10*P.fieldWidth) / P.spacing);
    nodes = hexGridNodes(phase, nRings, P.spacing, P.orientation);
    
    
    for f = 1:length(nodes(:,1))
        node = nodes(f, :);
        zField = getField(xx, yy, node(1), node(2), P.fieldWidth, P.pdfType);
        z = z + zField;
    end
    
end

z = z./sum(z(:));

end

function z = getField(xx, yy, cx, cy, width, type)

switch lower(type)
    case 'gaussian'
        % Width specifies standard deviation
        dd = hypot(xx-cx, yy-cy);
        z = normpdf(dd, 0, width);
        % z = rg.helpers.normpdf2(xx, yy, [cx cy], width);
    case 'tophat'
        % Width specifies step cut-off
        z = zeros(size(xx));
        dist = hypot(xx-cx, yy-cy);
        inField = dist<width;
        z(inField) = 1;
    case 'cone'
        % Width specifies edge of cone
        dist = hypot(xx-cx, yy-cy);
        z = max(0, 1-(dist/width));
    case 'cosine'
        % Width specifies edge of rectified cosine
        dist = hypot(xx-cx, yy-cy);
        z = cos(dist/width * (pi/2));
        % Truncate the function to zero at the field boundary
        z(dist > width) = 0;
end

end

function pointsW = wrapPhaseRhombus(points, orientation, origin)
% WRAPPHASERHOMBUS wraps x-y coordinates onto rhombus grid phase tile
%
% PW = WRAPPHASERHOMBUS(P) wraps the x/y points in P such that they fall
% within the rhombus phase tile of a grid with spacing 1. P must be a
% 2-column matrix where each row contains the x- and y-phase normalized to
% grid spacing (i.e. spacing = 1)
%
% The unit-rhombus corner positions are defined as below:
% [0,    0]
% [0.5, cos(pi/6)]
% [1.5, cos(pi/6)]
% [1,    0]
%

if nargin < 2 || isempty(orientation), orientation = 0; end
if nargin < 3 || isempty(orientation), origin = [0, 0]; end

r = tan(pi/3);
x = points(:, 1);
y = points(:, 2);

[x, y] = rotate2d(x, y, -orientation);
x = x + origin(1);
y = y + origin(2);

nYCycles = floor(y/(r/2));
yW = y - (nYCycles*(r/2)); % use the rounding already defined by floor()

x = x + nYCycles/2;
xShift = mod(yW/r, 1); % - floor(y/cos(pi/6));
xW = mod(x-xShift, 1) + xShift;

xW = xW - origin(1);
yW = yW - origin(2);
[xW, yW] = rotate2d(xW, yW, orientation);

pointsW = [xW, yW];

end