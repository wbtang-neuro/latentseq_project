function himg = plotAgent(x, y, P)

arguments
    x
    y
    P.agentType = "robot"
    P.height = nan
    P.width = nan
    % P.alignTo = "head"
    P.alpha = 1
    P.rotationDeg = 0      % not currently implemented
    P.axes = []
end

[x, y, P.height, P.width, P.rotationDeg] = parseParams(x, y, P.height, P.width, P.rotationDeg);
% n = numel(x);

if strcmpi(P.agentType, "robot")
    pngFile = "NicePng_robot-clipart-png_2940483.png";
    [z, ~, alpha] = imread(pngFile);
    refPoint = [0, 0.2];
elseif strcmpi(P.agentType, "rat")
    pngFile = "Rat_Top_by_GC.png";
    [z, ~, alpha] = imread(pngFile);
    z = rot90(z, -1);
    alpha = rot90(alpha, -1);
    refPoint = [0.5, 0];
end

z = flipud(z);
alpha = flipud(alpha) * P.alpha;
sz = size(z);

% npad = ceil(max(sz)/3);
% z = padarray(z, [npad, npad], 0);
% alpha = padarray(alpha, [npad, npad], 0);
% for i = 1:n
%     [zR{i}, alphaR{i}, szR(i, :)] = rotateImage(z, alpha, P.rotationDeg(i));
% end

if ~any(isnan(P.height))
    P.width = P.height .* sz(:,2)./sz(:,1);
elseif ~any(isnan(P.width))
    P.height = P.width .* sz(:,1)./sz(:,2);
else
    error("Either height or width must be specified");
end

n = numel(x);

ax = P.axes;
if isempty(ax)
    ax = gca;
end

for i = 1:n
    [xv, yv] = getAlignedCoords(x(i), y(i), P.width(i), P.height(i), sz, refPoint, P.rotationDeg(i));
    himg = image(ax, xv, yv, z, "alphaData", alpha);
end

end

% function [z, alpha, sz] = rotateImage(z, alpha, rot)
% 
% if rot ~= 0
%     z = imrotate(z, rot, "crop");
%     alpha = imrotate(alpha, rot, "crop");
% end
% sz = size(z, [1, 2]);
% 
% end

function varargout = parseParams(varargin)
n = max(cellfun(@numel, varargin));
varargout = varargin;
for i = 1:nargin
    arg = varargin{i}(:);
    if isscalar(arg)
        arg = repmat(arg, n, 1);
    end
    varargout{i} = arg;
end
end

function [xv, yv] = getAlignedCoords(x, y, w, h, sz, refPoint, rotationDeg)
% offsets don't yet account for effect of padding/cropping with rotation

xv0 = linspace(-1/2, 1/2, sz(2));
yv0 = linspace(-1/2, 1/2, sz(1));

% switch lower(alignTo)
%     case "eye"
%         offset = [0.32, 0.22];
%     case "antenna"
%         offset = [0, 0.48]; % unpadded
%         % offset = [0, 0.3];
%     case "head"
%         offset = [0, 0.2];
% end

% if rotationDeg ~= 0
    % refPoint = rotate2d(refPoint, -deg2rad(rotationDeg));
% end

% align the image such that (0, 0) is aligned to the target body part
xv0 = xv0 - refPoint(1);
yv0 = yv0 - refPoint(2);

% scale image to specified height and width
xv0 = xv0*w;
yv0 = yv0*h;

% shift the image to the specified location
xv = xv0 + x;
yv = yv0 + y;

end