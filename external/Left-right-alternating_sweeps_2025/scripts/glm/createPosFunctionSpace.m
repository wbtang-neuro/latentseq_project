function [fspace, decompX, basisPosMu] = createPosFunctionSpace(boxLims, useGpu, P)
if nargin < 2 || isempty(useGpu), useGpu = false; end
if nargin < 3 || isempty(P), P = glmParams(); end

% 2022-11-17 fixed mixup of decomp and basis func grid lims. The latter
% should be larger! But it wasn't before...

assert(isequal(size(boxLims), [2, 2]) && all(issorted(boxLims)), ...
    'Argument "boxLims" must be a 2*2 matrix in form: [Xl, Yl; Xh, Yh]');

% N.B. the supplied "boxLims" indicates the extent of the basis-func grid.
% The parameter "basisPosBoxPadding" specifies how far the *decomposition*
% grid extends beyond the edge of the BF grid.

lims = num2cell(boxLims, 1);
limsPad = cellfun(@(x) {x+[-1; 1]*P.basisPosBoxPadding}, lims); 
% basBox(:, 1) = limsBasis{1}([1, 1, 2, 2]);
% basBox(:, 2) = limsBasis{2}([1, 2, 2, 1]);
% inBasBox = inpolygon(gridCoords(:, 1), gridCoords(:, 2), basBox(:, 1), basBox(:, 2));
% basisPosMu = gridCoords(inBasBox, :);

sigma = P.basisPosSigma;
decompStep = P.basisPosDecompGridStep;
% if useGpu
%     sigma = gpuArray(sigma);
%     decompStep = gpuArray(decompStep);
%     basisPosMu = gpuArray(basisPosMu);
% end
gridDecomp = cellfun(@(x) {x(1): decompStep : x(2)}, lims);
[iiD, jjD] = ndgrid(gridDecomp{:});
decompX = [iiD(:), jjD(:)]; % y, x
decompX = cast(decompX, P.floatClass);

% Basis function grid - this extends slightly outside of the decomposition
% grid.
boxCenter = cellfun(@mean, lims);
basBoxSize = cellfun(@diff, limsPad);
nRings = max(basBoxSize) / P.basisPosSpacing;
basisPosMu = gridNodes(boxCenter, nRings, P.basisPosSpacing, 0);
D = pdist2(basisPosMu, decompX);
% dmin = min(D, [], 2);
box(:, 1) = lims{1}([1, 1, 2, 2]);
box(:, 2) = lims{2}([1, 2, 2, 1]);
validBas = inpolygon(basisPosMu(:, 1), basisPosMu(:, 2), box(:, 1), box(:, 2));
% basisPosMu = gridCoords(inBasBox, :);
% validBas = dmin <= (P.basisPosBoxPadding + 1e-5);
basisPosMu = basisPosMu(validBas, :);

F = createFcnSpace('gaus2', {basisPosMu, sigma*[1, 1]});
[FR, D] = F.dimReduce(decompX, SVD());
FR.rescale(1./std(D.Y));
% [FR, D] = F.dimReduce(decompX, PCA());
% FR.rescale(1./std(D.Y, [], 'all'));
FR.discardThresh(P.pcaThresh.pos, 'above', 2);
fspace = FR;

end

function [Points] = gridNodes(phase, nRings, spacing, orientation)
% GRIDNODES - generate XY coords of grid nodes
%
% This is an adaptation of Tor's code

% z = re^{i*theta}; x = Re|z|, y = Im|z|

% Example run: z = TSgridNodes3 (15, 20, 30);

% nRings = number of rings around center field
% spacing --- is spacing.
% orientation, given in radians

rInd = 1:nRings;
nInd = 6*nRings; % Total number of nodes, increases by 6 per ring
phi = spacing;

z = [];

for step = 1:length(rInd)
    
    a = rInd(step); % Ring number = unit distance along x-axis
    if nInd == 0 || step ~= length(rInd)
        b = 0:6*a-1; % # nodes in ring
    else
        b = 0:nInd-1;
    end
    
    c = mod(b,a);
    theta = atan((sqrt(3).*c)./(2*a-c)) + pi*(b-c)/(3*a) + orientation; % RG 21/03/2015
    r = phi .*sqrt((a-c).^2 +a.*c);
    z = [z r.*exp(1i .*theta)];
    
end

x = real(z)+phase(1);
y = imag(z)+phase(2);

Points = phase;
Points = [Points; x' y'];

end