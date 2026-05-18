function [L, Lall, pos] = posShiftCostFcn(params, spikeCounts, shiftData, interpD, prmInds, ...
    XA, XB, checkLlh, alphaPenalty, betaPenalty, gpuInterpClass)
% Optimization cost function for position-shifted GLM
%
% This cost function allows optimization for the GLM coefficients (beta)
% together with an extra set of parameters (alpha) which specifies the
% theta-modulated shifting of the animal's tracked position.
%
% N.B. the alpha parameters introduce a nonlinearity which mean that this
% isn't a classical GLM, therefore making the optimization more tricky.
%
% params    : vector of all parameters (alpha and beta)
%
% spikes    : spike counts (numeric vector)
%
% pos       : 2D position coordinates (2-column matrix, [X, Y])
%
% hd        : azimuth head direction (vector)
%
% DPosI     : position decomposition for interpolation (DRDecomp object)
%
% prmInds   : structure specifying indices of grouped parameters in the
%             'params' input. This must be a structure with the following
%             fields:
%                   'beta': this must be a structure where each field
%                   specifies the parameter indices of a group of GLM beta
%                   covariates, e.g. 'intercept', 'pos' or 'hd'.
%                   'alpha': this must be a structure with each field
%                   representing a group of alpha parameters.
%
%                   'alpha': similarly, this is a structure where each
%                   field specifies the indices of an alpha parameter
%                   group in the 'params' vector. The possible alpha
%                   parameter groups are:
%                       1) 'intercept', which defines the baseline AP
%                          position shift value, and
%                       2) additional groups which specify AP-shifts
%                          dependent on a group of GLM covariates (i.e. a
%                          group defined as a subfield of the 'beta' field
%                          in prmInds.
%
% XA        : matrix containing alpha covariates. Dimensions are
%             [time, covariates]. The number of covariates must
%             match the number of indices specified in prmInds. The AP
%             shift value for each time bin is calculated by
%             multiplying XA by the vector of alpha parameters. The order
%             of columns in XA must be consistent with the order of the
%             corresponding parameters. If there is an alpha intercept
%             term, there must be a column of ones in XA.
%
% XB        : matrix containing beta covariates (GLM design matrix).
%             Dimensions are [time, covariates].
%
% dt        : size of time bins
%
% limitToBox: specifies whether shifted positions should be 'cropped' if
%             they are shifted outside of the boundary of the box
%
% box       : coordinates of the box corners, as a 2-column matrix [x, y]

persistent alpha0 X0 firstRun SD XBTmp

% Called without args: clear persistent vars
if nargin == 0 || isempty(firstRun)
    firstRun = true;
    if nargin == 0
        return
    end
end

useGpu = isa(XB, 'gpuArray');
% gpuInterpClass = 'single';
% gpuInterpClass = underlyingType(XB);

if firstRun && ~isempty(shiftData)
    SD = shiftData;
    xg = double(interpD.XGrid{1});
    SD.gridStep = diff(xg([1, end])) / (numel(xg)-1);
    if useGpu
        SD.DPosIGpu = interpD.Y'; % should already be gpuArray
    end
end

params = params(:);

% interpGrid = interpD.XGrid;
% % if useGpu
% %     interpGrid = cellfun(@(x) {gpuArray(cast(x, gpuInterpClass))}, interpGrid);
% % end

alpha = params(prmInds.alpha);
beta = params(prmInds.beta);

% alpha = extractParams(params, prmInds.alpha);
% beta = extractParams(params, prmInds.beta);
doShift = ~isempty(alpha);

% Now we need to convert the shifted X/Y coords into the appropriate form
% for submitting to the GLM. The expensive process of decomposing the
% shifted positions can be avoided by instead reading off values from a
% precalculated interpolation grid.

pos = SD.pos;

if doShift
    
    % Check if we need to update XB for the position covariate, because either
    %   1) The shift parameters have changed
    %   2) The input dataset has changed (indicated by "firstRun" flag)
    if firstRun || ~isequal(alpha0, alpha)
        
        % Apply positional offsetting based on alpha params
        if ~isempty(alpha)
            deltaDist = XA * alpha;
            % shift by specified distance along HD line
            pos = posDistShift(pos, SD.angle, deltaDist);
        end
        
        icolXB = prmInds.betaPos;
        if useGpu
            if firstRun
                nc = numel(icolXB);
                nr = size(XB, 1);
                XBTmp = gpuArray.zeros(nc, nr, 4, gpuInterpClass);
            end
            posTmp = cast(pos, gpuInterpClass);
            XBTmp = fastBlerpGeneral(interpD.XGrid, SD.DPosIGpu, posTmp, XBTmp, SD.gridStep);
            XB(:, icolXB) = sum(XBTmp, 3)';
        else
            fastBlerpCpu(interpD.XGrid, interpD.Y, pos, XB, icolXB, SD.gridStep);
        end
        
%         if checkLlh
%             sz = cellfun(@numel, grid);
%             nc = size(interpD.Y, 2);
%             XB2 = zeros(nt, nc, 'like', XB);
%             for n = 1:nc
%                 Y = reshape(interpD.Y(:, n), sz);
%                 XB2(:, n) = interpn(grid{2}, grid{1}, Y, pos(:, 1), pos(:, 2));
%             end
%             if isempty(fig) || ~isvalid(fig)
%                 fig = figure();
%                 for n = 1:3
%                     subplot(3, 1, n, 'parent', fig);
%                 end
%             end
%             axs = findobj(fig, 'type', 'axes');
%             nRowPlt = min(1000, nt);
%             z1 = XB(1:nRowPlt, icolXB)';
%             z2 = XB2(1:nRowPlt, :)';
%             imgs = {z1, z2, z1-z2};
%             titles = {'Quick', 'Interpn', 'Difference'};
%             for n = 1:3
%                 ax = axs(n);
%                 imagesc(ax, 1:n, icolXB, imgs{n});
%                 xlabel(ax, 'Sample');
%                 ylabel(ax, 'Covariate #');
%                 colorbar(ax);
%                 title(ax, titles{n});
%             end
%             drawnow();
%         end
        
        X0 = XB;
        alpha0 = alpha;
        newX = true;
    else
        XB = X0;
        newX = false;
    end
else
    newX = false;
end

% Use the "persistent" GLM LLH function - it will only update beta values
% that have changed.
clearMem = firstRun || newX;
[L, Lall] = poissonRegLlhPersist(XB, spikeCounts, beta, clearMem);
if useGpu
    L = gather(L);
end
L = double(L);

if checkLlh
    llh2 = poissonRegLlh(XB, spikeCounts, beta, -1); % not normalized for length
    llh2 = gather(llh2) ./ numel(spikeCounts);
    llhDiff = L - llh2;
    if abs(llhDiff) > eps
        warning('LLH calculation may be incorrect: discrepancy of %.6e', llhDiff);
    end
end

% Apply L2 penalty terms for alpha and beta. LLH is already normalized for
% the length of the data, so we don't need to adjust the penalties for this.
% apen = penaltyTerm(alpha, alphaPenaltyType, alphaPenalty);
apen = penaltyTerm(alpha, 'L2', alphaPenalty);
bpen = penaltyTerm(beta(2:end), 'L2', betaPenalty);

L = L + apen + bpen;
firstRun = false;

end

function p = penaltyTerm(params, penaltyType, lambda)

if lambda
    if strcmpi(penaltyType, 'L1')
        psum = sum(abs(params));
    elseif strcmpi(penaltyType, 'L2')
        psum = sum(params.^2);
    else
        error('Penalty type must be "L1" or "L2"');
    end
    p = lambda * psum;
else
    p = 0;
end

end

function VOut = fastBlerpGeneral(grid, VIn, iX, VOut, gridStep)
[w, ir] = interp2Weights(grid, iX, gridStep);
VOut(:,:) = VIn(:, ir(:)) .* w(:)';
end

function VOut = fastBlerpCpu(grid, VIn, iX, VOut, icolumn, gridStep)

% grid: 2-element cell array containing image (y,x) pixel grid positions
% V: nt*nimage matrix

% Interpolate decomposed positions at the shifted position coordinates.
%
% Nearest-neighbor interpolation would be the most efficient, but this
% seems to cause problems for the optimization, so use bilinear
% instead, to avoid discrete jumps in the interpolated values.

ic = icolumn;
nc = numel(ic);
[w, ir] = interp2Weights(grid, iX, gridStep);

for n = 1:4
    rezero = n==1;
    colIndsInC = uint32(0 : nc-1);
    colIndsOutC = uint32(ic(:)-1)';
    rowIndsInC = uint32(ir(:, n)-1);
    if isa(VIn, "double")
        accumulateInPlace(VIn, VOut, colIndsInC, colIndsOutC, rowIndsInC, w(:, n), rezero);
    else
        accumulateInPlace_single(VIn, VOut, colIndsInC, colIndsOutC, rowIndsInC, w(:, n), rezero);
    end
end

end

function [w, inds] = interp2Weights(grid, iX, gridStep)

% grid: 2-element cell array containing image (y,x) pixel grid positions
% V: nt*nimage matrix

% Interpolate decomposed positions at the shifted position coordinates.
%
% Nearest-neighbor interpolation would be the most efficient, but this
% seems to cause problems for the optimization, so use bilinear
% instead, to avoid discrete jumps in the interpolated values.

% iX = gather(iX);

nGrid = cellfun(@numel, grid);
% gridRes = cellfun(@(x) x(2)-x(1), grid);

ni = size(iX, 1);

% For each 2-D point to be interpolated, calculate the indices of the
% nearest upper- and lower-bound neighbouring pixels, and how they
% should be weighted
subs = zeros(ni, 2, 2, 'like', iX);
giL = zeros(ni, 2, 'like', iX);

for d = 1:2
    iL = gridnn(grid{d}, iX(:, d), gridStep);
%     gcpu = gather(grid{d});
%     iL = gridnn(grid{d}, iX(:, d), 'floor');
    iU = iL+1;
    iU = min(iU, nGrid(d));
%     iU(iU > nGrid(d)) = nGrid(d); % limit upper bound to the grid size
    giL(:, d) = grid{d}(iL);
    subs(:, d, 1) = iL;
    subs(:, d, 2) = iU;
end

res = (iX - giL) ./ gridStep; % residual after subtracting lower bound
% res = (iX - giL) ./ gridRes(d(:))'; % residual after subtracting lower bound

% inds = uint32(zeros(ni, 4, 'like', iX));
inds = zeros(ni, 4, 'like', iX);
w = zeros(ni, 4, 'like', iX);
c = 0;

% Extract the values of the calculated upper-/lower-bound neighbour indices and weights
for i = 1:2 % Dim 1 (row)
    si = subs(:, 1, i);         % row subscripts of neighbor pixels (iter 1 = lower, iter2 = upper)
    wi = res(:, 1);             % weights for neighbor pixels
    if i==1, wi = 1-wi; end     % (convert lower/upper bound weights)
    
    for j = 1:2 % Dim 2 (column)
        c = c+1;
        % Get the column subscripts and weights (just as for rows)
        sj = subs(:, 2, j);
        wj = res(:, 2);
        if j==1, wj = 1-wj; end
        % Convert r/c subscripts to 1-D indices
        inds(:, c) = (sj-1) * nGrid(1) + si;
        % Combine row and column weights
        w(:, c) = wi .* wj;
    end
end

end