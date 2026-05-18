function BW = imregionalmaxStack(I, conn, allowEdgeMax)
% Find regional maxima in a stack of 2D images.
% 
% Returns results identical to imregionalmax for a single 2D image, but
% additionally works with images stacked along the third dimension.
% Super-fast on GPU!

if nargin < 2 || isempty(conn), conn = 8; end
if nargin < 3 || isempty(allowEdgeMax), allowEdgeMax = true; end

if conn==4
    % pixel r,c offsets, clockwise from left
    rshifts = [0, -1, 0, 1];
    cshifts = [-1, 0, 1, 0];
elseif conn==8
    rshifts = [0, -1, -1, -1, 0, 1, 1, 1];
    cshifts = [-1, -1, 0, 1, 1, 1, 0, -1];
else
    error("Value of argument 'conn' must be either 4 or 8'")
end

if allowEdgeMax
    padVal = -inf;
else
    padVal = nan;
end

% Generate vertically and horizontally shifted versions of images. If we
% need to calculate diagonal shifts too, we can start with these stored 
% arrays to avoid recomputation.
IshU = shiftImage(I, 1, -1, padVal);
IshD = shiftImage(I, 1, 1, padVal);
IshL = shiftImage(I, 2, -1, padVal);
IshR = shiftImage(I, 2, 1, padVal);

if isgpuarray(I)
    BW = gpuArray.true(size(I));
else
    BW = true(size(I));
end
    

for sh = 1:numel(cshifts)

    r = rshifts(sh);
    c = cshifts(sh);
    shifts = [r, c];

    if isequal(shifts, [0, -1])
        % Left
        Ish = IshL;
    elseif isequal(shifts, [-1, -1])
        % Up-left
        Ish = shiftImage(IshL, 1, -1, padVal); % shift up
    elseif isequal(shifts, [-1, 0])
        % Up
        Ish = IshU;
    elseif isequal(shifts, [-1, 1])
        % Up-right
        Ish = shiftImage(IshU, 2, 1, padVal); % shift right
    elseif isequal(shifts, [0, 1])
        % Right
        Ish = IshR;
    elseif isequal(shifts, [1, 1])
        % Down-right
        Ish = shiftImage(IshD, 2, 1, padVal); % shift right
    elseif isequal(shifts, [1, 0])
        % Down
        Ish = IshD;
    elseif isequal(shifts, [1, -1])
        % Down-left
        Ish = shiftImage(IshD, 2, -1, padVal); % shift left
    end

    BW = BW & I>Ish;

end

end

function I = shiftImage(I, dim, k, padVal)

% using circshift works much faster than using padarray, because it keeps
% the array size the same
I = circshift(I, k, dim);

if dim==1
    % I = I(idx, :, :);
    if k==1 % D shift
        I(1, :, :) = padVal;
    elseif k==-1 % U shift
        I(end, :, :) = padVal;
    end
elseif dim==2
    % I = I(:, idx, :);
    if k==1 % R shift
        I(:, 1, :) = padVal;
    elseif k==-1 % L shift
        I(:, end, :) = padVal;
    end
end


end
