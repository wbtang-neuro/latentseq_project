function [C,numberOfOverlapPixels] = normxcorr2_general_stack(T, A, nbatch, useGpu)
% Adapted to handle 3D arrays of 2D images stacked in third dimension.
% Runs in batches to minimize memory usage.

nstack = size(A, 3);
if nargin < 3 || isempty(nbatch), nbatch = min(1e3, nstack); end
if nargin < 4 || isempty(useGpu), useGpu = true; end

c = 0;
C = zeros(0, 0, nstack, underlyingType(T));

if useGpu
    T = gpuArray(single(T));
end

while c<nstack
    iend = min(c+nbatch, nstack);
    istack = c+1 : iend;
    As = A(:, :, istack);
    if useGpu
        As = gpuArray(single(As));
    end
    [Cs, numberOfOverlapPixels] = normxcorr2_general_stack_helper(T, As);
    if c==0
        C = zeros([size(Cs, [1, 2]), nstack], underlyingType(A));
    end
    C(:, :, istack) = gather(Cs);
    % C = cat(3, C, gather(Cs));
    c = c + nbatch;
end


end
