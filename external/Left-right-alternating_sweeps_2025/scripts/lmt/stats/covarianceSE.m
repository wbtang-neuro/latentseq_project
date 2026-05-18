function [A, dA] = covarianceSE(len, rho, x, z)
% Squared-exponential covariance function gradient
if nargin == 3, z = x; end

persistent dA3 A2

[~, D] = size(x);

Mx = x/len;
Mz = z/len;

useGpu = isgpuarray(x);

if isequal(Mx, Mz)
    % save on computation for calculating within-grid distances
    xD0 = pdist(Mx); % returns lower triangular matrix
    D = squareform(xD0);
    D = cast(D, "like", Mx);
    A = rho .* exp(-D.^2);
else
    if useGpu
        A2 = arrayfun(@(x1, x2, y1, y2) exp( - ( (x1-x2).^2 + (y1-y2).^2 ) ), Mx(:, 1), Mz(:, 1)', Mx(:, 2), Mz(:, 2)');
    else
        x1 = Mx(:, 1);
        x2 = Mz(:, 1)';
        y1 = Mx(:, 2);
        y2 = Mz(:, 2)';
        A2 = exp( - ( (x1-x2).^2 + (y1-y2).^2 ) );
    end
    A = rho .* A2;
end

if nargout > 1
    nx = size(x,1);
    nz = size(z,1);
    
    x = reshape(x, [nx, 1, D]);
    x = permute(x, [1, 3, 2]);
    
    zr = reshape(z, [1, nz, D]);
    zr = permute(zr, [1, 3, 2]);
    lensq = cast(len^2, "like", x);

    A3 = reshape(A2, nx, 1, nz);
    if useGpu
        dA3 = arrayfun(@(a, x, zr, e) a .* (x-zr) ./ e, A3, x, zr, lensq); % for gpuArray
    else
        dA3 = A3 .* (x-zr) ./ lensq;
    end
    dA = reshape(dA3, nx, D*nz);

end

    
end
