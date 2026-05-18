function [A, dA] = covarianceSECirc(len, rho, x, z)
% Squared-exponential covariance function gradient with circular wrapping
if nargin == 3, z = x; end

persistent dA2 A2

Asq = (circ_dist(x, z') ./ len) .^ 2;
A = rho .* exp(-Asq);

if nargout > 1
    [nx, D] = size(x);
    nz = size(z,1);
    % [~, D] = size(x);
    
    x = reshape(x, [nx, 1, D]);
    x = permute(x, [1, 3, 2]);
    
    zr = reshape(z, [1, nz, D]);
    zr = permute(zr, [1, 3, 2]);
    lensq = cast(len^2, "like", x);
    A2 = reshape(A, nx, 1, nz);

    if isgpuarray(x)
        dA2 = arrayfun(@(a, x, zr, e) a .* circ_dist(x, zr) ./ e, A2, x, zr, lensq);
    else
        dA2 = A2 .* circ_dist(x, zr) ./ lensq;
    end

    dA = reshape(dA2, nx, D*nz);

end

end