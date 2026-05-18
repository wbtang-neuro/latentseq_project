function [L, dL] = poissonLlh(X, Y, beta, returnSign, inpR, lambdaL1)
% Calculate log-likelihood gradients for multivariate Poisson regression

% X: regression design matrix (observations * variables)
% Y: response variable (observations * variables), as a vector
% beta: coefficients matrix (ncoeff * variables)
%
% returnSign: either +1 or -1, specifies whether to return LLH and
% inpR: input log-firing-rates (optional, leave empty if not used)
% derivatives with original sign (+1), or inverted sign (-1)

if nargin < 4 || isempty(returnSign), returnSign = 1; end
if nargin < 5, inpR = []; end
if nargin < 6 || isempty(lambdaL1), lambdaL1 = 0; end

useGpu = isgpuarray(X);

if useGpu
    beta = gpuArray(beta);
end

% Calculate LLH per spike
nt = size(X, 1);
nneur = numel(Y) / nt;
nbeta = size(X, 2);

beta = reshape(beta, nbeta, nneur);
E = X * beta;
E = E(:);

if ~isempty(inpR)
    E = E + inpR(:);
end

Yeta = Y'*E;
E = exp(E);
L = Yeta - sum(E); % LL

dY = reshape(Y-E, nt, nneur);
dL = X' * dY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate L1 / L2 penalty terms

Lpenalty = 0;
dLpenalty = 0;

if any(lambdaL1(:))
    [Lpenalty, dLpenalty] = calcL1Penalty(beta, lambdaL1);   
end

L = L + Lpenalty;
dL = dL + dLpenalty;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = gather(L);
dL = gather(dL(:));

if returnSign == -1
    L = -L;
    dL = -dL;
end

end

function [L, dL] = calcL1Penalty(beta, lambda)
vH = beta > 0;
vL = beta < 0;
blam = lambda.*beta;
L_L1H = -sum(blam(vH));
L_L1L = sum(blam(vL));
L = L_L1H + L_L1L;

dL = zeros(size(beta)) + lambda; % lambda can be scalar or vector
lambdaVec = dL + lambda;
dL(vH) = -lambdaVec(vH);
dL(vL) = lambdaVec(vL);
end
