function [L, Lall] = poissonRegLlhPersist(X, y, beta, clearMem)
% Version with memory of previous evaluations
%
% X: regression design matrix (observations * variables)
% y: response variable (column vector)
% yHat: prediction of y

if nargin < 4 || isempty(clearMem), clearMem = false; end

% Assume that X and y are fixed
persistent beta0 eta0 etaC Lall0

ny = numel(y);
useGpu = isa(X, 'gpuArray');

if clearMem || ~isequal(size(beta0), size(beta))
    % Initialize persistent vars
    betaUpdated = true(size(beta));
%     if useGpu
%         betaUpdated = gpuArray(betaUpdated);
%     end
else
    betaUpdated = beta0 ~= beta;
end

% betaUpdated(:) = true;

% The optimization will alter individual beta components
% separately, so it's usually not necessary to peform the full X*beta
% multiplication. If only a small number of beta elements change, it's
% much more efficient to only multiply the updated covariates and
% betas, accumulating the result into a running estimate of eta.

if all(betaUpdated)
    % All betas changed, do the full multiplication
    etaC = X*beta;
    eta0 = etaC;
    beta0 = beta;
else
    % Subset of betas changed, use in-place matrix-vector
    % multiplication to update eta0 and avoid expensive copying
    idx = find(betaUpdated);
    if isempty(idx)
        % Input is identical to last time: return stored Lall and L
        Lall = Lall0;
        L = sum(Lall) / ny;
        return;
    else
        w = beta(idx)-beta0(idx);
        if useGpu
    %         etaC = eta0 + w .* X(:, betaUpdated)
            XTmp = X(:, idx);
            Lall = arrayfun( @(y, x, eta0, w) double(exp(eta0 + w.*x) - y.*(eta0 + w.*x) ), y, XTmp, eta0, w);
            L = sum(Lall) / ny;
            Lall0 = Lall;
    %         L = sum( arrayfun( @(y, x, eta0, w, n) double(exp(eta0 + w.*x) - y.*(eta0 + w.*x) ) / n, y, XTmp, eta0, w, ny ) );
            return;
        else
            etaC(:) = eta0; % copy values, not reference
            indsC = uint16(idx-1);
            if isa(X, "double")
                mvmInPlace(X, etaC, w, indsC);
            else
                mvmInPlace_single(X, etaC, w, indsC); % N.B. 'w' is still double
            end
        end
    end
end

if useGpu
%     Lall = arrayfun(@(y, eta, n) double(exp(eta) - y.*eta) / n, y, etaC, ny);
    Lall = arrayfun(@(y, eta) double(exp(eta) - y.*eta), y, etaC);
else
    yHat = exp(double(etaC));
    Lall = yHat - double(y.*etaC);
end

% Total L is normalized by number of time bins
Lall0 = Lall;
L = sum(Lall) / ny;

% if useGpu
%     L = sum(arrayfun(@(y, eta, n) double(exp(eta) - y.*eta) / n, y, etaC, ny) );
% else
%     yHat = exp(etaC);
%     L = sum(yHat - y.*etaC) ./ ny;
% end


end