function [L, dL] =  StateSpaceModelsofSpikeTrains_ref_rg_inp( ...
    F, Y, invC, wX, logRInp, lambda)

useGpu = isgpuarray(F);

% PGPLVM cost function for F

[nt, nneur] = size(Y);
F2 = reshape(F,[],nneur); % log tuning curves (bin, neur)
logR = wX*F2; % all neuron predicted firing rates

% get full prediction of Y by adding external input
if ~isempty(logRInp)
    logR = logR + logRInp;
end

% Poisson log-likelihood
R = exp(logR);
if useGpu
    Y3 = reshape(Y, 1, nt, nneur);
    logR3 = reshape(logR, nt, 1, nneur);
    YR = squeeze(pagefun(@mtimes, Y3, logR3))';
else
    YR = sum(Y .* logR);
end

L_poiss = double(YR - sum(R));

% Quadratic term (gaussian prior)
F2d = double(F2);
invCd = double(invC);
L_gp = -.5*trace(invCd*(F2d*F2d')); % Anqui's original calculation

% L1 penalty
vH = F2 > 0;
vL = F2 < 0;
Flam = lambda .* F2d;
L_L1H = -sum(Flam(vH));
L_L1L = sum(Flam(vL));
L_L1 = L_L1H + L_L1L;

% LLH
L_poiss = gather(-L_poiss);
L_gp = gather(-L_gp);
L_L1 = gather(-L_L1);
L = sum(L_poiss) + L_gp + L_L1;

% GRADIENT
dL_pois = double(wX'*(Y-R));

% GP cost term
dL_gp = -invCd*F2d; % Original

% L1 penalty forces F values towards the mean
dL_L1 = zeros(size(dL_gp), "like", F);
lambdaMat = dL_L1 + lambda; % lambda can be scalar, vector or maxtrix

dL_L1(vH) = -lambdaMat(vH);
dL_L1(vL) = lambdaMat(vL);

dL = dL_pois + dL_gp + dL_L1;

dL = gather(-dL(:));

end