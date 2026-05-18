function [L, dL, L_poiss, L_prior] = logmargli_gplvm_se_block_ref_nogrid_rg_inp( ...
    uu, xgrid, Fk, yymat, BBwfun, covfun, nf, BBwTfun, ntr, logRinp, FkY, isCirc)

% PGPLVM log-likelihood and gradient for X (latent variables)

uux = reshape(uu,ntr,[],nf);
[ntseg, nneur] = size(yymat);

xxsamp = BBwfun(uux,0); % transform U -> X
xxsamp = cast(xxsamp, 'like', xgrid); % BBwfun can return any numeric class

xxsamp_mt = reshape(xxsamp,[],nf);

% If using GP in "difference" mode, the GP prior is applied to the
% difference between X and a predefined baseline value
if isCirc
    xxsamp_mt = circ_dist(xxsamp_mt, 0);
end

nt = size(xxsamp_mt, 1);
ng = size(xgrid, 1);

% Get weighting of X onto X grid, and dW/dX gradient
[Wxx, dWdxx] = covfun(xgrid,xxsamp_mt);
logR = Wxx'*Fk; % estimated log firing rate

if ~isempty(logRinp)
    logR = logR + logRinp;
end

R = exp(logR);

if isgpuarray(yymat)
    Y3 = reshape(yymat, 1, ntseg, nneur);
    logR3 = reshape(logR, ntseg, 1, nneur);
    YR = squeeze(pagefun(@mtimes, Y3, logR3))'; % 1 x nneur
else
    YR = sum(yymat .* logR);
end

L_poiss = double(YR - sum(R));

% The effect of rhoxx on L happens via this term: when rhoxx increases,
% u decreases in magnitude, thus allowing a larger change in X for a given
% change in L.
L_prior = -.5*(uu'*uu);

% return negative log-likelihood
L_poiss = double(gather(-L_poiss));
L_prior = double(gather(-L_prior));

L = sum(L_poiss) + L_prior;

% Poisson gradient
% Fk = covariance-adjusted F
% FkY = Fk*Y
dLp = FkY + Fk*R';
dLp = repmat(dLp,nf,1) .*  reshape(dWdxx,[],nt);
dLp = reshape(dLp,ng,[]);
dLp = double(gather(sum(dLp,1))); % sum across whole grid
dLp = reshape(dLp,nf,[])'; % (t, f)
dLp = reshape(dLp,ntr,[],nf); % (tr, t, f)
dLp = BBwTfun(dLp,0); % transform X -> U (?)

dL_prior =  double(uu);

dL = dLp(:) + dL_prior(:);

end