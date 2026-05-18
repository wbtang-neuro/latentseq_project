function [fspaces, decompX, info] = createDefaultFunctionSpaces(boxLims, useGpu, P)

if nargin < 2 || isempty(useGpu), useGpu = false; end
if nargin < 3 || isempty(P), P = glmParams(); end

info = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position (only if a box size is given)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin && ~isempty(boxLims)
    [fspaces.pos, decompX.pos, info.basisPosMu] = createPosFunctionSpace(boxLims, useGpu, P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Angular variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa = P.basisAngKappa;
mu = (1 : P.basisAngN)' / P.basisAngN * 2*pi;
mu = mu-pi; % added 2024-11-20 to fix plotting issues
mu = cast(mu, P.floatClass);
kappa = cast(kappa, P.floatClass);
F = createFcnSpace('vm', {mu, kappa});

nTh = P.basisAngN * 100;
th = (1:nTh)' / nTh * 2*pi;
th = th-pi; % added 2024-11-20 to fix plotting issues
if useGpu
    th = gpuArray(th);
end
th = cast(th, P.floatClass);

[FR, D] = F.dimReduce(th, SVD());
FR.rescale(1./std(D.Y));
FR.discardThresh(P.pcaThresh.hd, 'above', 2);
fspaces.hd = FR.copy();
fspaces.id = FR.copy();

FR = F.dimReduce(th, SVD());
FR.rescale(1./std(D.Y));
FR.discardThresh(P.pcaThresh.theta, 'above', 2);
fspaces.theta = FR.copy();

decompX.hd = th;
decompX.theta = th;

end