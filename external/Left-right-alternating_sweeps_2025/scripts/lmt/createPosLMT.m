function mdl = createPosLMT(Y, D)
% Create LMT for 2-D position

mdl = LmtModel(Y, [D.x, D.y]);
mdl.name = "pos";

% Specify theta phase range to use when aligning the latent and true
% positions
if ~isempty(D.theta)
    phir = [2, 5]; % reset phase
    phi = mod(D.theta - phir(1), 2*pi) + phir(1);
    vphi = phi > phir(1) & phi < phir(2);
    assert(numel(vphi)==mdl.nt);
    mdl.XAlignInds = vphi;
end

% set hparams and initialize

% Clamp pos X for longer, so that other models can remove residual
% variance before it's allowed to leak into pos X.
mdl.XclampIter = 40;

mdl.hparams.rhoff = 0.1;
mdl.hparams.rhoxx = 10; % 1000
mdl.hparams.lamff = 1; % remove this for single-module stuff
mdl.hparams.lenffR = 6;

mdl.hparamRanges.sigma = [3, 1, 100];
mdl.XboundPercentile = 0.01;

end