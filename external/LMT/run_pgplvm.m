function [result_la,setopt]=run_pgplvm(xx,yy,tgrid,nf,niter,xppca)

nt = numel(tgrid);
if (size(tgrid,1)~=nt)||(size(yy,1)~=nt)||(size(xx,1)~=nt)
    error(['Dimensions are not cool! tgrid, xx, and yy should have ' num2str(nt) ' rows.'])
end
    
% tgrid = [1:nt]';
    
setopt.tgrid = tgrid;
setopt.latentTYPE = 1; % kernel for the latent variable, 1. AR1, 2. SE

%% == 1. Compute baseline estimates ====

% Initialize the log of spike rates with the square root of spike counts.
ffmat = sqrt(yy);

% % Compute PPCA and show xx initialization
% if isempty(xppca)
%     xppca = pca(ffmat,nf);
% end
% setopt.xplds = xppca;%[xx xx];%xx; % %for initialization purpose

xplds = run_plds(yy,nf)';
setopt.xplds = xplds; % for initialization purpose

if nf==1
    xppcamat = align_xtrue(xppca,xx);
    setopt.xpldsmat = xppcamat;%[xx xx]; %xx; %  for plotting purpose
end

%% == 2. Compute P-GPLVM ====
% Set up options

setopt.lr = 0.95; % learning rate
setopt.ffTYPE = 2; % kernel for the tuning curve, 1. AR1, 2. SE

setopt.initTYPE = 1; % initialize latent: 1. use PLDS init; 2. use random init; 3. true xx
setopt.opthyp_flag = 1;
setopt.sigma2_init = 2; % initial noise variance
setopt.rhoxx = 10; % rho for Kxx
setopt.lenxx = 50; % len for Kxx
setopt.rhoff = 10; % rho for Kff
setopt.lenff = 50; % len for Kff


setopt.la_flag = 3; % 1. no la; 2. standard la; 3. decoupled la
setopt.hypid = [1,2,3,4]; % 1. rho for Kxx; 2. len for Kxx; 3. rho for Kff; 4. len for Kff; 5. sigma2 (annealing it instead of optimizing it)
setopt.niter = niter; % number of iterations
setopt.ffmat = ffmat;


% Compute P-GPLVM with Laplace Approximation
result_la = pgplvm_la_new(yy,nf,setopt,xx);
% result_la = pgplvm_la_wb(yy,nf,setopt,xx);

end