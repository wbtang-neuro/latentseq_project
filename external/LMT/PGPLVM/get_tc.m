function postm = get_tc(xx,ffmat,xgrid,rhoff,lenff)

ffTYPE = 2;
covfun = covariance_fun(rhoff,lenff,ffTYPE); % get the covariance function

[K0,~] = covfun(xx,xx);
% K0 = covfun(xx,xx);
% add sigma to diagonal of covariance matrix (adds white
% noise?)
sdiag = K0(1,1)*1*eye(size(K0));
K0 = K0 + sdiag;
invK = pdinv(K0);
% K = K0+1e-4*eye(size(xx,1));
% K = K0;
K1 = covfun(xgrid,xx);

% postm = K1*(K\ffmat);
postm = K1*invK*ffmat;

