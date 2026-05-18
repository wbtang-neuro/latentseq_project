function [f, df, hessian] = ln_poisson_model_cheeseboard(param,data,modelType,bin_nums)

X = data{1}; % subset of A
Y = data{2}; % number of spikes

% compute the firing rate
u = X * param;
rate = exp(u);

% roughness regularizer weight - note: these are tuned using the sum of f,
% and thus have decreasing influence with increasing amounts of data
b_pos = 8e0; b_hd = 5e1;b_gd1 = 5e1; b_gd2 = 5e1; b_gd3 = 5e1;

% start computing the Hessian
rX = bsxfun(@times,rate,X);       
hessian_glm = rX'*X;

%% find the P, H, S, or T parameters and compute their roughness penalties

% initialize parameter-relevant variables
J_pos = 0; J_pos_g = []; J_pos_h = []; 
J_hd = 0; J_hd_g = []; J_hd_h = [];  
J_gd1 = 0; J_gd1_g = []; J_gd1_h = [];
J_gd2 = 0; J_gd2_g = []; J_gd2_h = [];  
J_gd3 = 0; J_gd3_g = []; J_gd3_h = [];  

% find the parameters
numPos = bin_nums(1); numhd = bin_nums(2); numgd1 = bin_nums(3);numgd2 = bin_nums(4);numgd3 = bin_nums(5);
[param_pos,param_hd,param_gd1,param_gd2,param_gd3] = find_param(param,modelType,numPos,numhd,numgd1,numgd2,numgd3);

% compute the contribution for f, df, and the hessian
if ~isempty(param_pos)
    [J_pos,J_pos_g,J_pos_h] = rough_penalty_2d(param_pos,b_pos);
end

if ~isempty(param_hd)
    [J_hd,J_hd_g,J_hd_h] = rough_penalty_1d_circ(param_hd,b_hd);
end

if ~isempty(param_gd1)
    [J_gd1,J_gd1_g,J_gd1_h] = rough_penalty_1d_circ(param_gd1,b_gd1);
end

if ~isempty(param_gd2)
    [J_gd2,J_gd2_g,J_gd2_h] = rough_penalty_1d_circ(param_gd2,b_gd2);
end

if ~isempty(param_gd2)
    [J_gd3,J_gd3_g,J_gd3_h] = rough_penalty_1d_circ(param_gd3,b_gd3);
end

%% compute f, the gradient, and the hessian 

f = sum(rate-Y.*u) + J_pos + J_hd + J_gd1 + J_gd2 + J_gd3;
df = real(X' * (rate - Y) + [J_pos_g; J_hd_g; J_gd1_g;J_gd2_g;J_gd3_g]);
hessian = hessian_glm + blkdiag(J_pos_h,J_hd_h,J_gd1_h,J_gd2_h,J_gd3_h);

%% smoothing functions called in the above script
function [J,J_g,J_h] = rough_penalty_2d(param,beta)

    numParam = numel(param);
    D1 = spdiags(ones(sqrt(numParam),1)*[-1 1],0:1,sqrt(numParam)-1,sqrt(numParam));
    DD1 = D1'*D1;
    M1 = kron(eye(sqrt(numParam)),DD1); M2 = kron(DD1,eye(sqrt(numParam)));
    M = (M1 + M2);
    
    J = beta*0.5*param'*M*param;
    J_g = beta*M*param;
    J_h = beta*M;
    
function [J,J_g,J_h] = rough_penalty_1d_circ(param,beta)
    
    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    
    % to correct the smoothing across first and last bin
    DD1(1,:) = circshift(DD1(2,:),[0 -1]);
    DD1(end,:) = circshift(DD1(end-1,:),[0 1]);
    
    J = beta*0.5*param'*DD1*param;
    J_g = beta*DD1*param;
    J_h = beta*DD1;

function [J,J_g,J_h] = rough_penalty_1d(param,beta)

    numParam = numel(param);
    D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
    DD1 = D1'*D1;
    J = beta*0.5*param'*DD1*param;
    J_g = beta*DD1*param;
    J_h = beta*DD1;
   
%% function to find the right parameters given the model type
function [param_pos,param_hd,param_gd1,param_gd2,param_gd3] = find_param(param,modelType,numPos,numhd,numgd1,numgd2,numgd3)

param_pos = []; param_hd = []; param_gd1 = []; param_gd2 = []; param_gd3 = [];

if all(modelType == [1 1 0 0 0])
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numhd);
elseif all(modelType == [1 1 1 1 1]) 
    param_pos = param(1:numPos);
    param_hd = param(numPos+1:numPos+numhd);
    param_gd1 = param(numPos+numhd+1:numPos+numhd+numgd1);
    param_gd2 = param(numPos+numhd+numgd1+1:numPos+numhd+numgd1+numgd2);
    param_gd3 = param(numPos+numhd+numgd1+numgd2+1:numPos+numhd+numgd1+numgd2+numgd3);
end
    


