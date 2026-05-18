function [testFit,trainFit,param_mean] = fit_model_cheeseboard(A,dt,spiketrain,filter,modelType,numFolds,bin_nums)

%% Description
% This code will section the data into 10 different portions. Each portion
% is drawn from across the entire recording session. It will then
% fit the model to 9 sections, and test the model performance on the
% remaining section. This procedure will be repeated 10 times, with all
% possible unique testing sections. The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.


%% Initialize matrices and section the data for k-fold cross-validation

[~,numCol] = size(A);
% sections = numFolds;
sections = numFolds*5;


% divide the data up into 5*num_folds pieces
edges = round(linspace(1,numel(spiketrain)+1,sections+1));

% initialize matrices
testFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of test data
trainFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of train data
paramMat = nan(numFolds,numCol);

%% perform k-fold cross validation
count = 0;
% seeds = [1,2,3,4,5];
% for run = 1:5
    rng(123);
    p = randperm(numel(spiketrain));
    for k = 1:numFolds
        count = count +1;
        % fprintf('\t\t- Cross validation fold %d of %d\n', k, numFolds);
        
        % get test data from edges - each test data chunk comes from entire session
        % test_ind  = edges(k):edges(k+1)-1;
        % % get test data from edges - each test data chunk comes from entire session
        test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
            edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
            edges(k+4*numFolds):edges(k+4*numFolds+1)-1]   ;
        % 
        % test_ind  = p(test_ind);
        test_spikes = spiketrain(test_ind); %test spiking
        [~,sortids] = sort(test_ind);
        test_spikes = test_spikes(sortids);
        smooth_spikes_test = conv(test_spikes,filter,'same'); %returns vector same size as original
        smooth_fr_test = smooth_spikes_test./dt;
        test_A = A(test_ind,:);
        test_A = test_A(sortids,:);
        
        % training data
        train_ind = setdiff(1:numel(spiketrain),test_ind);
        train_ind  = p(train_ind);
        train_spikes = spiketrain(train_ind); %test spiking
        [~,sortids] = sort(train_ind);
        train_spikes = train_spikes(sortids);
        smooth_spikes_train = conv(train_spikes,filter,'same'); %returns vector same size as original
        smooth_fr_train = smooth_spikes_train./dt;
        train_A = A(train_ind,:);
        train_A = train_A(sortids,:);
        
        
        opts = optimset('Gradobj','on','Hessian','on','Display','off','Algorithm','trust-region');
        
        data{1} = train_A; data{2} = train_spikes;
        if k == 1
            init_param = 1e-3*randn(numCol, 1);
        else
            init_param = param;
        end
        [param] = fminunc(@(param) ln_poisson_model_cheeseboard(param,data,modelType,bin_nums),init_param,opts);
        
        %%%%%%%%%%%%% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%
        % compute the firing rate
        fr_hat_test = exp(test_A * param)/dt;
        smooth_fr_hat_test = fr_hat_test; %returns vector same size as original
        
        % compare between test fr and model fr
        sse = sum((smooth_fr_hat_test-smooth_fr_test).^2);
        sst = sum((smooth_fr_test-mean(smooth_fr_test)).^2);
        varExplain_test = 1-(sse/sst);
        
        % compute correlation
        correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,'type','Pearson');
        
        % compute llh increase from "mean firing rate model" - NO SMOOTHING
        r = exp(test_A * param); n = test_spikes; meanFR_test = nanmean(test_spikes); 
        
        log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
        log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
        log_llh_test = (-log_llh_test_model + log_llh_test_mean);
        log_llh_test = log_llh_test/log(2); % convert from nats to bits
        
        % compute MSE
        mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
        
        % fill in all the relevant values for the test fit cases
        testFit(count,:) = [varExplain_test correlation_test log_llh_test mse_test sum(n) numel(test_ind) ];
        
        %%%%%%%%%%%%% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%
        % compute the firing rate
        fr_hat_train = exp(train_A * param)/dt;
        smooth_fr_hat_train = fr_hat_train; %returns vector same size as original
        
        % compare between test fr and model fr
        sse = sum((smooth_fr_hat_train-smooth_fr_train).^2);
        sst = sum((smooth_fr_train-mean(smooth_fr_train)).^2);
        varExplain_train = 1-(sse/sst);
        
        % compute correlation
        correlation_train = corr(smooth_fr_train,smooth_fr_hat_train,'type','Pearson');
        
        % compute log-likelihood
        r_train = exp(train_A * param); n_train = train_spikes; meanFR_train = nanmean(train_spikes);   
        log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
        log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
        log_llh_train = (-log_llh_train_model + log_llh_train_mean);
        log_llh_train = log_llh_train/log(2); % convert from nats to bits
        
        % compute MSE
        mse_train = nanmean((smooth_fr_hat_train-smooth_fr_train).^2);
        
        trainFit(count,:) = [varExplain_train correlation_train log_llh_train mse_train sum(n_train) numel(train_ind)];
    
        % save the parameters
        paramMat(count,:) = param;
    
    end
% end
param_mean = nanmean(paramMat);

return
