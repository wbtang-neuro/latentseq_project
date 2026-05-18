clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script fits LN models on single-cell demo data. 
% It identifies goal-direction and non-goal-direction cells via
% cross-validation.
% Reference: Hardcastle et al., Neuron, 2017
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%% define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

prefix = 'hp18_day16_20250420';
epoch = 3;
recnames = string([prefix,'-EP',num2str(epoch)]);
%% set parameters
savedata = 0; % save result?
speedthresh = 10; % apply speed threshold? in cm/s
bin = 0.2; % 200 ms
n_pos_bins = 20; % number of position bins
n_dir_bins = 12; % number of heading/goal-direction angle bins
numFolds = 10; %  10-fold cross-validation
% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
%% load data
load([firdir,prefix,'-EP',num2str(epoch),'.mat']) % load session data
epochtimes = [Dsession.t(1), Dsession.t(end)]; % epoch start and end time
hpnum = length(Dsession.units); % number of CA1 cells
%% bin running periods
RUN = vec2list((Dsession.speed > speedthresh),Dsession.t); % generate [start end] list of running epochs
taskIntervals = SplitIntervals(RUN,'pieceSize',bin); % 200 ms bins in running behavior
%% find goal-approaching segments
goaldir = zeros(size(Dsession.gd));
for goal = 1:3
    gd_current = Dsession.gd(:,goal);
    goaldir(:,goal) = atan2(sin(gd_current - Dsession.hd), cos(gd_current - Dsession.hd));
end
%% get [x,y,hd,gd] labels for each time bin
behavioral_labels_ep = [];
for i = 1:length(taskIntervals(:,1))
    timerange = taskIntervals(i,:);
    validids = find( Dsession.t >= timerange(1) & Dsession.t <=timerange(2) );
    if ~isempty(validids)
        behavioral_labels_ep = [behavioral_labels_ep;nanmean(Dsession.x(validids)),nanmean(Dsession.y(validids)),nanmean(Dsession.hd(validids)),...
            nanmean(goaldir(validids,:))];%[x,y,hd,gd]
    else
        behavioral_labels_ep = [behavioral_labels_ep;nan,nan,nan,nan,nan,nan];
    end
end
%% bin spikes, get spike matrix
spikes_CA1 = [];
for i = 1:hpnum
    timestamps = Dsession.units(i).spikeTimes;
    id = ones(size(timestamps))*i; % the id for that unit
    spikes_CA1 = [spikes_CA1; timestamps id]; % add [timestamps id] for this unit to the matrix
end
spikes_CA1 = sortrows(spikes_CA1,1); % sort spikes accoring to their time

id = spikes_CA1(:,2);

% Shift spike times to start at 0, and list bins unless explicitly provided
m = min([min(spikes_CA1(:,1)) min(taskIntervals(:))]);
spikes_CA1(:,1) = spikes_CA1(:,1) - m;
taskIntervals = taskIntervals - m;

% Create spike count matrix
nBins = size(taskIntervals,1);
if isempty(nBins), return; end
n = zeros(nBins,hpnum); % T x nCell
for unit = 1:hpnum
    temp = CountInIntervals(spikes_CA1(id==unit,1),taskIntervals);
    n(:,unit) = temp;
end
total_n = sum(n);
%% exclude invalid bins
validid = find(~isnan(behavioral_labels_ep(:,1)));

spike_matrix = n(validid,:);
behavioral_labels_ep = behavioral_labels_ep(validid,:); %[x,y,hd,gd]
%% exclude cells that don't have enough spike counts
validid = find(total_n > 100) ; % 100 spikes in total
spike_matrix = spike_matrix(:,validid);
%% fit LN model on single-unit response
cellIDs = 1:hpnum;
% compute position matrix
[posgrid, posxVec,posyVec] = feauture2D_map(behavioral_labels_ep(:,1:2), n_pos_bins);

% compute heading dir matrix
[hdgrid, hdVec] = feature_map(behavioral_labels_ep(:,3), n_dir_bins);

% compute goal-direction matrix (3 goals in total)
[gd1grid,gd1Vec] = feature_map(behavioral_labels_ep(:,4), n_dir_bins);
[gd2grid,gd2Vec] = feature_map(behavioral_labels_ep(:,5), n_dir_bins);
[gd3grid,gd3Vec] = feature_map(behavioral_labels_ep(:,6), n_dir_bins);

numModels = 2;% 1 full model, 2 goal-direction (gd) ablated model

testFit = cell(numModels,1);
trainFit = cell(numModels,1);
param = cell(numModels,1);
A = cell(numModels,1);
modelType = cell(numModels,1);

% Full model
A{1} = [ posgrid hdgrid gd1grid gd2grid gd2grid]; modelType{1} = [1 1 1 1 1];
% Ablated model
A{2} = [ posgrid hdgrid]; modelType{2} = [1 1 0 0 0]; % goal-direction (gd) ablated

% cell loop
for cellno = 1:numel(spike_matrix(1,:))
    fprintf('\t- EP-%d Fitting Cell %d of %d\n', epoch,cellno, numel(spike_matrix(1,:)));
    spiketrain = spike_matrix(:,cellno);
    fprintf('\t- Fitting the full model\n');
    [testFit_full,trainFit_full,param_full] = fit_model_cheeseboard(A{1},bin,spiketrain,filter,modelType{1},numFolds,[n_pos_bins^2,n_dir_bins,n_dir_bins,n_dir_bins,n_dir_bins]);
    LLH_values_full = testFit_full(:,3);%full model
    if any(~isnan(LLH_values_full))
        pval_baseline(1) = signrank(LLH_values_full,[],'tail','right');
        if pval_baseline(1) < 0.05 % only take the cell that show significant selectivity
            fprintf('\t- Fitting the second model\n');
            [testFit,trainFit,param]= fit_model_cheeseboard(A{2},bin,spiketrain,filter,modelType{2},numFolds,[n_pos_bins^2,n_dir_bins,n_dir_bins,n_dir_bins,n_dir_bins]);

            LLH_values = testFit(:,3);%second model
            pval_baseline = signrank(LLH_values,[],'tail','right');

            if pval_baseline < 0.05
                pvals = ones(1,2);
                [pvals(1),~] = signrank(testFit_full(:,1),testFit(:,1),'tail','right');
                [pvals(2),~] = signrank(testFit_full(:,3),testFit(:,3),'tail','right');

                p_max = max(pvals);
                p_min = min(pvals);
                fprintf('\t- Pval(max) = %f\n',p_max);
                fprintf('\t- Pval(min) = %f\n',p_min);
                if p_max < 0.05 % full model is sig. better
                    selected_model = 1; % full model
                elseif  p_min > 0.05 % full model is never better
                    selected_model = 2; % second model
                else % unclear
                    selected_model = NaN;
                end
                
            % re-set if selected model is not above baseline
            else
                selected_model = NaN;
            end

            selected_model_all(cellno) = selected_model;
            UID_select(cellno)  =  validid(cellno);
            model_info{cellno}.testFit_full = testFit_full;
            model_info{cellno}.trainFit_full = trainFit_full;
            model_info{cellno}.param_full = param_full;
            model_info{cellno}.testFit = testFit;
            model_info{cellno}.trainFit = trainFit;
            model_info{cellno}.param = param;
        else
            selected_model_all(cellno) = nan;
            UID_select(cellno)  =  validid(cellno);
            model_info{cellno} = [];
        end
    else
        selected_model_all(cellno) = nan;
        UID_select(cellno)  =  validid(cellno);
        model_info{cellno} = [];
    end
end
LNmodel.selected_model = selected_model_all;
LNmodel.UID = UID_select;  
LNmodel.model_info = model_info;
%% save result?
if savedata
    save(sprintf('%s%s.LNmodel_allVar-EP%02d.mat', resultdir, prefix,epoch),'LNmodel')
end




