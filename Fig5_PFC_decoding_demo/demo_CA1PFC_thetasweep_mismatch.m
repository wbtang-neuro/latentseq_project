clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script applies Bayesian decoding on PFC cell firing from demo data. 
% It compares the decoding mismatch between CA1 and PFC during goal vs. 
% lateral theta sweeps.
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%%  define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

animalprefix = 'HP18';
prefix = 'hp18_day16_20250420'; % demo session with good numbers of cells in both CA1 and PFC
epoch = 3;
%% load data
recnames = string([prefix,'-PFC-EP',num2str(epoch)]);
load([firdir,prefix,'-EP',num2str(epoch),'.mat']) % load hpc session info
dec = load([resultdir,prefix,'-EP',num2str(epoch),'.mat']); % load hpc decoding result
%% get CA1 sweep info
chk = dec.chk;
sweeps = dec.sweeps;
sd = [sweeps.hpfPosDirection]';
hd = chk.hd;
gd = chk.gd;
goalIDs  = Dsession.goalID(chk.iCen);
%% calculate CA1 goal and heading angles relative to sweeps 
egosd = circ_dist(sd, hd); %heading direction

goalsd = zeros(length(egosd),3); % goal direction
for i = 1:3
    goalsd(:,i) = circ_dist(sd, gd(:,i));
end
goaldir = nan(size(egosd));
for i = 1:length(goalIDs)
    if goalIDs(i)~=0
        goaldir(i) = goalsd(i,goalIDs(i));
    end
end
%% remove noisy sweeps
vswp = chk.speed>S.minSpeed*100 & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5;
egosd(~vswp) = nan;
goaldir(~vswp) = nan;
%% CA1 goal-directed sweeps vs lateral sweeps
g_sweepIDs = find(abs(goaldir) <= 10/180*pi);% less than 10 deg, goal sweeps
ng_sweepIDs = find(abs(goaldir) >= 40/180*pi);% larger than 40 deg, lateral sweeps
%% Bayesian decoding with PFC cells
res_pfc = runPosDecoding_wb_cheeseboard_gd(save=0, rec=recnames,fld = resultdir,validbins_posxy = validbins_xy,validbins = validbins);
decpos_PFC = gsmooth(res_pfc.decpos, 0.8);
decpos = gsmooth(dec.decpos, 0.8);
%% gather CA1-PFC decoding mismatch
% goal sweeps
gdecpos_gather = [];
for s = g_sweepIDs'
    swpIDs= sweeps(s).iStart:(sweeps(s).iStop);
    swp = decpos(swpIDs, :);
    swp_PFC = decpos_PFC(swpIDs, :);
    err = swp-swp_PFC;
    gdecpos_gather = [gdecpos_gather;hypot(err(:,1),err(:,2))]; % mismatch
end
% lateral sweeps
ngdecpos_gather = [];
for s = ng_sweepIDs'
    swpIDs= sweeps(s).iStart:(sweeps(s).iStop);
    swp = decpos(swpIDs, :);
    swp_PFC = decpos_PFC(swpIDs, :);
    err = swp-swp_PFC;
    ngdecpos_gather = [ngdecpos_gather;hypot(err(:,1),err(:,2))]; % mismatch
end
%% bar plot of the result
% mean
means = [mean(ngdecpos_gather), mean(gdecpos_gather)];
% SEM
sems = [std(ngdecpos_gather)/sqrt(length(ngdecpos_gather)), std(gdecpos_gather)/sqrt(length(gdecpos_gather))];

% plot bars
figure
b = bar(means);
hold on
b.FaceColor = 'flat';
b.CData(1,:) = [0 0 1];   % first bar (blue)
b.CData(2,:) = [1 0 0];         % second bar (red)

% add error bars
errorbar(1:2, means, sems, 'k', 'linestyle', 'none', 'LineWidth', 1.5)

set(gca,'XTick',1:2,'XTickLabel',{'Lateral sweeps','Goal sweeps'})
ylabel('CA1-PFC decoding mismatch (m)')
ylim([0.3,0.6])