clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script applies Bayesian decoding on goal-directed navigation demo data. 
% It identifies theta sweeps and plots the polar histogram of the theta 
% sweep directions relative to heading vs. to upcoming goal
% with the code base from Vollan et al., Nature, 2024
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%%  define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

prefix = 'hp18_day16_20250420';
epoch = 3;
%% load data
recnames = [prefix,'-EP',num2str(epoch)];
load([firdir,recnames,'.mat']) % load Dsession info
%% theta sweep detection
%---detect theta sweeps
res_hc = runPosDecoding_wb_cheeseboard_gd(save=0, rec=string(recnames),fld = resultdir,validbins_posxy = validbins_xy,validbins = validbins);
dec = res_hc;
%---Alternatively, to save computing time, load the existing result
% dec = load([resultdir,prefix,'-EP',num2str(epoch),'.mat']); % load hpc decoding result
%% get theta sweep info
chk = dec.chk;
sweeps = dec.sweeps;
sd = [sweeps.hpfPosDirection]';
hd = chk.hd;
gd = chk.gd;
goalIDs  = Dsession.goalID(chk.iCen);

% sweep direction relative to the heading direction
egosd = circ_dist(sd, hd);

% sweep direction relative to the upcoming goal
goalsd = zeros(length(egosd),3);
for i = 1:3
    goalsd(:,i) = circ_dist(sd, gd(:,i));
end
goaldir = nan(size(egosd));
for i = 1:length(goalIDs)
    if goalIDs(i)~=0
        goaldir(i) = goalsd(i,goalIDs(i));
    end
end

% shuffle the upcoming goal identity and calculate the angle again
goalIDs_shuf = zeros(size(goalIDs));
goalIDs_shuf(goalIDs == 1) = 3;
goalIDs_shuf(goalIDs == 2) = 1;
goalIDs_shuf(goalIDs == 3) = 2;
goaldir_shuf = nan(size(egosd));
for i = 1:length(goalIDs)
    if goalIDs(i)~=0
        goaldir_shuf(i) = goalsd(i,goalIDs_shuf(i));
    end
end

% remove noisy sweeps, using Vollan 2025 criteria 
vswp = chk.speed>S.minSpeed*100 & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5; 
egosd(~vswp) = nan;
goaldir(~vswp,1) = nan;
goaldir_shuf(~vswp,1) = nan;
%% theta sweep statistics
[prcAltern,pAltern] = computeAlternationPercent(egosd);
prcswp = 100*sum(vswp)./sum(chk.speed>S.minSpeed*100);
mulen = mean([dec.sweeps(vswp).length]);

% Print stats
fprintf("Percent sweep: %.2f\n", prcswp);
fprintf("Sweep length: %.2f\n", mulen*100);
fprintf("Alternation: %.2f\n", 100*prcAltern);
%% plot the polar histograms of sweep directions relative to heading, goal, and shuffled goal
figure("WindowStyle","normal")
tl = tiledlayout(1,3, 'TileIndexing','rowmajor');
plotEgoHistPolar(tl = tl, egodir=egosd,alpha=1);
title("Heading")
plotEgoHistPolar(tl = tl, egodir=goaldir,alpha=1);
title("To goal")
plotEgoHistPolar(tl = tl, egodir=goaldir_shuf,alpha=1);
title("To shuffled goal")



