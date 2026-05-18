clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script applies Bayesian decoding on random foraging demo data. 
% It identifies theta sweeps and plots the polar histogram of the theta 
% sweep directions relative to heading, as well as the averaged sweeps.
% with the code base from Vollan et al., Nature, 2024
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%% define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

prefix = 'hp18_day22_20250426';
epoch = 2;
%% load data
recnames = [prefix,'-EP',num2str(epoch)];
load([firdir,recnames,'.mat']) % load Dsession info
%% theta sweep detection
%---detect theta sweeps
res_hc = runPosDecoding_wb_cheeseboard(save=0, rec=string(recnames),fld = resultdir,validbins_posxy = validbins_xy,validbins = validbins);
dec = res_hc;
%---Alternatively, to save computing time, load the existing result
% dec = load([resultdir,prefix,'-EP',num2str(epoch),'.mat']); % load hpc decoding result
%% get theta sweep info
chk = dec.chk;
sweeps = dec.sweeps;
sd = [sweeps.hpfPosDirection]';
hd = chk.hd;
egosd = circ_dist(sd, hd); % sweep direction relative to the heading direction

% remove noisy sweeps, using Vollan 2024 criteria 
vswp = chk.speed>S.minSpeed*100 & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5; 
egosd(~vswp) = nan;
%% theta sweep statistics
[prcAltern,pAltern] = computeAlternationPercent(egosd);
prcswp = 100*sum(vswp)./sum(chk.speed>S.minSpeed*100);
mulen = mean([dec.sweeps(vswp).length]);

% Print stats
fprintf("Percent sweep: %.2f\n", prcswp);
fprintf("Sweep length: %.2f\n", mulen*100);
fprintf("Alternation: %.2f\n", 100*prcAltern);
%% plot the polar histogram of sweep directions relative to heading 
figure("WindowStyle","normal")
tl = tiledlayout(1,1, 'TileIndexing','rowmajor');
plotEgoHistPolar(tl = tl, egodir=egosd,alpha=1);
set(gca, "FontSize", 10)
title("Heading")
%% plot averaged sweeps
figure("WindowStyle","normal")
clear cols;
fds = ["right", "left"];
cols.right = [1,0,0, .4];
cols.left = [0,0,1, .4];
clf

ax = nexttile; 
xline(0); yline(0);
xlim([-.08,.08]);
ylim([-.05,.15]);
axis square off
plot([.08,.08], [0,.1]+.05, 'k','LineWidth',2)
text(.09, 0.1, "10 cm")

title(ax, "Averaged sweeps")
D = Dsession;
[tmp] = plotAvgSweeps(D,dec);
for fd = fds
    plot(ax, tmp.(fd)(:, 1), tmp.(fd)(:, 2), '-',Color=cols.(fd))
    scatter(ax, tmp.(fd)(end, 1), tmp.(fd)(end, 2), 20,cols.(fd)(1:3), 'MarkerFaceAlpha',.7)
end



