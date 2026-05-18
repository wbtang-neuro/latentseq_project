clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script applies Bayesian decoding on goal-directed navigation demo data. 
% It identifies replay sweeps and plots the polar histogram of the replay 
% sweep directions relative to goal locations
% with the code base from Vollan et al., Nature, 2024
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%%  define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

animalprefix = 'HP18';
prefix = 'hp18_day16_20250420';
epoch = 3;
%% define paramters
speedthresh = 5; % cm/s, only take replay during immobility periods
%% load data
recnames = [prefix,'-EP',num2str(epoch)];
load([firdir,recnames,'.mat']) % load Dsession info
%% replay sweep detection
%---detect replay sweeps
res_hc = runPosDecoding_replay_cheeseboard_gd(save=0, rec=string(recnames),fld = resultdir,validbins_posxy = validbins_xy,validbins = validbins);
dec = res_hc;
% %---Alternatively, to save computing time, load the existing result
% dec = load([resultdir,prefix,'-EP',num2str(epoch),'_replay','.mat']); % load hpc decoding result
%% get replay sweep info
dec = res_hc;
chk = dec.chk;
sweeps = dec.sweeps;
sd = [sweeps.hpfPosDirection]';
hd = chk.hd;
pos = [Dsession.x, Dsession.y];
possm = dec.possm;
decpos = gsmooth(dec.decpos, 0.8);
%% get place-field decoding prob 
units  = Dsession.units;
sc = reconstructSpikeCounts({units.spikeInds}, Dsession.nt);
for u = 1:numel(units)
    sc(:,u) =locsmooth(sc(:,u),1/0.01,0.05);
end

nu = numel(units);
vbins = validbins;
for u = 1:nu
    tc = units(u).rmf;
    rmap = tc.z+eps;
    rmap(isnan(rmap))=eps;
    rmap = imgaussfilt(rmap,1);
    tuning_pf(:, u) = rmap(vbins);
end
prob = decodeBayes(sc, tuning_pf, 0.01,1); % using rate maps
%% remove noisy sweeps using Vollan's criteria
vswp = chk.speed<speedthresh & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5;
ss = find(vswp);% valid sweeps
replaysweep_idx = [];
for s = ss'
    trng = [sweeps(s).tStart-5, sweeps(s).tStart+5];
    vt = restrictq(Dsession.t, trng);
    inds = find(vt);
    [~,startID] = min(abs(Dsession.t - sweeps(s).tStart));
    swpinds = sweeps(s).iSweep;
    swpIDs= sweeps(s).iStart:(sweeps(s).iStop);
    swp = decpos(swpIDs, :);
    if length(swpinds) >2 && max(dec.spread(swpinds)) < 0.3 && ((swp(1,1) < 1.5) || (swp(end,1) < 1.5)) % pick replay that swept onto the maze
        replaysweep_idx = [replaysweep_idx;s];
    end
end
%% calculate replay sweep direction relative to the closest goal
goalxy = Dsession.goalxy;
delta_all = []; %angle, reset
for s = replaysweep_idx'
    swpIDs= sweeps(s).iStart:(sweeps(s).iStop);
    swp = decpos(swpIDs, :); % decoded location
    swp_start = swp(1,:);
    swp_stop = swp(end,:);
    if swp_start(1) < 1.5 || swp_stop(1) < 1.5 % sweep out the home box
        goaldist_start = swp_start - goalxy; 
        goaldist_start = hypot(goaldist_start(:, 1), goaldist_start(:, 2)); % distance to goal positions
        [mindist_start, mindist_goal_start] = min(goaldist_start);
        
        goaldist_stop = swp_stop - goalxy;
        goaldist_stop = hypot(goaldist_stop(:, 1), goaldist_stop(:, 2)); % distance to goal positions
        [mindist_stop, mindist_goal_stop] = min(goaldist_stop);
    
        % deal with forward and reverse replay
        if mindist_start < mindist_stop % start point is close to goal
            current_goal = mindist_goal_start;
            dpos = swp_stop - swp_start;
            sweepDirection = atan2(dpos(2), dpos(1));
            dpos_goal = swp_stop - goalxy(mindist_goal_start,:);
            sweepgd = atan2(dpos_goal(2), dpos_goal(1));
        else % end point is close to goal
            current_goal = mindist_goal_stop;
            dpos = swp_start - swp_stop;
            sweepDirection = atan2(dpos(2), dpos(1));
            dpos_goal = swp_start - goalxy(mindist_goal_stop,:);
            sweepgd = atan2(dpos_goal(2), dpos_goal(1));
        end
        delta = atan2(sin(sweepDirection - sweepgd), cos(sweepDirection - sweepgd));
        delta_all = [delta_all;delta];
    end
end
%% plot the angles
figure("WindowStyle","normal")
tl = tiledlayout(1,1, 'TileIndexing','rowmajor');
plotEgoHistPolar(tl = tl, egodir=delta_all,alpha=1);