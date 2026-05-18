clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script calculates theta modulation of goal-direction vs. 
% non-goal-direction cell response during goal sweeps from demo data. 
% see demo_LNmodel_GoalCells.m for goal-direction cell detection
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%% define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

prefix = 'hp18_day16_20250420';
epoch = 3;

recnames = string([prefix,'-EP',num2str(epoch)]);
%% define parameters
binsize = 1; % ms
nstd = 15; % std of the gaussian smooth window
g1 = gaussian(nstd, 3*nstd+1); % gaussian smooth window
%% load data
load([firdir,prefix,'-EP',num2str(epoch),'.mat'])
dec = load([resultdir,prefix,'-EP',num2str(epoch),'.mat']); % load hpc decoding result
hpnum = length(Dsession.units);
%% load the LN model info
cellID = 1:hpnum;
load(sprintf('%s%s.LNmodel_allVar-EP%02d.mat', resultdir, prefix,epoch))
%% plot firing rates of goal-direction vs. non-goal-direction cells over theta phases during lateral and goal sweeps
figure("WindowStyle","normal",'Position',[580,680,415,240])
tl = tiledlayout(1,1);
for model = 1:2
    validid = find(LNmodel.selected_model == model); % 1 = goal-direction cells; 2 = non-goal-directional cells
    CellIDs = cellID(LNmodel.UID(validid));
    %% firing rate trace for (non-)goal-direction cells
    totaltime = Dsession.t(end)- Dsession.t(1);
    MUAtime = Dsession.t(1)*1000:binsize:Dsession.t(end)*1000;
    % Get the spike times
    MUA_spikes_mat = [];
    cellcount = 0;
    for cell = CellIDs
        if ~isempty(Dsession.units(cell).spikeTimes)
            spikeu = Dsession.units(cell).spikeTimes*1000;  % in ms;
            valid_id = find(spikeu >=  Dsession.t(1)*1000 & spikeu<=  Dsession.t(end)*1000);
            spikeu = spikeu(valid_id);
            if length(spikeu) > 10 % at least 10 spikes
                % Get firing rate of neuron
                histspks_all = histc(spikeu,MUAtime);
                MUA_spikes_mat = [MUA_spikes_mat,histspks_all];
            end
        end
    end
    MUA_spikes = sum(MUA_spikes_mat,2);

    % gaussian smooth
    MUA_spikes = smoothvect(MUA_spikes, g1);
    MUA_PYR = [MUAtime'/1000,MUA_spikes];
    %% calculate goal and heading angles relative to sweeps
    chk = dec.chk;
    sweeps = dec.sweeps;
    sd = [sweeps.hpfPosDirection]';
    hd = chk.hd;
    hdsm = chk.hdsm;
    gd = chk.gd;
    goalIDs  = Dsession.goalID(chk.iCen);

    egosd = circ_dist(sd, hd); %heading direction

    goalsd = zeros(length(egosd),3); % goal direction
    for i = 1:3
        goalsd(:,i) = circ_dist(sd, gd(:,i));
    end

    goaldir = nan(length(egosd),2);
    for i = 1:length(goalIDs)
        if goalIDs(i)~=0
            goaldir(i,1:2) = [goalsd(i,goalIDs(i)),goalIDs(i)];
        end
    end
    %% remove noisy sweeps
    vswp = chk.speed>S.minSpeed*100 & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5;
    egosd(~vswp) = nan;
    goaldir(~vswp) = nan;
    %% detect goal-directed sweeps 
    g_sweepIDs = find(abs(goaldir(:,1)) <= 10/180*pi);% less than 10 deg to goal, goal sweeps
    %% get theta modulation of MUA firing rate (goal sweeps)
    theta_phase = [];
    MUA_trace = [];
    for s = g_sweepIDs'
        tStart =  sweeps(s).tStart;
        thetaID = find(Dsession.thetaChunks.tStart(1:end-1) <= tStart & Dsession.thetaChunks.tStart(2:end) > tStart);
        timeID = find(Dsession.t >= Dsession.thetaChunks.tStart(thetaID) & Dsession.t < Dsession.thetaChunks.tStart(thetaID+1));
        theta_phase = [theta_phase;Dsession.t(timeID),Dsession.theta(timeID)];

        timeID = find(MUA_PYR(:,1) >= Dsession.thetaChunks.tStart(thetaID) & MUA_PYR(:,1) < Dsession.thetaChunks.tStart(thetaID+1));
        MUA_trace = [MUA_trace;MUA_PYR(timeID,:)];
    end
    theta_phase(:,2) = mod(theta_phase(:,2), 2*pi);
    theta_phase(:,2) = theta_phase(:,2)./(2*pi);
    %% get circular rate map (goal sweeps)
    MUA_phase_map = Map(theta_phase,MUA_trace,'smooth',1.5,'minTime',0,...
        'nBins',51,'maxGap',0.1,'mode','discard','maxDistance',5,'type','cl');
    %% plot results
    if model == 1
        pax = polaraxes(tl);
        pax.ThetaZeroLocation = "left";
        pax.ThetaTick = [0, 90, 180, 270];
        pax.ThetaTickLabel = ["","behind","", "ahead"];
        theta = (MUA_phase_map.x*2*pi);
        hold on
        aa = MUA_phase_map.z';
        polarplot(theta,aa/0.01,'r-','LineWidth',2);

        rmax = max(aa/0.01);
    else
        aa = MUA_phase_map.z';
        rmin = min(aa/0.01);
        polarplot(theta,aa/0.01,'k-','LineWidth',2);
    end

end
rlim([rmin rmax]); % Sets the r-axis limits
legend('GD cells','NGD cells')
