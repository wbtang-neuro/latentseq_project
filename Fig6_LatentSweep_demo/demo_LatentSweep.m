clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% This script applies Bayesian decoding  on C1-C2 demo data using latent 
% states (from LMT). It also plots the example sweeps shown in Fig. 6.
% LMT references: Wu et al., 2017 and Luo et al., 2024
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%% define data
sweepsSetup
firdir = 'E:\Tang_latentseq_NN_2026\sample_data/'; % folder for data
resultdir = 'E:\Tang_latentseq_NN_2026\sample_results/'; % folder for saved results

prefix = 'hp18_day14_20250418';
epoch = 3; % C1 epoch
%% load sweep results
load([firdir,prefix,'-EP',num2str(epoch),'.mat']) % session info
dec = load([resultdir,prefix,'-EP',num2str(epoch),'.mat']); % hpc decoding result
%% get sweep info
chk = dec.chk;
sweeps = dec.sweeps;
sd = [sweeps.hpfPosDirection]';
hd = chk.hd;
egosd = circ_dist(sd, hd);
vswp = chk.speed> 15 & [sweeps.nvalid]'>3 & [sweeps.straight]'>.5;% remove noisy sweeps
egosd(~vswp) = nan;
vswp(find(abs(egosd) > pi/2)) = false; % remove reverse sweeps (minority)
%% get single-cell tuning curves to latent variables (states)
load(fullfile(resultdir,  [prefix,'.lmt.mat']));% load LMT results
gridbound = [];
for i=1:3
    gridbound = [gridbound; min(result_la.xxsamp(:,i)) max(result_la.xxsamp(:,i))];
end
xgrid = gen_grid(gridbound,30,nf);

fftc_init = get_tc(result_la.xxsamp,result_la.ffmat,xgrid,result_la.rhoff,result_la.lenff); % get single-cell latent tunings from the LMT model
tuning = exp(fftc_init) + eps;
epIDs1 = find(timevec >= Dsession.t(1) & timevec <= Dsession.t(end));
epIDs2 = find(timevec >= Dsession.t(end)+261.425); % remove the first few learning trials
%% get spike counts
units  = Dsession.units;
sc = reconstructSpikeCounts({units.spikeInds}, Dsession.nt);
for u = 1:numel(units)
    sc(:,u) =locsmooth(sc(:,u),1/0.01,0.05);
end
%% get decoding prob of latent states
prob_lmt = decodeBayes(sc, tuning, 0.01,1); % using LMT tunings

maxDist = (max(result_la.xxsamp(:,1))-min(result_la.xxsamp(:,1)))*0.12; % max jump distance
sweeps = sweeps(vswp);
decpos = gsmooth(dec.decpos, 0.8);
%% plot example sweeps
for s = [588, 640,888,932] % plot example sweeps shown in Fig. 6
% for s = 1:length(sweeps) % Alternatively, if you prefer, you can plot all the sweeps in this session.
    trng = [sweeps(s).tStart-5, sweeps(s).tStart+5];
    vt = restrictq(Dsession.t, trng);
    inds = find(vt);
    swpinds = sweeps(s).iSweep;
    prob_sweep = prob_lmt(swpinds,:)';
    npos = numel(xgrid(:,1));
    [maxprob, imx] = max(prob_sweep);
    posmax = xgrid(imx,:);
    pdiff = diff(posmax);
    % Copmute distance between points
    nextdist = sqrt(pdiff(:, 1).^2 + pdiff(:, 2).^2 + pdiff(:, 3).^2);
    isbad = find(nextdist > maxDist); % remove sweeps that have decoding exceed the max jump distance

    if isempty(isbad) % no bad points
        %% spatial space
        figure(1),clf
        vt = restrictq(timevec, trng);
        inds = find(vt);
        % plot latent states for C1 and C2
        scatter(result_la.xxsamp(epIDs1, 1), result_la.xxsamp(epIDs1, 2), 25,[0.6,0.6,0.6], 'filled', 'MarkerFaceAlpha', 0.5); %C1
        scatter(result_la.xxsamp(epIDs2, 1), result_la.xxsamp(epIDs2, 2), 25,[255,212,121]/255, 'filled', 'MarkerFaceAlpha', 0.5); %C2
        axis square
        % plot the trajectory
        plot(result_la.xxsamp(vt, 1), result_la.xxsamp(vt, 2), 'Color',[1,.5,.5]*.3, 'LineWidth',2);
        scatter(result_la.xxsamp(inds(1), 1), result_la.xxsamp(inds(1), 2), 100, 'sg','MarkerEdgeColor','k');
        % lmt prob
        for ii=1:3
            lmt_bins{ii} = linspace(gridbound(ii,1),gridbound(ii,2),50);
        end
        [Xpos,Ypos,Zpos] = meshgrid(1:50,1:50,1:50);
        index_xyz = [Xpos(:),Ypos(:),Zpos(:)];

        % plot the sweep
        xgrid_interp = gen_grid(gridbound,50,nf);
        fftc_init = get_tc(result_la.xxsamp,result_la.ffmat,xgrid_interp,result_la.rhoff,result_la.lenff);
        tuning_interp = exp(fftc_init) + eps;
        prob_lmt_seg = decodeBayes(sc(swpinds-10:swpinds+10,:), tuning_interp, 0.01,1); % using LMT tunings
        prob_lmt_seg = prob_lmt_seg(11:10+length(swpinds),:);

        problmt_gather = nan(50,50,length(swpinds));
        scalar = linspace(0,1,length(swpinds));
        for n  = 1:length(swpinds)
            prob_current = prob_lmt_seg(n,:)';
            prob_tmp = zeros(50,50);
            for i = 1:length(prob_current)
                currentInd = index_xyz(i,:);
                prob_tmp(currentInd(1),currentInd(2)) = prob_tmp(currentInd(1),currentInd(2)) + prob_current(i);
            end
            threshold = prctile(prob_tmp(:),99.5);
            prob_tmp(prob_tmp < threshold) = NaN;
            problmt_gather(:,:,n) = prob_tmp;
        end

        % Time-weighted sum
        weighted_sum = nansum(problmt_gather .* reshape(scalar, 1, 1, []), 3);
        % Total posterior sum
        total_sum = nansum(problmt_gather, 3);
        % Mean elapsed time (NaN where no valid data)
        mean_prob = nan(size(total_sum));
        valid_mask = total_sum > 0;
        mean_prob(valid_mask) = weighted_sum(valid_mask) ./ total_sum(valid_mask);

        h = imagesc(lmt_bins{1},lmt_bins{2},mean_prob);
        % Make NaNs transparent
        set(h, 'AlphaData', (~isnan(mean_prob))*0.7)
        coloarlabels = colormap(cool(n));

        colormap(coloarlabels)
        axis square
        % plot decoded positions
        plot(posmax(:, 1), posmax(:, 2),'cyan','linewidth',2)
        coloarlabels = colormap(cool(length(posmax(:, 1))));
        sswp = scatter(posmax(:, 1), posmax(:, 2), 50*ones(size(posmax(:, 1))),coloarlabels,'filled','MarkerEdgeColor','k');

        axis off
        disp('Press any key to plot the next sweep:')
        pause
    end
end