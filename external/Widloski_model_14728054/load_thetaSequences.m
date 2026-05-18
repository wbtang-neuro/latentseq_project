% clear all

phaseOffset_vec = pi;%linspace(0,2*pi,20);
for kk = 1:length(phaseOffset_vec)
kk
rat = 3;

load Experiment_Information
load Analysis_Information
load clusters
load behavior
load spikeDensity
% load binDecoding_2, decoder_binDecoding_slow = decoder_binDecoding;
load binDecoding_15, decoder_binDecoding_slow = decoder_binDecoding;
load binDecoding_05, 
load mazeTimes

cd LFP_merged.LFP/
if rat == 1, load LFP11.mat, end
if rat == 2, load LFP4.mat, end
if rat == 3, load LFP4.mat, end
cd ..

%LFP
    LFP_full = LFP_Data; clear LFP_data

%positions
    times = load_timeBins(Times_day,shiftSizeDecoding,windowSizeDecoding_replay);
    positions_full = load_positions_full(Run_Times,Sleep_Times,mean(times,2),positions);
    positions_full(isnan(positions_full(:,4)),4) = 0;
    
%spike density
    spikeDensity_full = spikeDensity;
        
sessionNum_decoder = 2;
sessionNum = sessionNum_decoder;

    lickingTimes = load_fieldVec([session(sessionNum).trials],'lickingTimes',2);

    %bin decoding
        speedRangeThr = [10 inf];
        sequence_posteriorSpreadThr = 15;
        sequence_jumpThr = 20;
        [x_NaN,x,timeBins,posteriorSpread,spikeDensity,positions] = compute_filtering_binDecoding(decoder_binDecoding(sessionNum_decoder),sequence_posteriorSpreadThr,speedRangeThr,sequence_spikeDensityThr,positions_full,spikeDensity_full,spikeDensityStepSize,sequence_jumpThr,Times_day,binSize);
        [x_NaN_slow,x_slow,timeBins_slow,posteriorSpread_slow,spikeDensity_slow,positions_slow] = compute_filtering_binDecoding(decoder_binDecoding_slow(sessionNum_decoder),sequence_posteriorSpreadThr,speedRangeThr,sequence_spikeDensityThr,positions_full,spikeDensity_full,spikeDensityStepSize,sequence_jumpThr,Times_day,binSize);
            x_NaN_slow = compute_dataInterpolation(x_NaN_slow,x_NaN(:,1));
            x_slow = compute_dataInterpolation(x_slow,x_NaN(:,1));

            % v_slow = [diff(x_NaN_slow(:,2)), diff(x_NaN_slow(:,3))]/mode(diff(x_NaN_slow(:,1))); v_slow = [v_slow(1,:); v_slow];
            % speed = sqrt(v_slow(:,1).^2 + v_slow(:,2).^2);
            % angle = atan2(v_slow(:,2),v_slow(:,1));
            % x_NaN_slow = [x_NaN_slow speed angle];

    %rat position
        positions_scaled = compute_locsToBins(positions,x_edges,y_edges);

    %LFP
    key = 0;
    % phaseOffset = 0; %!!!!!!!!!!!!!!!
    phaseOffset = phaseOffset_vec(kk);
    while key == 0
        %raw LFP (time, raw LFP, smoothed LFP)
        LFP = compute_dataTemporalConcatenation(LFP_full,load_timeBounds(Run_Times{sessionNum})); [~,ind_unique,~] = unique(LFP(:,1)); LFP = LFP(ind_unique,:);
        LFP_smoothed = conv(LFP(:,2),setUp_gaussFilt([1 1000],0.1*LFPSampRate),'same');
        LFP = [LFP(:,1:2) LFP_smoothed];

        %theta oscillations (time, theta LFP, theta phase, theta power)
        [LFP_filt_thetaPhase,LFP_filt_thetaPower,LFP_filt_theta] = compute_filteredLFP(ThetaFreqRange,LFP(:,2),LFPSampRate); 
        LFP_filt_thetaPower_smoothed = conv(LFP_filt_thetaPower,setUp_gaussFilt([1 1000],100),'same');
        LFP_theta = [LFP(:,1) LFP_filt_theta LFP_filt_thetaPhase LFP_filt_thetaPower_smoothed];
        LFP_theta(:,3) = mod(LFP_theta(:,3)-phaseOffset,2*pi);
        LFP_theta = compute_dataInterpolation(LFP_theta,x(:,1),[]);

            %from LFP_filtered
            % LFP_filtered_sub = compute_dataInterpolation(LFP_filtered,x(:,1),[]);
            % LFP_theta = [LFP_filtered_sub(:,1) LFP_filtered_sub(:,3) LFP_filtered_sub(:,8)];

        %Find  theta demarcations and remove theta cycles associated with slow movement
            [~,theta_locs] = findpeaks(LFP_theta(:,3),'minpeakheight',2*pi-0.4);
            theta_timeBounds = [LFP_theta(theta_locs(1:end-1)+1,1) LFP_theta(theta_locs(2:end),1)];
            theta_time = theta_timeBounds(:,1);

        theta_rat = compute_dataInterpolation(positions_scaled,theta_time); 
        theta_time(theta_rat(:,4)<speedThr,:) = [];
        theta_timeBounds(theta_rat(:,4)<speedThr,:) = [];
        theta_rat(theta_rat(:,4)<speedThr,:) = [];
            x_NaN_removed = x_NaN; x_NaN_removed(isnan(x_NaN(:,2)),:) = [];
            x_NaN_slow_removed = x_NaN_slow; x_NaN_slow_removed(isnan(x_NaN_slow(:,2)),:) = [];
        x_NaN_theta = compute_dataInterpolation(x_NaN_removed,theta_time);
        x_NaN_slow_theta = compute_dataInterpolation(x_NaN_slow_removed,theta_time);
        numThetaCycles = length(theta_time);

            % histogram(theta_timeBounds(:,2)-theta_timeBounds(:,1),linspace(0,0.25,20),'FaceColor','k','EdgeColor','none','Normalization','probability')
            % axis square, xlabel('theta period (sec)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
            % title(compute_round(nanmean(theta_timeBounds(:,2)-theta_timeBounds(:,1)),100)), axis tight

            %Sample uniformly within the theta cycle 
                LFP_theta_phase = [theta_time,2*pi*(0:length(theta_time)-1)'];
                numBins = 64; 
                phaseShiftSize = 2*pi/numBins;
                phase = (0:phaseShiftSize:max(LFP_theta_phase(:,2))-phaseShiftSize)';
                times_theta_unifPhaseSamp = compute_dataInterpolation([LFP_theta_phase(:,2) LFP_theta_phase(:,1)],phase,[]); times_theta_unifPhaseSamp = times_theta_unifPhaseSamp(:,2);

            %resample uniformly according to theta phase
                spikeDensity_phase = compute_dataInterpolation(spikeDensity,times_theta_unifPhaseSamp,[]);
                spikeDensity_theta = reshape(spikeDensity_phase(:,end),numBins,size(spikeDensity_phase,1)/numBins)';
                x_NaN_phase = compute_dataInterpolation(x_NaN,times_theta_unifPhaseSamp,[]);
                x_NaN_slow_phase = compute_dataInterpolation(x_NaN_slow,times_theta_unifPhaseSamp,[]);
                positions_phase = compute_dataInterpolation(positions_scaled,times_theta_unifPhaseSamp,5);

                distanceToRat_phase = sqrt((x_NaN_phase(:,2)-positions_phase(:,2)).^2 + (x_NaN_phase(:,3)-positions_phase(:,3)).^2);
                distanceToRat_theta = reshape(distanceToRat_phase(:,end),numBins,size(distanceToRat_phase,1)/numBins)';

                ind_thetaCycles = find(theta_rat(:,4)>10);
                ind_thetaCycles(ind_thetaCycles>size(distanceToRat_phase,1)/numBins) = [];

                key = 1; %!!!!!!!!!!!!!!!!!!!

                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % if phaseOffset~=0 && key == 0
                    % subplot(121)
                    % plot(180/pi*linspace(0,2*pi,64),nanmean(distanceToRat_theta(ind_thetaCycles,:)),'k','linewidth',4)
                    % xlabel('phase (deg)'), ylabel({'distance from the rat'}), axis square, set(gca,'fontsize',14,'xtick',[0 180 360]), set(gcf,'Color','w'), axis tight
                    % 
                    % subplot(122)
                    % plot(180/pi*linspace(0,2*pi,64),nanmean(spikeDensity_theta(ind_thetaCycles,:)),'k','linewidth',4)
                    % xlabel('phase (deg)'), ylabel({'spike density','(z-scored)'}), axis square, set(gca,'fontsize',14,'xtick',[0 180 360]), set(gcf,'Color','w'), axis tight
                    % 
                    % key = 1;
                    % keyboard
                % end
                % 
                % %compute phase offset
                % dd = nanmean(distanceToRat_theta(ind_thetaCycles,:));
                % xx = linspace(0,2*pi,64);
                % [pk,loc] = find(dd==max(dd));
                % phaseOffset = xx(loc);

    end

    %theta sequence properties
        x_NaN_snippet_origin_rot_full = nan(numThetaCycles,400,2);
        theta_seqDir = nan(numThetaCycles,3);
        theta_seqDist = nan(numThetaCycles,1);
        theta_seqDistFromRat = nan(numThetaCycles,1);
        theta_seqDistFromRep = nan(numThetaCycles,1);
        for i = 2:numThetaCycles
            % i/numThetaCycles

            x_NaN_snippet = compute_dataTemporalConcatenation(x_NaN,theta_timeBounds(i,:));
            LFP_theta_snippet = compute_dataTemporalConcatenation(LFP_theta,theta_timeBounds(i,:));

            %angle relative to rat heading
                x_NaN_snippet_cut = x_NaN_snippet;
                    ind_nan = find(LFP_theta_snippet(:,3)<pi | isnan(x_NaN_snippet_cut(:,2)) | isnan(x_NaN_snippet_cut(:,3)));
                    % ind_nan = find(isnan(x_NaN_snippet_cut(:,2)) | isnan(x_NaN_snippet_cut(:,3)));
                    x_NaN_snippet_cut(ind_nan,:) = [];
                distanceToRat_snippet = sqrt((x_NaN_snippet_cut(:,2)-theta_rat(i,2)).^2 + (x_NaN_snippet_cut(:,3)-theta_rat(i,3)).^2 );
                distanceToRep_snippet = sqrt((x_NaN_snippet_cut(:,2)-x_NaN_slow_theta(i,2)).^2 + (x_NaN_snippet_cut(:,3)-x_NaN_slow_theta(i,3)).^2 );
                if isempty(x_NaN_snippet_cut) || isnan(x_NaN_slow_theta(i,2)), continue, end

                % %method 1
                %     delta = diff(x_NaN_snippet_cut(:,2:3),[],1);
                %         delta(prod(delta,2)==0,:) = [];
                %     theta_seqDir(i,1) = circ_mean(atan2(delta(:,2),delta(:,1)));
                %     theta_seqDir(i,2) = circ_dist(theta_rat(i,5),theta_seqDir(i,1));
                %     theta_seqDir(i,3) = circ_dist(theta_seqDir(i-1,1),theta_seqDir(i,1));
                % 
                % %method 2
                %     ind_max = find(distanceToRat_snippet==max(distanceToRat_snippet)); ind_max = ind_max(1);
                %     ind_min = find(distanceToRat_snippet==min(distanceToRat_snippet)); ind_min = ind_min(1);
                %     delta = diff([theta_rat(i,2:3); x_NaN_snippet_cut(ind_max,2:3)],[],1);
                %     theta_seqDir(i,4) = atan2(delta(:,2),delta(:,1));
                %     theta_seqDir(i,5) = circ_dist(theta_rat(i,5),theta_seqDir(i,4));
                %     theta_seqDir(i,6) = circ_dist(theta_seqDir(i-1,4),theta_seqDir(i,4));

                %method 3
                    ind_max = find(distanceToRep_snippet==max(distanceToRep_snippet)); ind_max = ind_max(1);
                    delta = diff([x_NaN_slow_theta(i,2:3); x_NaN_snippet_cut(ind_max,2:3)],[],1);
                    theta_seqDir(i,1) = atan2(delta(:,2),delta(:,1));
                    theta_seqDir(i,2) = circ_dist(theta_rat(i,5),theta_seqDir(i,1));
                    theta_seqDir(i,3) = circ_dist(theta_seqDir(i-1,1),theta_seqDir(i,1));

            %sequence distance
                theta_seqDist(i) = compute_sequenceDistance(x_NaN_snippet(:,2:3));
                theta_seqDistFromRat(i) = max(distanceToRat_snippet);
                theta_seqDistFromRep(i) = max(distanceToRep_snippet);

            %rotate relative to rat heading
                % x_NaN_snippet_origin = x_NaN_snippet(:,2:3) - repmat(theta_rat(i,2:3),size(x_NaN_snippet,1),1);
                x_NaN_snippet_origin = x_NaN_snippet(:,2:3) - repmat(x_NaN_slow_theta(i,2:3),size(x_NaN_snippet,1),1);
                M = [cos(-theta_rat(i,5)) -sin(-theta_rat(i,5)); sin(-theta_rat(i,5)) cos(-theta_rat(i,5))];
                x_NaN_snippet_origin_rot = (M*x_NaN_snippet_origin')';
                x_NaN_snippet_origin_rot_full(i,:,:) = compute_vecBuffered(x_NaN_snippet_origin_rot,400);

        end

        % save('R1D4S3.mat','theta_rat','theta_seqDir','theta_seqDist','x_NaN_snippet_origin_rot_full')

        % keyboard

        % load R1D3S3.mat
        % load R1D4S3.mat

            speedThr_thetaCycle = 20;
            ind_thetaCycles = find(theta_rat(:,4)>speedThr_thetaCycle);

            %temporal regularity
                subplot(131)
                delT = squareform(pdist(theta_rat(ind_thetaCycles,1))); delT = delT(:);
                delX = squareform(pdist(sign(theta_seqDir(ind_thetaCycles,2)))); delX = delX(:)/2;
                edges_time = linspace(0,2,60); centers_time = edges_time(1:end-1) + mean(diff(edges_time))/2;
                data = nan(3,length(edges_time)-1);
                for i = 1:length(edges_time)-1
                    ind = find(delT>edges_time(i) & delT<=edges_time(i+1));
                    data(1,i) = nanmean(delX(ind));
                    data(2,i) = nanstd(delX(ind))./sqrt(length(ind));
                    data(3,i) = length(ind);
                end
                shadedErrorBar(centers_time,100*data(1,:),100*data(2,:),'lineprops',{'k','linewidth',2})
                xlabel('time between events (sec)'), ylabel('% alternation'), set(gca,'fontsize',14)
                hline(50,'k--'), axis square, box on, axis tight

            %angular regularity
                subplot(132)
                edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
                h1 = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges_angles,'normalization','probability');
                h2 = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges_angles,'normalization','probability');
                h1 = smooth(h1,3); h2 = smooth(h2,3);
                yyaxis left, plot(180/pi*centers_angles,h2,'k','linewidth',4)
                    ylabel('fraction')
                yyaxis right, plot(180/pi*centers_angles,h1,'r','linewidth',4)
                xlabel('angle (deg)'), axis square, set(gca,'fontsize',14,'xtick',[0 90 180]), set(gcf,'color','w'), xlim(180*[0 1])
                ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; 
                legend('prev. sweep','rat heading')
                ax = gca;
                ax.YAxis(1).Color = 'k';
                ax.YAxis(2).Color = 'r';
                set(gca,'xtick',0:60:180)

            %sequence distance
                subplot(133)
                edges_distance = linspace(0,40,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
                histogram(theta_seqDist(ind_thetaCycles),edges_distance,'facecolor','k','edgecolor','none','Normalization','probability');
                axis square, xlabel('distance (cm)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
   
            keyboard
            
        %relative angle vs rat speed
            [edges_speed,centers_speed] = load_timeBins([20 60],2,10);
            a = cool(length(edges_speed));
            H1 = []; H2 = [];
            for i = 1:length(centers_speed)
                ind_thetaCycles = find(theta_rat(:,4)>edges_speed(i,1) & theta_rat(:,4)<=edges_speed(i,2));
                h1 = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges_angles,'normalization','probability');
                h2 = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges_angles,'normalization','probability');
                h1 = smooth(h1,3)'; h2 = smooth(h2,3)';
                H1 = [H1; h1/max(h1)]; H2 = [H2; h2/max(h2)];
            end
            % for i = 1:2
            %     subplot(1,2,i)
            %     if i==1, h = H1; else, h = H2; end
            %     h = pcolor(180/pi*centers_angles,centers_speed,h); set(h,'edgecolor','none')
            %     cb = colorbar(); ylabel(cb,'fraction','FontSize',14,'Rotation',270)
            %     xlabel('angle (deg)'), ylabel('speed')
            %     axis square, set(gca,'fontsize',14), set(gcf,'color','w')
            % end

            subplot(4,5,kk)
            h = pcolor(180/pi*centers_angles,centers_speed,H2); set(h,'edgecolor','none')
            cb = colorbar(); ylabel(cb,'fraction','FontSize',14,'Rotation',270)
            xlabel('angle (deg)'), ylabel('speed')
            axis square, set(gca,'fontsize',14), set(gcf,'color','w')
            drawnow

end


    % %movement bouts
    %     positions_sub = positions_scaled;
    %     positions_sub(positions_scaled(:,4)<20,:) = nan;
    %     [boundaries,lengths] = compute_allSequences_NaNseparated(positions_sub(:,1));     
    %     times_movementBouts = [positions_scaled(boundaries(:,1),1) positions_scaled(boundaries(:,2),1)];
    %     times_movementBouts(times_movementBouts(:,2)-times_movementBouts(:,1)<1,:) = [];
    % 
    %     for i = 1:size(times_movementBouts,1)
    % 
    %         Times_snippet = times_movementBouts(i,:);
    %         Times_snippet = [Run_Times{sessionNum_decoder}(1,1)+13*60+58.24 Run_Times{sessionNum_decoder}(1,1)+13*60+58.24+1.22];
    %         % Times_snippet = [Run_Times{sessionNum_decoder}(1,1)+2*60+47.5 Run_Times{sessionNum_decoder}(1,1)+2*60+47.5+1.2];
    %         % Times_snippet = [Run_Times{sessionNum_decoder}(1,1)+4*60+26.6 Run_Times{sessionNum_decoder}(1,1)+4*60+26.6+2];
    %         % Times_snippet = [Run_Times{sessionNum_decoder}(1,1)+7*60+7.4 Run_Times{sessionNum_decoder}(1,1)+7*60+7.4+1];
    %         % Times_snippet = [Run_Times{sessionNum_decoder}(1,1)+10*60+50.6 Run_Times{sessionNum_decoder}(1,1)+10*60+50.6+1.5];
    % 
    %         positions_snippet = compute_dataTemporalConcatenation(positions_scaled,Times_snippet);
    %         spikeDensity_snippet = compute_dataTemporalConcatenation(spikeDensity,Times_snippet);
    %         LFP_snippet = compute_dataTemporalConcatenation(LFP,Times_snippet);
    %         x_NaN_snippet = compute_dataTemporalConcatenation(x_NaN,Times_snippet);
    %         x_NaN_slow_snippet = compute_dataTemporalConcatenation(x_NaN_slow,Times_snippet);
    %         LFP_theta_snippet = compute_dataTemporalConcatenation(LFP_theta,Times_snippet);
    %         [theta_time_snippet,theta_time_snippet_nan] = compute_dataTemporalConcatenation([theta_time theta_time],Times_snippet);
    %             ind_thetaCycles = find(~isnan(theta_time_snippet_nan(:,2)));
    %             theta_time_snippet = theta_time_snippet(:,2); 
    %         x_NaN_theta_snippet = compute_dataTemporalConcatenation(x_NaN_theta,Times_snippet);
    %         x_NaN_slow_theta_snippet = compute_dataTemporalConcatenation(x_NaN_slow_theta,Times_snippet);
    %         x_NaN_snippet_origin_rot_full_sub = [reshape(x_NaN_snippet_origin_rot_full(ind_thetaCycles,:,1),[],1) reshape(x_NaN_snippet_origin_rot_full(ind_thetaCycles,:,2),[],1)];
    % 
    %             sequenceJumps_snippet = [compute_sequenceJumps(x_NaN_snippet(:,2:3)); nan];
    %             posteriorSpread_snippet = compute_dataTemporalConcatenation([x(:,1) posteriorSpread],Times_snippet);
    %             distanceToRat_snippet = sqrt((x_NaN_snippet(:,2)-positions_snippet(:,2)).^2 + (x_NaN_snippet(:,3)-positions_snippet(:,3)).^2 );
    %             distanceToRep_snippet = sqrt((x_NaN_snippet(:,2)-x_NaN_slow_snippet(:,2)).^2 + (x_NaN_snippet(:,3)-x_NaN_slow_snippet(:,3)).^2 );
    % 
    %         subplot(1,3,[1 2])
    %         plot(positions_snippet(:,2),positions_snippet(:,3),'k:','linewidth',2)
    %         hold on, plot(x_NaN_slow_snippet(:,2),x_NaN_slow_snippet(:,3),'r:','linewidth',2), hold off
    %         hold on, cplot(x_NaN_snippet(:,2),x_NaN_snippet(:,3),1:size(x_NaN_snippet,1),'linewidth',5), hold off
    %         % plot_mazeProperties, 
    %         axis square, axis equal, set(gca,'fontsize',14)
    %         time_string = load_timeString((Times_snippet(1)-Run_Times{sessionNum}(1,1))/60);
    % 
    %         for j = 1:length(ind_thetaCycles)
    %             ii = ind_thetaCycles(j);
    %             hold on, quiver(x_NaN_slow_theta(ii,2),x_NaN_slow_theta(ii,3),cos(theta_seqDir(ii,1)),sin(theta_seqDir(ii,1)),6,'r','linewidth',2,'maxHeadSize',0.6), hold off
    %         end
    %         xl = xlim; yl = ylim;
    %         text(xl(1),yl(2),time_string,'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','K','fontsize',14)
    % 
    %         subplot(233)
    %         plot(theta_time_snippet-Times_snippet(1),180/pi*theta_seqDir(ind_thetaCycles,2),'r','linewidth',3)
    %         hold on, plot(theta_time_snippet-Times_snippet(1),180/pi*theta_seqDir(ind_thetaCycles,3),'k','linewidth',3), hold off
    %         ylabel({'relative angle'}), xlabel('time (sec)'), set(gca,'fontsize',14)
    %         hline(0,'k--'), xlim([0 diff(Times_snippet)]), axis square, %pbaspect([4 2 1]), 
    %         legend('rat heading','prev sweep')
    % 
    %         subplot(236)
    %         plot(x_NaN_snippet_origin_rot_full_sub(:,2),x_NaN_snippet_origin_rot_full_sub(:,1),'k.','markersize',16)
    %         axis square, set(gca,'fontsize',14), set(gcf,'color','w')
    %         vline(0,'k--'), hline(0,'k--')
    %         xlabel('x (cm)'), ylabel('y (cm)')
    % 
    % 
    %         % subplot(333)
    %         % yyaxis left, plot(spikeDensity_snippet(:,1)-Times_snippet(1),spikeDensity_snippet(:,3),'k','linewidth',2)
    %         %     xlabel('time (sec)'), ylabel({'spike density','(z-scored)'}), 
    %         %     axis tight, set(gca,'fontsize',14), xlim([0 diff(Times_snippet)]), axis square, %pbaspect([4 2 1]), 
    %         %     yl = ylim;
    %         %     hold on, vline_efficient(theta_time_snippet - Times_snippet(1),yl,[]), hold off
    %         % yyaxis right, plot(positions_snippet(:,1)-Times_snippet(1),positions_snippet(:,4),'b','linewidth',1)
    %         %     ylabel('rat speed (cm/s)')
    %         % 
    %         %     ax = gca;
    %         %     ax.YAxis(1).Color = 'k';
    %         %     ax.YAxis(2).Color = 'b';
    % 
    % 
    %         % for j = 1:length(ind_thetaCycles)
    %         %     ii = ind_thetaCycles(j);
    %         %     x_NaN_snippet_sub = compute_dataTemporalConcatenation(x,theta_timeBounds(ii,:));
    %         %     x_NaN_slow_snippet_sub = compute_dataTemporalConcatenation(x_NaN_slow,theta_timeBounds(ii,:));
    %         %     positions_snippet_sub = compute_dataTemporalConcatenation(positions_snippet,theta_timeBounds(ii,:));
    %         %     LFP_theta_snippet_sub = compute_dataTemporalConcatenation(LFP_theta,theta_timeBounds(ii,:));
    %         % 
    %         %     x_NaN_snippet_cut = x_NaN_snippet_sub;
    %         %     % x_NaN_snippet_cut(LFP_theta_snippet_sub(:,3)<pi,2:3) = nan;
    %         %     x_NaN_snippet_cut(isnan(x_NaN_snippet_cut(:,2)) | isnan(x_NaN_snippet_cut(:,3)),:) = [];
    %         % 
    %         %     distanceToRat_snippet_sub = sqrt((x_NaN_snippet_cut(:,2)-theta_rat(ii,2)).^2 + (x_NaN_snippet_cut(:,3)-theta_rat(ii,3)).^2 );
    %         %     distanceToRep_snippet = sqrt((x_NaN_snippet_cut(:,2)-x_NaN_slow_theta(ii,2)).^2 + (x_NaN_snippet_cut(:,3)-x_NaN_slow_theta(ii,3)).^2 );
    %         % 
    %         %     % %method 1
    %         %     %     delta = diff(x_NaN_snippet_cut(:,2:3),[],1);
    %         %     %         delta(prod(delta,2)==0,:) = [];
    %         %     %     theta_seqDir(ii,1) = circ_mean(atan2(delta(:,2),delta(:,1)));
    %         %     %     theta_seqDir(ii,2) = circ_dist(theta_rat(ii,5),theta_seqDir(ii,1));
    %         %     %     theta_seqDir(ii,3) = circ_dist(theta_seqDir(ii,1),theta_seqDir(ii-1,1));
    %         %     % 
    %         %     % %method 2
    %         %     %     ind_max = find(distanceToRat_snippet_sub==max(distanceToRat_snippet_sub)); ind_max = ind_max(1);
    %         %     %     delta = diff([theta_rat(ii,2:3); x_NaN_snippet_cut(ind_max,2:3)],[],1);
    %         %     %     theta_seqDir(ii,1) = atan2(delta(:,2),delta(:,1));
    %         %     %     theta_seqDir(ii,2) = circ_dist(theta_rat(ii,5),theta_seqDir(ii,1));
    %         %     %     theta_seqDir(ii,3) = circ_dist(theta_seqDir(ii,1),theta_seqDir(ii-1,1));
    %         % 
    %         %     %method 3
    %         %         ind_max = find(distanceToRep_snippet==max(distanceToRep_snippet)); ind_max = ind_max(1);
    %         %         delta = diff([x_NaN_slow_theta(ii,2:3); x_NaN_snippet_cut(ind_max,2:3)],[],1);
    %         %         theta_seqDir(ii,1) = atan2(delta(:,2),delta(:,1));
    %         %         theta_seqDir(ii,2) = circ_dist(theta_rat(ii,5),theta_seqDir(ii,1));
    %         %         theta_seqDir(ii,3) = circ_dist(theta_seqDir(ii-1,1),theta_seqDir(ii,1));
    %         % 
    %         %     subplot(1,3,[1 2])
    %         %     plot(x_NaN_snippet(:,2),x_NaN_snippet(:,3),'color',0.9*ones(1,3),'linewidth',1)
    %         %     hold on, plot(positions_snippet(:,2),positions_snippet(:,3),'color','k','linewidth',1), hold off
    %         %     hold on, cplot(x_NaN_snippet_cut(:,2),x_NaN_snippet_cut(:,3),distanceToRep_snippet,'linewidth',4), hold off, caxis([0 5])
    %         %     hold on, quiver(positions_snippet_sub(1,2),positions_snippet_sub(1,3),cos(theta_seqDir(ii,1)),sin(theta_seqDir(ii,1)),10,'r-','linewidth',3), hold off
    %         % 
    %         %     keyboard
    %         % end
    % 
    %         keyboard
    %     end


    % save('thetaEvents.mat','decoder_replay')
    



            
