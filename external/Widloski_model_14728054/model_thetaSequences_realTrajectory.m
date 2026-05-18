use_theta = 1;
use_feedforward = 1;
use_periodic = 1;

spatialDim = 2;

use_plot = 0;

data_params = struct;

%parameters

    %network parameters
        D = 2;
        n = D*64; %num of neurons per spatialDim
 
        big = 2*n; %padding for convolutions
        N = n^spatialDim;
        dt = 0.5/1000; %step size
        m = 4;

        % v_max = 0.35;
        T_theta_max = 0.12;

        % T_theta_max = 0.12/2;
        % T_replay_min = 0.4;

        tau_s = 20/1000; %synaptic time constant
        beta_0 = 600; %uniform excitation    
        w_rec = 24; %strength of recurrent exc.
        sig_rec = 0.14;%0.02; %width of recurrent exc.
        w_inh = w_rec; %global inhibition
        w_adapt = 40; %amplitude of adaptation 
        tau_adapt = 0.8; %adaptation time constant 
        w_FF = 400; %amplitude of feedforward exc.
        sig_FF = 0.03; %0.001; %width of feedforward exc.

        T_plot = T_theta_max;

        %place field centers
            [X,Y] = meshgrid(D*(1:n)/n,D*(1:n)/n); 

        %activity envelope
            A = ones(n,n);

%rat trajectory

        %increasing velocity
            T_run = 10; %duration of run between stopping periods (sec)
            if use_plot == 1
                numStops = 4; %number of stops
            else
                numStops = 10;
            end

            T = numStops*(T_run); %total simulation time
            L = [2+8 98-8]; %size of box (cm)
            % L = [2 98]; %size of box (cm)
            dt_randomWalk = 0.1; %otherwise, takes too long to generate random walk
            speed_sigma = 0; %step size variance
            theta_sigma = 0.02*2*pi; %orientation variance (each step)

            v_lap = [ones(T_run/dt_randomWalk,1)];

            v = repmat(v_lap,1,numStops); v = v(:);
                v_max_repmat = repmat(linspace(0,2,numStops),length(v_lap),1);
                v = repmat(v_lap,1,numStops).*v_max_repmat; 
                v = v(:);

            x = load_randomWalk(L,T,100*v,speed_sigma,theta_sigma,dt_randomWalk,[],spatialDim);
            x(:,2:3) = x(:,2:3)/100;
            v = [diff(x(:,2)), diff(x(:,3))]/mode(diff(x(:,1))); v = [v(1,:); v];
            speed = sqrt(v(:,1).^2 + v(:,2).^2);
            angle = atan2(v(:,2),v(:,1));
            x = [x speed angle]; x = [0 x(1,2:end); x];
            x = compute_dataInterpolation(x,(dt:dt:T)');
            v = x(:,4);
            x = x(:,2:3);

        % %Real trajectory
        % load behavior_Billy3_07
        % load Experiment_Information.mat
        % T = 1300;
        %     positions_sub = compute_dataTemporalConcatenation(positions,[25742 27041]);
        %     x = positions_sub(:,2); x = x-min(x); x = (x/max(x)*0.8 + 0.1); 
        %     y = positions_sub(:,3); y = y-min(y); y = (y/max(y)*0.8 + 0.1); 
        %     speed = positions_sub(:,4)/100;
        %     t = positions_sub(:,1) - positions_sub(1,1);
        %     positions_sub = compute_dataInterpolation([t x y speed],(0:dt:T)');
        % 
        %     x = positions_sub(:,2:3);
        %     v = positions_sub(:,4);

        % %correlated random walk
        % 
        %     %simulation duration
        %     T_transition = 1;
        %     T_run = 6; %duration of run between stopping periods (sec)
        %     T_stop = 1; %duration of stopping periods
        %     if use_plot == 1
        %         numStops = 4; %number of stops
        %     else
        %         numStops = 200;
        %     end
        % 
        %     T = numStops*(T_run + T_stop + 2*T_transition); %total simulation time
        %     L = [2+8 98-8]; %size of box (cm)
        %     % L = [2 98]; %size of box (cm)
        %     dt_randomWalk = 0.1; %otherwise, takes too long to generate random walk
        %     speed_sigma = 0; %step size variance
        %     theta_sigma = 0.02*2*pi; %orientation variance (each step)
        % 
        %     v_transition = (1+sin(2*pi*(dt_randomWalk:dt_randomWalk:T_transition)'/2/(T_transition) - pi/2))/2;
        % 
        %     %fixed max speed
        %     % v_lap = v_max*[zeros(T_stop/dt_randomWalk,1); v_transition; ones(T_run/dt_randomWalk,1); v_transition(end:-1:1)];
        % 
        %     %variable max speed
        %     v_lap = [zeros(T_stop/dt_randomWalk,1); v_transition; ones(T_run/dt_randomWalk,1); v_transition(end:-1:1)];
        %     v = repmat(v_lap,1,numStops); v = v(:);
        %         v_max_repmat = repmat(abs(normrnd(v_max,0.1,1,numStops)),length(v_lap),1);
        %         v = repmat(v_lap,1,numStops).*v_max_repmat; 
        %         v = v(:);
        % 
        %     x = load_randomWalk(L,T,100*v,speed_sigma,theta_sigma,dt_randomWalk,[],spatialDim);
        %     x(:,2:3) = x(:,2:3)/100;
        %     v = [diff(x(:,2)), diff(x(:,3))]/mode(diff(x(:,1))); v = [v(1,:); v];
        %     speed = sqrt(v(:,1).^2 + v(:,2).^2);
        %     angle = atan2(v(:,2),v(:,1));
        %     x = [x speed angle]; x = [0 x(1,2:end); x];
        %     x = compute_dataInterpolation(x,(dt:dt:T)');
        %     v = x(:,4);
        %     x = x(:,2:3);
        % 
        %     keyboard

        %plot
            % subplot(211)
            % plot(dt:dt:80,v(1:80/dt),'k','linewidth',2)
            % xlabel('time (sec)'), ylabel('rat speed (cm/s)'), set(gca,'fontsize',14), set(gcf,'color','w'), pbaspect([4 1 1])
            % 
            % subplot(212)
            % cplot(x(:,1),x(:,2),v,'linewidth',2), colormap(flipud(parula)),
            % axis square, set(gca,'fontsize',14), xlabel('x (m)'), ylabel('y (m)')
            % cb = colorbar(); ylabel(cb,'rat speed','Rotation',270)
            % 
            % keyboard

    % plot(x)
    % hold on
    % plot(abs(diff(v)/dt))
    % hold off
    % keyboard


%frequency of sequence generation
    %ITI sampling
        %nonstationary distribution: sample theta ITI's based on running speed
            % time = 0:dt:T;
            % sigma_theta = 1.2;%0.5;
            % f_theta = v/max(v)*(1/T_theta_max) + (1-v/max(v))*(1/T_replay_min);
            % isi = 1./normrnd(f_theta(1),sigma_theta);
            % while isi<0
            %     isi = 1./normrnd(f_theta(1),sigma_theta);
            % end
            % t = 1+isi;
            % Theta = t;
            % while t<T
            %     % isi = 1./normrnd(interp1(dt:dt:T,f_theta,t),sigma_theta);
            %     % while isi<0
            %     %     isi = 1./normrnd(interp1(dt:dt:T,f_theta,t),sigma_theta);
            %     % end
            % 
            %     ind = round(t)/dt-5/dt+1:round(t)/dt+5/dt;
            %         ind(ind<=0) = [];
            %         ind(ind>length(time)) = [];
            %     isi = 1./normrnd(interp1(time(ind),f_theta(ind),t),sigma_theta);
            %     while isi<0
            %         isi = 1./normrnd(interp1(time(ind),f_theta(ind),t),sigma_theta);
            %     end
            % 
            %     t = t+isi;
            %     Theta = [Theta; t];
            % end
            % Theta(isnan(Theta)) = [];


        %stationary distribution
            sigma_theta = 0.6;%1.2;%;
            ITI = 1./normrnd((1/T_theta_max),sigma_theta,ceil(T/T_theta_max),1);
                ITI(ITI<=0) = 1/T_theta_max;
            Theta = cumsum(ITI);

            % histogram(ITI,linspace(0,0.25,20),'FaceColor','k','EdgeColor','none','Normalization','probability')
            % axis square, xlabel('theta period (sec)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
            % xlim([0 0.25])
            % keyboard

        %theta modulation
            theta = [Theta(1:end-1) Theta(1:end-1) + 0.5*diff(Theta)];
            % theta_ff = [Theta(1:end-1) - 0.1*diff(Theta) Theta(1:end-1) + 0.1*diff(Theta)];
            theta_ff = [Theta(1:end-1) - 0.04*diff(Theta) Theta(1:end-1) + 0.04*diff(Theta)];
            
                %fix negative theta windows
                ind = find(theta(:,1)<0); if ~isempty(ind), if ind>1, theta(ind,1) = theta(ind-1,2); else, theta(ind,1) = 0; end, end
                ind = find(theta_ff(:,1)<0); if ~isempty(ind), if ind>1, theta_ff(ind,1) = theta_ff(ind-1,2); else, theta_ff(ind,1) = 0; end, end
            
            theta_ind = floor(theta/dt);
            theta_ff_ind = floor(theta_ff/dt);
            g_theta = zeros(T/dt,1);
            g_theta_FF = zeros(T/dt,1);
            for i = 1:size(theta,1)
                g_theta(theta_ind(i,1):theta_ind(i,2)) = 1;
                g_theta_FF(theta_ff_ind(i,1):theta_ff_ind(i,2)) = 1;
            end
            g_theta = g_theta(1:T/dt);
            g_theta_FF = g_theta_FF(1:T/dt);

            %used to modulate facilitation
            g_rest = zeros(T/dt,1);
            g_rest(v<0.01) = 1;


    %plot
        % subplot(211)
        % plot(dt:dt:1,g_theta(1:1/dt),'k','linewidth',2)
        % hold on, plot(dt:dt:1,g_theta_FF(1:1/dt),'r','linewidth',2), hold off
        % legend('recurrent','place'), ylabel('modulation'), xlabel('time (sec)'), box on, set(gca,'fontsize',14), axis tight, pbaspect([4 1 1]), set(gcf,'color','w')
        % 
        % subplot(223)
        % cplot(x(:,1),x(:,2),g_theta_FF,'linewidth',2), colormap((parula))
        % axis square, set(gca,'fontsize',14), xlabel('x (m)'), ylabel('y (m)')
        % cb = colorbar(); ylabel(cb,'Place input modulation','Rotation',270)
        % 
        % subplot(224)
        % cplot(x(:,1),x(:,2),g_theta,'linewidth',2), colormap((parula))
        % axis square, set(gca,'fontsize',14), xlabel('x (m)'), ylabel('y (m)')
        % cb = colorbar(); ylabel(cb,'Recurrent input modulation','Rotation',270)
        % 
        % keyboard

%synaptic weight matrics
    w = mvnpdf([X(:) Y(:)],D*(n/2+1)/n*[1 1],sig_rec^2*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w;% - w_inh; 
    % w = 0.5*mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec^2*eye(2)/2) - 0.95*mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec^2*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w;
    % plot(w(:,end/2)), title(sum(w(:,end/2)))
 
    % I1 = 1/2*(1+sin(1*2*pi*X*3 + 0*2*pi*Y*3)) ; 
    % I2 = 1/2*(1+sin(cos(-120)*-2*pi*X*3 + sin(-120)*2*pi*Y*3)) ; 
    % I3 = 1/2*(1+sin(cos(120)*-2*pi*X*3 + sin(120)*2*pi*Y*3)) ; 
    % imagesc(I1+I2+I3)

    if use_periodic == 1
        w_ft=fft2(fftshift(w));
    else
        w_ft=fft2((w),big,big);
    end

%initialize vectors
    activityCOM = nan(T/dt,spatialDim);
    activitySpread = nan(T/dt,1);
    r = zeros(n,n); a = zeros(n,n);
    spk_count = zeros(N,1);

%simulation
    for t=1:T/dt
        if mod(t,10*T_plot/dt)==0, t/(T/dt), end
        
        %combined input conductance
            if use_periodic == 1
                g_rec = real(ifft2(fft2(r).*w_ft));
            else
                g_rec = real(ifft2(fft2(r,big,big).*w_ft)); g_rec = g_rec(n/2+1:big-n/2,n/2+1:big-n/2);
            end
            g_FF = reshape(mvnpdf([X(:) Y(:)],x(t,:),sig_FF^2*eye(2)),n,n); g_FF = w_FF*g_FF/max(g_FF(:));
            if use_feedforward==0, g_FF = 0; end
            g_feedbackInh = -w_inh*sum(r(:));
            g_unifInput = beta_0;
            g_adapt = -w_adapt*a;

            if use_theta==1
                G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF;
            else
                G = A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_FF;
            end

        %pass conductance variables through nonlinearity to generate output rates, F
            F = 0.*(G<0) + G.*(G>=0);   %rectified-linear transfer function
            
        %track bump properties
            activitySpread(t) = compute_imageSpread(F,2);
      
            %circular center of mass
                bin_centers = linspace(-pi,pi,n);
                center_of_mass = nan(1,spatialDim);
                for i = 1:spatialDim
                    hist_data = nanmean(r,i);
                    hist_data = hist_data / sum(hist_data);   
                    weighted_mean_sin = sum(dot(hist_data,sin(bin_centers)));
                    weighted_mean_cos = sum(dot(hist_data,cos(bin_centers)));            
                    center_of_mass(i) = (atan2(weighted_mean_sin, weighted_mean_cos)+pi)*n/(2*pi);
                end
                activityCOM(t,:) = center_of_mass;

        %update neural activities 
            if use_spiking == 1
                %spikes generated from inhomogeneous Poisson process
                % spk = poissrnd(F*dt);
    
                %subdivide interval m times and take mth spike (for generating spikes with CV = 1/sqrt(m) )
                spk_sub = poissrnd(repmat(F(:),1,m)*dt);
                spk_count = spk_count+sum(spk_sub,2);
                spk = floor(spk_count/m);
                spk_count = rem(spk_count,m);
                if spatialDim==2, spk = reshape(spk,n,n); end
                
                %update firing rates
                r = r + spk - r*dt/tau_s;
                
                %update adaptation dynamics
                a = a + spk - a*dt/tau_adapt;
                
            else
                %update firing rates
                r = r + F*dt - r*dt/tau_s;
                
                %update adaptation dynamics
                a = a + F*dt - a*dt/tau_adapt;
            end
          
        %Plot
        if use_plot==1 &&  mod(t,T_plot/dt)==0 && t>1/dt    
            % subplot(331), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('firing rate'), clim([0 3])
            % subplot(334), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('adaptation (-)'), clim([0 30])
            % subplot(337), plot(dt:dt:T,v), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
            subplot(231), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('Firing rate'), clim([0 3])
            subplot(234), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('Adaptation'), clim([0 30])
            subplot(1,3,[2 3]), 
                tt = t-1/dt:t;
                activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread(:,1)) | [compute_sequenceJumps(activityCOM_sub); 0]>1,:) = nan;
                plot((activityCOM_sub(tt,1)),(activityCOM_sub(tt,2)),'.','color',0.8*ones(1,3),'markersize',8)
                % hold on, cplot((activityCOM_sub(t-1/dt+1:t,1)),(activityCOM_sub(t-1/dt+1:t,2)),1:1/dt,'.','markersize',12), hold off
                hold on, cplot((activityCOM_sub(tt,1)),(activityCOM_sub(tt,2)),tt,'.','markersize',12), hold off
                hold on, plot((n*x(t,1)-0.5)/D,(n*x(t,2)-0.5)/D,'ro','linewidth',2), hold off
                hold on, plot((n*x(tt,1)-0.5)/D,(n*x(tt,2)-0.5)/D,'k:','linewidth',2), hold off
                axis([0 n 0 n]), set(gca,'xtick',[],'ytick',[]), 
                axis square, title('Center-of-mass trajectory')
                % camroll(-90)

            % set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
            % xl = xlim; yl = ylim; text(xl(1)+0.01*diff(xl),yl(1),strcat('elapsed time:',{' '},num2str(t*dt,'%0.2f'),' sec'),'HorizontalAlignment','left','VerticalAlignment','top','fontsize',14);
            % writeVideo(V,getframe(gcf));
            drawnow
        end
            
    end
    % keyboard
    
    if use_plot == 1
        t = dt*(1:T/dt)';
        spk_cell = cell(N,1);
        for i = 1:N
            spk_cell{i} = t(spk_mat(i,:)==1);
            spk_cell{i} = spk_cell{i} + 0.1*(2*rand(size(spk_cell{i}))-1)/2;
        end
        plot_raster_times(spk_cell)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %rat position
            t = dt*(dt:T/dt)';
            angle = [atan2(diff(x(:,2)),diff(x(:,1))); nan];
            rat = [t x v angle];

        %Theta cycles
            [pks,locs] = findpeaks(abs(diff(g_theta)));
            theta_bounds = [locs(1:2:end-1)+1 locs(2:2:end)+1];
            theta_timeBounds = t(theta_bounds);
            theta_time = nanmean(theta_timeBounds,2);
            numThetaCycles = size(theta_bounds,1);

            [pks,locs] = findpeaks(abs(diff(g_theta_FF)));
            theta_bounds_FF = [locs(1:2:end-1)+1 locs(2:2:end)+1]; 

        %rat properties
            id_angleColumns = 5;
            theta_rat = compute_dataInterpolation(rat,theta_time,id_angleColumns);

            [pks,locs] = findpeaks(abs(diff(g_theta_FF)));
            theta_bounds_FF = [locs(1:2:end-1)+1 locs(2:2:end)+1]; 
            theta_timeBounds_FF = t(theta_bounds_FF);
            theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),id_angleColumns);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if use_plot == 0
            data_params.dt = dt;
            data_params.spatialDim = spatialDim;
            data_params.n = n;
            data_params.rat = rat;
            data_params.theta_bounds = theta_bounds;
            data_params.theta_bounds_FF = theta_bounds_FF;
            data_params.activityCOM = activityCOM;
            data_params.T_theta_max = T_theta_max;
            data_params.w_adapt = w_adapt;
            data_params.tau_adapt = tau_adapt;
          
            % save('data_params_2D_increasingSpeed','data_params','-v7.3')
        end

        keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
        % load data_params_2D_realTrajectory_variableTheta
        % load data_params_2D_increasingSpeed
        load data_params_2D_realTrajectory

        currentDir = pwd;
        % cd /home/john/Documents/analysis_ThetaSequences/
        cd /media/john/WorkingData/Data/BarrierMaze/analysis_ThetaSequences/
            load thetaSequences
        cd(currentDir)
        positions = data(1).x;
        
        dt = data_params.dt;
        n = data_params.n;
        rat = data_params.rat; 
            rat(:,1) = rat(:,1) + positions(1);
            rat(:,4) = 100*rat(:,4);
        theta_bounds = data_params.theta_bounds;
        theta_bounds_FF = data_params.theta_bounds_FF;
        activityCOM = data_params.activityCOM;
        T_theta_max = data_params.T_theta_max;
        w_adapt = data_params.w_adapt;
        tau_adapt = data_params.tau_adapt;

        speedThr_runningPeriods = 3;
        durationThr_runningPeriods = 5;
        distanceThr_thetaCycle = 0;

        %Theta cycles
            t = rat(:,1);
            theta_timeBounds = t(theta_bounds);
            theta_time = theta_timeBounds(:,1);%nanmean(theta_timeBounds,2);
            numThetaCycles = size(theta_bounds,1);

        %find running timeBounds
            [FilterA,FilterB]=butter(2,0.02); 
            speed_smooth = compute_butterFilter(positions(:,4),FilterA,FilterB);
            angle_smooth = compute_butterFilter(unwrap(positions(:,9)),FilterA,FilterB);
    
            positions_sub = positions; positions_sub(speed_smooth<speedThr_runningPeriods,:) = nan;
            [boundaries,lengths] = compute_allSequences_NaNseparated(positions_sub(:,1));     
            rat_running_timeBounds = [positions_sub(boundaries(:,1),1) positions_sub(boundaries(:,2),1)];
            durations = rat_running_timeBounds(:,2)-rat_running_timeBounds(:,1);
            rat_running_timeBounds(durations<durationThr_runningPeriods,:) = [];
            numRunningPeriods = size(rat_running_timeBounds,1);

        %rat properties
            theta_rat = compute_dataInterpolation(rat,theta_time,5);
            theta_timeBounds_FF = t(theta_bounds_FF);
            theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),5);

            %angular variability
            [FilterA,FilterB]=butter(2,0.02); 
            % speed_smooth = compute_butterFilter(positions(:,4),FilterA,FilterB);
            % angle_smooth = compute_butterFilter(unwrap(positions(:,9)),FilterA,FilterB);
            % VTE = abs(diff(angle_smooth)./diff(positions(:,1)));
            % theta_VTE = compute_dataInterpolation([positions(1:end-1,1) VTE],theta_time);
            % theta_VTE = theta_VTE(:,2);

            % HD = unwrap(positions(:,9)); HD_diff = abs(diff(HD));
            % zldphi = cumsum(HD_diff); 
            % zldphi_diff = diff(zldphi);
            % zldphi_diff_smooth = filtfilt(FilterA,FilterB,zldphi_diff);
            % zldphi_diff_smooth = [zldphi_diff_smooth(1); zldphi_diff_smooth; zldphi_diff_smooth(end)];
            % theta_VTE = compute_dataInterpolation([positions(:,1) zldphi_diff_smooth],theta_time);
            % theta_VTE = theta_VTE(:,2);

            HD = unwrap(positions(:,9)); 
            HD_diff = abs(diff(HD));
            HD_diff_smooth = filtfilt(FilterA,FilterB,HD_diff);
            zldphi = compute_zscore(HD_diff_smooth);
            zldphi = [zldphi; zldphi(end)];
            theta_VTE = compute_dataInterpolation([positions(:,1) zldphi],theta_time);
            theta_VTE = theta_VTE(:,2);

        %Rat running and stopping periods
            numRunningPeriods = size(rat_running_timeBounds,1);

        %time since rat running
            theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
            theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);

        %sequence direction
            theta_seqDir = nan(numThetaCycles,3);
            theta_seqDist = nan(numThetaCycles,2);
            activityCOM_rot_full = nan(numThetaCycles,400,2);
            activityCOM_rot = nan(size(activityCOM));
            for i = 2:numThetaCycles
                theta_sub = activityCOM(theta_bounds(i,1):theta_bounds(i,2),:);
                rat_sub = rat(theta_bounds(i,1):theta_bounds(i,2),:);

                if sqrt((theta_sub(end,1)-theta_sub(1,1))^2 + (theta_sub(end,2)-theta_sub(1,2))^2)>0.9*n, continue, end

                distanceToRat_snippet = sqrt((theta_sub(:,1)-n*theta_rat(i,2)-0.5).^2 + (theta_sub(:,2)-n*theta_rat(i,3)-0.5).^2 );

                % %method 1
                %     delta = diff(theta_sub,[],1);
                %         delta(prod(delta,2)==0,:) = [];
                %     theta_seqDir(i,1) = circ_mean(atan2(delta(:,2),delta(:,1)));
                %     theta_seqDir(i,2) = circ_dist(theta_rat(i,5),theta_seqDir(i,1));
                %     theta_seqDir(i,3) = circ_dist(theta_seqDir(i-1,1),theta_seqDir(i,1));

                %method 2
                    % ind_min = find(distanceToRat_snippet==min(distanceToRat_snippet)); ind_min = ind_min(1);
                    % ind_max = find(distanceToRat_snippet==max(distanceToRat_snippet)); ind_max = ind_max(1);
                    % ind_min = 1;
                    % ind_max = length(distanceToRat_snippet);
                    % delta = diff([n*theta_rat(i,2:3)-0.5; theta_sub(end,:)],[],1);
                    delta = diff([theta_sub(1,:); theta_sub(end,:)],[],1);
                    theta_seqDir(i,1) = atan2(delta(:,2),delta(:,1));
                    theta_seqDir(i,2) = circ_dist(theta_rat(i,5),theta_seqDir(i,1));
                    theta_seqDir(i,3) = circ_dist(theta_seqDir(i-1,1),theta_seqDir(i,1));

                %sequence distance
                    theta_seqDist(i,1) = compute_sequenceDistance(unwrap(theta_sub))/n;
                    theta_seqDist(i,2) = sqrt((theta_sub(end,1)-theta_sub(1,1))^2 + (theta_sub(end,2)-theta_sub(1,2))^2)/n;

                %rotate theta sequence by rat heading direction
                    theta_sub_origin = theta_sub - repmat(n*theta_rat(i,2:3)-0.5,size(theta_sub,1),1);
                    M = [cos(-theta_rat(i,5)) -sin(-theta_rat(i,5)); sin(-theta_rat(i,5)) cos(-theta_rat(i,5))];
                    theta_sub_origin_rot = (M*theta_sub_origin')';

                    activityCOM_rot(theta_bounds(i,1):theta_bounds(i,2),:) = theta_sub_origin_rot;
                    activityCOM_rot_full(i,:,1) = compute_vecBuffered(theta_sub_origin_rot(:,1),400);
                    activityCOM_rot_full(i,:,2) = compute_vecBuffered(theta_sub_origin_rot(:,2),400);

                    % if theta_seqDist(i,1)>0.01
                    %     hold on, plot(theta_sub(:,1),theta_sub(:,2),'k.')
                    %     hold on, plot(n*theta_rat(i,2)-0.5,n*theta_rat(i,3)-0.5,'ro')
                    %     keyboard
                    % end
            end

        linewidth = 4;
        fontsize = 16;
        gap = 0.05;

   %angle variability vs rat speed
        % [edges_speed,centers_speed] = load_timeBins([5 60],2,8);
        [edges_speed,centers_speed] = load_timeBins([5 60],3,3);
        edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
        H1_data = nan(length(centers_speed),2);
        H2_data = nan(length(centers_speed),2);
        P_data = nan(length(centers_speed),2);
        for i = 1:length(centers_speed)
            ind_thetaCycles = find(theta_rat(:,4)>edges_speed(i) & theta_rat(:,4)<=edges_speed(i+1));
            % if length(ind_thetaCycles)<20, continue, end

            H1_data(i,1) = nanmean(abs(theta_seqDir(ind_thetaCycles,2)));
            H1_data(i,2) = nanstd(abs(theta_seqDir(ind_thetaCycles,2)))./sqrt(length(ind_thetaCycles));
            H2_data(i,1) = nanmean(abs(theta_seqDir(ind_thetaCycles,3)));
            H2_data(i,2) = nanstd(abs(theta_seqDir(ind_thetaCycles,3)))./sqrt(length(ind_thetaCycles));
        
            P_data(i,1) = nanmean(abs(sign(theta_seqDir(ind_thetaCycles,3)) + sign(theta_seqDir(ind_thetaCycles+1,3)))==0);
        end

        subplot(221)
        plot(0.01*centers_speed,100*P_data(:,1),'k','linewidth',linewidth)
        xlabel('rat speed (m/sec)'), ylabel('% alternation')
        axis square, axis tight;  box on, 
        yl = ylim; xl = xlim;
        xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
        ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
        xlim([0 0.6]), set(gca,'xtick',0:0.2:0.6), 

        subplot(222)
        plot(0.01*centers_speed,180/pi*H1_data(:,1),'k','linewidth',linewidth)
        xlabel('rat speed (m/sec)'), ylabel('|sweep angle| (deg)')
        axis square, axis tight;  box on, 
        yl = ylim; xl = xlim;
        xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
        ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
        xlim([0 0.6]), set(gca,'xtick',0:0.2:0.6), 

     %angle variance vs rat angle variability
        % [edges_VTE,centers_VTE] = load_timeBins([-1.8 1.8],0.2,1);
        [edges_VTE,centers_VTE] = load_timeBins([-1.8 1.8],0.3,0.3);
        edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
        H1_data = nan(length(centers_VTE),2);
        H2_data = nan(length(centers_VTE),2);
        P_data = nan(length(centers_VTE),2);
        numData = nan(length(centers_VTE),1);
        for i = 1:length(centers_VTE)
            ind_thetaCycles = find(theta_VTE>edges_VTE(i,1) & theta_VTE<=edges_VTE(i,2));
            numData(i) = length(ind_thetaCycles);
            % if length(ind_thetaCycles)<100, continue, end

            H1_data(i,1) = nanmean(abs(theta_seqDir(ind_thetaCycles,2)));
            H1_data(i,2) = nanstd(abs(theta_seqDir(ind_thetaCycles,2)))./sqrt(length(ind_thetaCycles));
            H2_data(i,1) = nanmean(abs(theta_seqDir(ind_thetaCycles,3)));
            H2_data(i,2) = nanstd(abs(theta_seqDir(ind_thetaCycles,3)))./sqrt(length(ind_thetaCycles));

            P_data(i,1) = nanmean(abs(sign(theta_seqDir(ind_thetaCycles(1:end-1),3)) + sign(theta_seqDir(ind_thetaCycles(1:end-1)+1,3)))==0);
        end

        subplot(223)
        plot(centers_VTE,100*P_data(:,1),'k','linewidth',linewidth)
        xlabel('rat VTE'), ylabel('% alternation')
        axis square, axis tight;  box on, 
        yl = ylim; xl = xlim;
        xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
        ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
        xlim([-2 2]), set(gca,'xtick',-2:1:2), 

        subplot(224)
        plot(centers_VTE,180/pi*H1_data(:,1),'k','linewidth',linewidth)
        xlabel('rat VTE'), ylabel('|sweep angle| (deg)')
        axis square, axis tight;  box on, 
        yl = ylim; xl = xlim;
        xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
        ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
        xlim([-2 2]), set(gca,'xtick',-2:1:2), 


        set(gcf,'color','w')
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontsize); 

%%







            %temporal regularity
                % edges_time = linspace(0,2,60); centers_time = edges_time(1:end-1) + mean(diff(edges_time))/2;
                [edges_time,centers_time] = load_timeBins([0 2],0.06,0.06);
                delT = squareform(pdist(theta_rat(ind_thetaCycles,1))); delT = delT(:);
                delX = squareform(pdist(sign(theta_seqDir(ind_thetaCycles,2)))); delX = delX(:)/2;
                delX(delT>max(edges_time(:))) = []; delT(delT>max(edges_time(:))) = [];
                data_pairwise_opposite = nan(3,length(centers_time));
                for i = 1:length(centers_time)
                    ind = find(delT>edges_time(i) & delT<=edges_time(i,2));
                    data_pairwise_opposite(1,i) = nanmean(delX(ind));
                    data_pairwise_opposite(2,i) = nanstd(delX(ind))./sqrt(length(ind));
                    data_pairwise_opposite(3,i) = length(ind);
                end

            %Evolution across theta sequences
                speedThr_thetaCycle = 0;
                maxSequence = 40;
                % [timeBins_evolution,centers_timeBins_evolution] = load_timeBins([0 6],0.1,1);
                [timeBins_evolution,centers_timeBins_evolution] = load_timeBins([0.3 5],0.2,0.6);
                data_relativeAngleToRat_evolution = nan(numRunningPeriods,size(timeBins_evolution,1));
                data_relativeAngleToPrev_evolution = nan(numRunningPeriods,size(timeBins_evolution,1));
                data_ratSpeed_evolution = nan(numRunningPeriods,size(timeBins_evolution,1));
                data_pairwise_opposite_evolution  = nan(numRunningPeriods,size(timeBins_evolution,1));
                data_relativeAngleToRat_evolution_mat = nan(numRunningPeriods,maxSequence);
                data_relativeAngleToPrev_evolution_mat = nan(numRunningPeriods,maxSequence);
                for i = 1:size(rat_running_timeBounds,1)
                    Times_snippet = rat_running_timeBounds(i,:);
                    data_sub = compute_dataTemporalConcatenation([theta_time theta_timeSinceRunning theta_seqDir theta_rat(:,4) theta_seqDist],Times_snippet);
                    data_sub(data_sub(:,6)<speedThr_thetaCycle,:) = [];
                    % data_sub(data_sub(:,7)<distanceThr_thetaCycle,:) = [];
                    if size(data_sub,1)<5, continue, end
        
                    for j =  1:size(timeBins_evolution,1)
                        data_subsub = compute_dataTemporalConcatenation(data_sub,rat_running_timeBounds(i,1)+timeBins_evolution(j,:));
        
                        data_relativeAngleToRat_evolution(i,j) = nanmean(abs(data_subsub(:,4)));
                        data_relativeAngleToPrev_evolution(i,j) = nanmean(abs(data_subsub(:,5)));
        
                        data_ratSpeed_evolution(i,j) = nanmean(abs(data_subsub(:,6)));
        
                        data_sub_sub = sign(data_subsub(:,5));
                        data_sub_sub = abs(data_sub_sub(1:end-1)+data_sub_sub(2:end))==0;
                        data_pairwise_opposite_evolution(i,j) = nanmean(data_sub_sub);
                    end
                    data_relativeAngleToRat_evolution_mat(i,:) = compute_vecBuffered(data_sub(:,4),maxSequence)';
                    data_relativeAngleToPrev_evolution_mat(i,:) = compute_vecBuffered(data_sub(:,5),maxSequence)';
                end



        speedThr_thetaCycle = 20;
        fontsize = 16;
        linewidth = 4;
        gap = 0.05;

        %theta period
            subplot(321)
            theta_durations = theta_timeBounds(:,2) - theta_timeBounds(:,1);
            % edges_thetaDuration = linspace(0,0.25,20); centers_thetaDuration = edges_thetaDuration(1:end-1) + mean(diff(edges_thetaDuration))/2;
            edges_thetaDuration = 0.05:0.01:0.2; centers_thetaDuration = edges_thetaDuration(1:end-1) + mean(diff(edges_thetaDuration))/2;
            h = histcounts(2*theta_durations(ind_thetaCycles),edges_thetaDuration,'Normalization','probability');
            plot(centers_thetaDuration,h,'k','linewidth',linewidth)
            xlabel('theta period (sec)'), ylabel('fraction'), 
            axis square, axis tight;  box on, 
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])

        %sequence distance
            ind_thetaCycles = find(100*theta_rat(:,4)>speedThr_thetaCycle);
            subplot(322)
                edges_distance = 0:2:45; centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
            h = histcounts(100*theta_seqDist(ind_thetaCycles,1),edges_distance,'Normalization','probability');
            plot(0.01*centers_distance,h,'k','linewidth',linewidth)
            xlabel('sweep distance (m)'), ylabel('fraction'), 
            axis square, axis tight;  box on, 
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])

        %temporal regularity
            subplot(323)
            plot(centers_time,100*data_pairwise_opposite(1,:),'k','linewidth',linewidth)
            xlabel('time between sweeps (sec)'), ylabel('% alternation')
            axis square, axis tight;  box on, 
            hline(50,'k--'),
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
            ylim([10 90]);

        %angular regularity
            subplot(324)
            h1 = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges_angles,'normalization','probability');
            h2 = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges_angles,'normalization','probability');
            h1 = smooth(h1,3); h2 = smooth(h2,3);
            yyaxis right, 
                plot(180/pi*centers_angles,h2,'k','linewidth',linewidth)
                set(gca,'xtick',[0 60 120 180]), 
                axis square, axis tight;  box on, 
                yl = ylim; xl = xlim;
                xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
                ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])

            yyaxis left, 
                plot(180/pi*centers_angles,h1,'r','linewidth',linewidth)
                set(gca,'xtick',[0 60 120 180]), 
                ylabel('fraction')
                axis square, axis tight;  box on, 
                yl = ylim; xl = xlim;
                xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
                ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])

                legend('prev. sweep','rat heading','location','southwest')

            ax = gca;
            ax.YAxis(1).Color = 'r';
            ax.YAxis(2).Color = 'k';
            set(gca,'xtick',0:60:180)

        %temporal regularity across time
            subplot(325)
            yyaxis left, plot(centers_timeBins_evolution,100*nanmean(data_pairwise_opposite_evolution),'k','linewidth',linewidth)
            ylabel('% alternation'), 
            axis square, axis tight;  box on, 
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
            xlim([0 5])
            
            yyaxis right, plot(centers_timeBins_evolution,nanmean(data_ratSpeed_evolution),'b','linewidth',linewidth)
            xlabel({'time since start','of movement bout (sec)'}), ylabel('rat speed (m/sec)')
            axis square, axis tight;  box on, 
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])
            xlim([0 5])

            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'b';
            set(gca,'xtick',0:2:4)

        %angular regularity across time
            subplot(326)
            yyaxis left, plot(centers_timeBins_evolution,180/pi*(nanmean(data_relativeAngleToRat_evolution)),'k','linewidth',linewidth)
            ylabel('|sweep angle| (deg)')
            axis square, axis tight;  box on, 
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])

            yyaxis right, plot(centers_timeBins_evolution,nanmean(data_ratSpeed_evolution),'b','linewidth',linewidth)
            xlabel({'time since start','of movement bout (sec)'}), ylabel('rat speed (m/sec)')
            axis square, axis tight;  box on, 
            yl = ylim; xl = xlim;
            xlim([xl(1)-diff(xl)*gap xl(2)+diff(xl)*gap])
            ylim([yl(1)-diff(yl)*gap yl(2)+diff(yl)*gap])

            ax = gca;
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Color = 'b';

        set(gcf,'color','w'),
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', fontsize); 



%%



    %     %find running timeBounds
    %         [FilterA,FilterB]=butter(2,0.02); 
    %         speed_smooth = compute_butterFilter(positions(:,4),FilterA,FilterB);
    % 
    %         positions_sub = positions; positions_sub(speed_smooth<speedThr_runningPeriods,:) = nan;
    %         [boundaries,lengths] = compute_allSequences_NaNseparated(positions_sub(:,1));     
    %         rat_running_timeBounds = [positions_sub(boundaries(:,1),1) positions_sub(boundaries(:,2),1)];
    %         durations = rat_running_timeBounds(:,2)-rat_running_timeBounds(:,1);
    %         rat_running_timeBounds(durations<durationThr_runningPeriods,:) = [];
    %         numRunningPeriods = size(rat_running_timeBounds,1);
    % 
    %         % T = 20000;
    %         % posSampRate = 1/dt;
    %         % data_ratSpeed_evolution = nan(numRunningPeriods,T);
    %         % for i = 1:numRunningPeriods
    %         %     positions_running = compute_dataTemporalConcatenation([rat(:,1),rat(:,4)],rat_running_timeBounds(i,:));
    %         %     TT = min(T,size(positions_running,1));
    %         %     data_ratSpeed_evolution(i,1:TT) = positions_running(1:TT,2);
    %         % end
    %         % plot(1/posSampRate*(1:T),nanmean(data_ratSpeed_evolution))
    % 
    %     %time since rat running
    %         theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
    %         theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
    %         theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);
    % 
    %     %Evolution across theta sequences
    %         speedThr_thetaCycle = 0;
    %         maxSequence = 40;
    %         % [thetaSequenceNum_evolution,centers_thetaSequenceNum_evolution] = load_timeBins([1 30],1,3); 
    %         [thetaSequenceNum_evolution,centers_thetaSequenceNum_evolution] = load_timeBins([1 maxSequence],1,5); 
    %         data_relativeAngleToRat_evolution = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_relativeAngleToPrev_evolution = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_ratSpeed_evolution = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_pairwise_opposite_evolution  = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_relativeAngleToRat_evolution_mat = nan(numRunningPeriods,maxSequence);
    %         data_relativeAngleToPrev_evolution_mat = nan(numRunningPeriods,maxSequence);
    %         for i = 1:size(rat_running_timeBounds,1)
    %             Times_snippet = rat_running_timeBounds(i,:);
    %             data_sub = compute_dataTemporalConcatenation([theta_time theta_timeSinceRunning theta_seqDir theta_rat(:,4) 100*theta_seqDist(:,1)],Times_snippet);
    %             data_sub(data_sub(:,6)<speedThr_thetaCycle,:) = [];
    %             data_sub(data_sub(:,7)<distanceThr_thetaCycle,:) = [];
    %             if size(data_sub,1)<10, continue, end
    % 
    %             data_sub = compute_vecBuffered(data_sub,maxSequence);
    %             for j =  1:length(thetaSequenceNum_evolution)
    %                 data_relativeAngleToRat_evolution(i,j) = nanmean(abs(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),4)));
    %                 data_relativeAngleToPrev_evolution(i,j) = nanmean(abs(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),5)));
    % 
    %                 data_ratSpeed_evolution(i,j) = nanmean(abs(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),6)));
    % 
    %                 data_sub_sub = sign(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),5));
    %                 data_sub_sub = abs(data_sub_sub(1:end-1)+data_sub_sub(2:end))==0;
    %                 data_pairwise_opposite_evolution(i,j) = nanmean(data_sub_sub);
    %             end
    %             data_relativeAngleToRat_evolution_mat(i,:) = compute_vecBuffered(data_sub(:,4),maxSequence)';
    %             data_relativeAngleToPrev_evolution_mat(i,:) = compute_vecBuffered(data_sub(:,5),maxSequence)';
    %         end
    % 
    % 
    % %plot
    % %%%%%%%%%%%%%%%%%%%%%%
    % 
    % speedThr_thetaCycle = 20;
    % ind_thetaCycles = find(theta_rat(:,4)>speedThr_thetaCycle);
    % 
    % %sequence distance
    %     ind_thetaCycles = find(theta_rat(:,4)>speedThr_thetaCycle);
    %     subplot(331)
    %     edges_distance = linspace(0,40,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
    %     histogram(100*theta_seqDist(ind_thetaCycles,1),edges_distance,'facecolor','k','edgecolor','none','Normalization','probability');
    %     axis square, xlabel('distance (cm)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
    %     xlim([min(edges_distance) max(edges_distance)])
    % 
    % %theta period
    %     subplot(332)
    %     theta_durations = theta_timeBounds(:,2) - theta_timeBounds(:,1);
    %     edges_thetaDuration = linspace(0,0.25,20); centers_thetaDuration = edges_thetaDuration(1:end-1) + mean(diff(edges_thetaDuration))/2;
    %     histogram(theta_durations(ind_thetaCycles),edges_thetaDuration,'facecolor','k','edgecolor','none','Normalization','probability');
    %     axis square, xlabel('period (s)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
    %     xlim([min(edges_thetaDuration) max(edges_thetaDuration)])
    % 
    % %temporal regularity
    %     subplot(334)
    %     delT = squareform(pdist(theta_rat(ind_thetaCycles,1))); delT = delT(:);
    %     delX = squareform(pdist(sign(theta_seqDir(ind_thetaCycles,2)))); delX = delX(:)/2;
    %     edges_time = linspace(0,2,60); centers_time = edges_time(1:end-1) + mean(diff(edges_time))/2;
    %     data_sub = nan(3,length(edges_time)-1);
    %     for i = 1:length(edges_time)-1
    %         ind = find(delT>edges_time(i) & delT<=edges_time(i+1));
    %         data_sub(1,i) = nanmean(delX(ind));
    %         data_sub(2,i) = nanstd(delX(ind))./sqrt(length(ind));
    %         data_sub(3,i) = length(ind);
    %     end
    %     shadedErrorBar(centers_time,100*data_sub(1,:),100*data_sub(2,:),'lineprops',{'k','linewidth',2})
    %     % plot(centers_time,100*data_sub(1,:),'k','linewidth',4)
    %     xlabel('time between events (sec)'), ylabel('% alternation')
    %     hline(50,'k--'), axis square, box on, axis tight
    %     ylim([20 80])
    % 
    % %angular regularity
    %     subplot(335)
    %     edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
    %     h1 = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges_angles,'normalization','probability');
    %     h2 = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges_angles,'normalization','probability');
    %     h1 = smooth(h1,3); h2 = smooth(h2,3);
    %     yyaxis left, 
    %         plot(180/pi*centers_angles,h2,'r','linewidth',3)
    %         ylabel('fraction')
    %     yyaxis right, 
    %         plot(180/pi*centers_angles,h1,'k','linewidth',3)
    %     xlabel('angle (deg)'), axis square, set(gca,'xtick',[0 90 180]), set(gcf,'color','w'), xlim(180*[0 1])
    %     ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; 
    %     legend('prev. sweep','rat heading','location','southwest')
    %     ax = gca;
    %     ax.YAxis(1).Color = 'k';
    %     ax.YAxis(2).Color = 'r';
    %     set(gca,'xtick',0:60:180)
    % 
    % %angular regularity across time
    %     subplot(337)
    %     yyaxis left, shadedErrorBar(centers_thetaSequenceNum_evolution,180/pi*nanvar(data_relativeAngleToRat_evolution),180/pi*nanstd(data_relativeAngleToRat_evolution)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
    %     % yyaxis left, plot(centers_thetaSequenceNum_evolution,180/pi*nanmean(data_relativeAngleToRat_evolution),'k-','linewidth',4)
    %         % hold on, plot(centers_thetaSequenceNum_evolution,180/pi*nanmean(data_relativeAngleToPrev_evolution),'r-','linewidth',4), hold off
    %     ylabel('angle variance (deg)')
    %     % yyaxis right, plot(centers_thetaSequenceNum_evolution,nanmean(data_ratSpeed_evolution),'b','linewidth',4)
    %     yyaxis right, shadedErrorBar(centers_thetaSequenceNum_evolution,0.01*nanmean(data_ratSpeed_evolution),0.01*nanstd(data_ratSpeed_evolution)./sqrt(numRunningPeriods),'lineprops',{'b','linewidth',2})
    %     xlabel('sweep number'), ylabel('rat speed')
    %     axis square, box on, set(gcf,'color','w'),
    %     ax = gca;
    %     ax.YAxis(1).Color = 'k';
    %     ax.YAxis(2).Color = 'b';
    %     % ax.YAxis(1).Limits = [5 43];
    %     % ax.YAxis(2).Limits = [0 0.2];
    % 
    % %temporal regularity across time
    %     subplot(338)
    %     yyaxis left, shadedErrorBar(centers_thetaSequenceNum_evolution,100*nanmean(data_pairwise_opposite_evolution),100*nanstd(data_pairwise_opposite_evolution)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
    %     % yyaxis left, plot(centers_thetaSequenceNum_evolution,100*nanmean(data_pairwise_opposite_evolution),'k','linewidth',4)
    %     ylabel('% alternation'), 
    %     % yyaxis right, plot(centers_thetaSequenceNum_evolution,nanmean(data_ratSpeed_evolution),'b','linewidth',4)
    %     yyaxis right, shadedErrorBar(centers_thetaSequenceNum_evolution,0.01*nanmean(data_ratSpeed_evolution),0.01*nanstd(data_ratSpeed_evolution)./sqrt(numRunningPeriods),'lineprops',{'b','linewidth',2})
    %     xlabel('sweep number'), ylabel('rat speed')
    %     axis square, box on, set(gcf,'color','w'),
    %     ax = gca;
    %     ax.YAxis(1).Color = 'k';
    %     ax.YAxis(2).Color = 'b';
    %     % ax.YAxis(1).Limits = [5 55];
    %     % ax.YAxis(2).Limits = [0 0.18];
    % 
    % %relative angle vs rat speed
    %     [edges_speed,centers_speed] = load_timeBins([0 100],2.5,10);
    %     edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
    %     a = cool(length(edges_speed));
    %     H1 = []; H2 = [];
    %     H1_data = nan(length(centers_speed),2);
    %     H2_data = nan(length(centers_speed),2);
    %     for i = 1:length(centers_speed)
    %         ind_thetaCycles = find(theta_rat(:,4)>edges_speed(i,1) & theta_rat(:,4)<=edges_speed(i,2));
    %         % if length(ind_thetaCycles)<20, continue, end
    %         H1_data(i,1) = nanmean(abs(theta_seqDir(ind_thetaCycles,2)));
    %         H1_data(i,2) = nanstd(abs(theta_seqDir(ind_thetaCycles,2)))./sqrt(length(ind_thetaCycles));
    %         H2_data(i,1) = nanmean(abs(theta_seqDir(ind_thetaCycles,3)));
    %         H2_data(i,2) = nanstd(abs(theta_seqDir(ind_thetaCycles,3)))./sqrt(length(ind_thetaCycles));
    % 
    %         h1 = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges_angles,'normalization','probability');
    %         h2 = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges_angles,'normalization','probability');
    %         h1 = smooth(h1,3)'; h2 = smooth(h2,3)';
    %         H1 = [H1, h1'/max(h1)]; H2 = [H2, h2'/max(h2)];
    %     end
    %     subplot(339)
    %     % h = pcolor(centers_speed,180/pi*centers_angles,H2); set(h,'edgecolor','none')
    %     % plot(centers_speed,180/pi*H1_data(:,1),'k','linewidth',4)
    %     shadedErrorBar(0.01*centers_speed,180/pi*H1_data(:,1),180/pi*H1_data(:,2),'lineprops',{'k','linewidth',2})
    %     % hold on, plot(centers_speed,180/pi*H2_data(:,1),'r','linewidth',4), hold off
    %     xlabel('rat speed'), ylabel({'relative angle to','rat heading (deg)'})
    %     axis square, set(gca,'fontsize',14), set(gcf,'color','w')
    %     box on
    % 
    %     set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14); 
