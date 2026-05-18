clc
clear all
close all
%% Script Overview (Wenbo Tang, March 10, 2026):
% this is the original model, demonstrating L-R alternating sweeps during 
% random foraging, kindly shared by John Widloski
% Reference: Widloski et al., Cell Rep, 2025
%% add the code basepath
addpath(genpath('E:\Tang_latentseq_NN_2026\'))
%% set the model
use_spiking = 1; %set to 1 for stochastic, 0 for deterministic
load('position_info.mat') % load previous generated random walk positions, for reproducibility 
colormap('cool')
%parameters
    n = 64; %num of neurons per spatialDim
    big = 2*n; %padding for convolutions
    N = n^2;
    dt = 0.5/1000; %step size
    m = 4; %coefficient of variance of spiking: CV = 1/sqrt(m)

    tau_s = 20/1000; %synaptic time constant
    beta_0 = 600; %uniform excitation    
    w_rec = 24; %strength of recurrent exc.
    sig_rec = 0.14; %width of recurrent exc.
    w_inh = w_rec; %global inhibition
    w_adapt = 60; %amplitude of adaptation 
    tau_adapt = 0.8; %adaptation time constant 
    w_FF = 400;%2000; %amplitude of feedforward exc.
    sig_FF = 0.03; %0.001; %width of feedforward exc.
    v_max = 0.4;
    T_theta_max = 0.12;

    [X,Y] = meshgrid((1:n)/n,(1:n)/n); %place field centers

%generate rat trajectory as correlated random walk
    T_transition = 1;
    T_run = 6; %duration of run between stopping periods (sec)
    T_stop = 1; %duration of stopping periods
    numStops = 1; %number of running cycles
    T = numStops*(T_run + T_stop + 2*T_transition); %total simulation time

    L = [2+8 98-8]; %size of box (cm)
    dt_coarse = 0.1; %use coarse time bins, takes random walk generation (below) takes too long
    
    %speed profile for 1 running cycle
    v_transition = (1+sin(2*pi*(dt_coarse:dt_coarse:T_transition)'/2/(T_transition) - pi/2))/2;
    v_lap = v_max*[zeros(T_stop/dt_coarse,1); v_transition; ones((T_run + T_transition)/dt_coarse,1)];
    v = repmat(v_lap,1,numStops); v = v(:);

    spatialDim = 2;
    speed_sigma = 0; %step size variance
    theta_sigma = 0.02*2*pi; %orientation variance (each step)
    % x = load_randomWalk(L,T,100*v,speed_sigma,theta_sigma,dt_coarse,[],spatialDim);
    % use previously generated random walk instead for reproducibility
    x = x_walk;
    x(:,2:3) = x(:,2:3)/100; %rescale positions between 0 and 1
    v = [diff(x(:,2)), diff(x(:,3))]/mode(diff(x(:,1))); v = [v(1,:); v];
    speed = sqrt(v(:,1).^2 + v(:,2).^2);
    angle = atan2(v(:,2),v(:,1));
    x = [x speed angle]; x = [0 x(1,2:end); x];

    %resample trajectory at higher temporal frequency
    x = compute_dataInterpolation(x,(dt:dt:T)');
    v = x(:,4);
    x = x(:,2:3);

%theta modulation
    %sample theta periods
    sigma_theta = 1.2;
    ITI = 1./normrnd((1/T_theta_max),sigma_theta,ceil(T/T_theta_max),1);
        ITI(ITI<=0) = 1/T_theta_max;
    Theta = cumsum(ITI);
    
    %internal and feedforward modulation
    theta_int = [Theta(1:end-1) Theta(1:end-1) + 0.5*diff(Theta)];
    theta_ff = [Theta(1:end-1) - 0.1*diff(Theta) Theta(1:end-1) + 0.1*diff(Theta)];

        %fix negative theta windows
        ind = find(theta_int(:,1)<0); if ~isempty(ind), if ind>1, theta_int(ind,1) = theta_int(ind-1,2); else, theta_int(ind,1) = 0; end, end
        ind = find(theta_ff(:,1)<0); if ~isempty(ind), if ind>1, theta_ff(ind,1) = theta_ff(ind-1,2); else, theta_ff(ind,1) = 0; end, end
    
    %define theta conductances across time
    theta_ind = floor(theta_int/dt);
    theta_ff_ind = floor(theta_ff/dt);
    g_theta_int = zeros(T/dt,1);
    g_theta_FF = zeros(T/dt,1);
    for i = 1:size(theta_int,1)
        g_theta_int(theta_ind(i,1):theta_ind(i,2)) = 1;
        g_theta_FF(theta_ff_ind(i,1):theta_ff_ind(i,2)) = 1;
    end
    g_theta_int = g_theta_int(1:T/dt);
    g_theta_FF = g_theta_FF(1:T/dt);
    
        %for plotting
        g_theta_nan = g_theta_int;
        g_theta_nan(g_theta_int==0) = nan;

%synaptic weight matrics
    w = mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec^2*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w;% - w_inh; 

%initialize vectors
    activityCOM = nan(T/dt,spatialDim);
    activitySpread = nan(T/dt,1);
    r = zeros(n,n); a = zeros(n,n);
    spk_count = zeros(N,1);

%% simulation
    for t=1:(T-2)/dt
        
        %nput conductances
        g_rec = real(ifft2(fft2(r).*fft2(fftshift(w))));
        g_feedbackInh = -w_inh*sum(r(:));
        g_unifInput = beta_0;
        g_adapt = -w_adapt*a;
        g_FF = reshape(mvnpdf([X(:) Y(:)],x(t,:),sig_FF^2*eye(2)),n,n); g_FF = w_FF*g_FF/max(g_FF(:));

        %combined input conductance
        G = g_theta_int(t)*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF;

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
                
                spk = reshape(spk,n,n);
                
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
        if mod(t,T_theta_max/dt)==0   
            ax1 = subplot(321); imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('Firing rate'), colormap(ax1,"parula");
            if max(max(r)) == 0
                clim([0 2])
            else 
                clim([0  max(max(r))])
            end
            aa = g_rec + g_adapt + g_unifInput + g_feedbackInh;
            ax2 = subplot(322); imagesc(aa), set(gca,'ydir','normal'), axis off, axis square, title('Excitation'),  clim([-1000 1000]),colormap(ax2,"parula");
            ax3 = subplot(3,1,[2 3]); colormap(ax3,'cool');
                activityCOM_sub = g_theta_nan.*activityCOM; 
                activityCOM_sub(isnan(activitySpread(:,1)) | [compute_sequenceJumps(activityCOM_sub); 0]>1,:) = nan;
                [lo,hi]= findcontiguous(find(~isnan(activityCOM_sub(:,1))));
                for ss = 1:length(lo) 
                    seg = lo(ss):hi(ss);
                    cplot((activityCOM_sub(seg,1)+3),(activityCOM_sub(seg,2))+3,seg-lo(ss)+1,'.','markersize',12), hold on
                end
                hold on, plot((n*x(t,1)+2.5),(n*x(t,2)+2.5),'ko','linewidth',2), hold off
                hold on, plot((n*x(1:t,1)+2.5),(n*x(1:t,2)+2.5),'k','linewidth',2), hold off
                axis([0 n 0 n]), set(gca,'xtick',[],'ytick',[]), 
                axis square, 
                xlabel('neuron'), ylabel('neuron')
            set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
            clim([0 150])
            drawnow
        end
            
    end
    

