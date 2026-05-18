use_spiking = 1;
use_theta = 1;
use_feedforward = 1;
use_envelope = 0;
use_periodic = 1;
use_facilitation = 1;

spatialDim = 2;

use_plot = 1;

% V = VideoWriter('thetaSequences_2.avi','Motion JPEG AVI'); V.Quality = 75;
% open(V);

if spatialDim ==1
    tau_adapt_vec = 0.4:0.4:4; %[0.8,3];
    w_adapt_vec = 0.5:0.5:5; %2:2:16;
    w_facil_vec = 1:2;
    tau_facil_vec = 8:8:24; 
else
    v_max_vec = 0.1:0.1:2;%0.15:0.05:0.6; %0.1:0.1:1;
    T_theta_max_vec = 0.08:0.02:0.26;
    w_adapt_vec = 10:10:100;
    tau_adapt_vec = 0.4:0.4:4; 
    sig_rec_vec = 0.4:0.2:2.2;
    sig_FF_vec = 0.4:0.2:2.2;
end

data_params = struct;
% load data_params_2D_bigEnvironment_ratSpeed_2

for kk = 1:1
    if spatialDim == 1
        if kk == 1, L = length(tau_adapt_vec);
        elseif kk == 2, L = length(w_adapt_vec);
        elseif kk == 3, L = length(w_facil_vec);
        elseif kk == 4, L = length(tau_facil_vec); 
        end
    else
        L = 10;
    end

for kkk = 4:10
    [kk,kkk]

%parameters
    %simulation duration
    if spatialDim == 1
        T_transition = 0.35;%10;%2
        T_run = 6;%0;%8; %duration of run between stopping periods (sec)
        T_stop = 12; %duration of stopping periods
    else
        T_transition = 1;
        T_run = 6; %duration of run between stopping periods (sec)
        T_stop = 1; %duration of stopping periods
    end
    if use_plot == 1
        numStops = 4; %number of stops
    else
        if spatialDim == 1
            numStops = 400;
        else
            numStops = 200;
        end
    end
    T = numStops*(T_run + T_stop + 2*T_transition); %total simulation time
    % T = numStops*(T_run); %total simulation time
    % keyboard

    %network parameters
        if spatialDim==2, n = 64; %num of neurons per spatialDim
        elseif spatialDim==1, n = 64;
        end
        big = 2*n; %padding for convolutions
        N = n^spatialDim;
        dt = 0.5/1000; %step size
        if spatialDim==1
            m = 1; %CV = 1/sqrt(m)
        else
            m = 4;
        end
    
        if spatialDim==1
            tau_s = 30/1000; %synaptic time constant
            beta_0 = 400; %uniform excitation    
            w_rec = 60; %amplitude of recurrent exc.
            sig_rec = 0.1; %width of recurrent exc.
            w_inh = 67; %global inhibition
            w_adapt = 4; %amplitude of adaptation
            tau_adapt = 3; %adaptation time constant !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            w_FF = 300; %amplitude of feedforward exc.
            sig_FF = 0.01; %width of feedforward exc.
            w_facil = 0; %amplitude of facilitation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tau_facil = 16;
            v_max = 0.1;
            T_theta_max = 0.4; 
            T_replay_min = 0.8;

            if kk == 1, tau_adapt = tau_adapt_vec(kkk);
            elseif kk == 2, w_adapt = w_adapt_vec(kkk);
            elseif kk == 3, w_facil = w_facil_vec(kkk);
            elseif kk == 4, tau_facil = tau_facil_vec(kkk); w_facil = 1;
            end
        else
            tau_s = 20/1000; %synaptic time constant
            beta_0 = 600; %uniform excitation    
            w_rec = 24; %strength of recurrent exc.
            sig_rec = 0.14; %width of recurrent exc.
            w_inh = w_rec; %global inhibition
            w_adapt = 60; %amplitude of adaptation 
            tau_adapt = 0.8; %adaptation time constant 
            w_FF = 400;%2000; %amplitude of feedforward exc.
            sig_FF = 0.03; %0.001; %width of feedforward exc.
            v_max = 0.3;
            T_theta_max = 0.12;

            if kk == 1, v_max = v_max_vec(kkk);
            elseif kk == 2, T_theta_max = T_theta_max_vec(kkk);
            elseif kk == 3, w_adapt = w_adapt_vec(kkk);
            elseif kk == 4, tau_adapt = tau_adapt_vec(kkk);
            elseif kk == 5
                sig_rec = sig_rec*sig_rec_vec(kkk); 
                w_rec = w_rec/sig_rec_vec(kkk); 
                w_inh = w_inh/sig_rec_vec(kkk);
            elseif kk == 6
                sig_FF = sig_FF*sig_FF_vec(kkk); 
                w_FF = w_FF/sig_FF_vec(kkk); 
            end
        end

        T_plot = T_theta_max;%/12;

        %place field centers
        if spatialDim == 2
            [X,Y] = meshgrid((1:n)/n,(1:n)/n); 
        else
            X = (1:n)'/n;
        end

        %activity envelope
        if use_envelope == 1
            kappa = 0.4; %controls width of main body of envelope
            a0 = 40;    %contrls steepness of envelope
            A = zeros(n,1);
            for i = 1:n
                r = abs(i-n/2);
                if r<kappa*n 
                    A(i) = 1;
                else
                    A(i) = exp(-a0*((r-kappa*n)/((1-kappa)*n))^2);
                end
            end
            if spatialDim == 2
                A = A*A';
            end
        else
            if spatialDim == 1
                A = ones(n,1);
            else
                A = ones(n,n);
            end
        end

%rat trajectory
    if spatialDim == 1
        %rat max speed
        % v_max = 0.1;

        v_transition = (1+sin(2*pi*(dt:dt:T_transition)'/2/(T_transition) - pi/2))/2;
        v_lap = v_max*[zeros(T_stop/dt,1); v_transition; ones(T_run/dt,1); v_transition(end:-1:1)];
        v = repmat(v_lap,1,numStops); v = v(:);
            if use_plot == 1, v = circshift(v,-(T_stop + T_transition + T_run/2)/dt); end %uncomment this if sim starts with rat running full speed
        x = mod(cumsum(v)*dt,1);
        acc = [abs(diff(v)); 0];

        % g_transition = 1-acc/max(acc);

        v_thr = 0.03;
        d = 0*(v>v_thr) + 1*(v<v_thr);
        [~,locs] = findpeaks(abs(diff(d)));
        g_transition = zeros(1,T/dt);
        g_transition(locs) = 1;
        g_transition = conv(g_transition,setUp_gaussFilt([1,20/dt],4000/dt),'same');
        g_transition = 1-g_transition/max(g_transition);
        % g_transition = ones(1,T/dt);

        % subplot(211)
        % plot(dt:dt:T,v,'k.','linewidth',2)
        % 
        % subplot(212)
        % plot(dt:dt:T,g_transition,'b','linewidth',2)
        % vline(dt*locs,'g--')
        % keyboard

        % v = v_max*ones(T/dt,1);
        % x = mod(cumsum(v)*dt,1);

        % v = zeros(T/dt,1);
        % x = 0.5*ones(T/dt,1);
    else
        %rat max speed
        % v_max = 0.4; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        % v_lap = v_max*(1+sin(2*pi*(dt:dt:T_run)'/(2*T_run) + pi/2))/2;
        % v_lap = [v_lap; zeros(T_stop/dt,1); v_lap(end:-1:1)];
        % v = repmat(v_lap,1,numStops); v = v(:);
        % x = [0.5*ones(T/dt,1) mod(0.2+cumsum(v)*dt,1)];

        v = v_max*ones(T/dt,1);
        % v = v_max*(1+sin(2*pi*(dt:dt:T)/40))/2;
        x = [0.5*ones(T/dt,1) mod(0.1+cumsum(v)*dt,1)];
        % x = [x(:,2) x(:,1)];

        %correlated random walk
            L = [2+8 98-8]; %size of box (cm)
            % L = [2 98]; %size of box (cm)
            dt_randomWalk = 0.1; %otherwise, takes too long to generate random walk
            speed_sigma = 0; %step size variance
            theta_sigma = 0.02*2*pi; %orientation variance (each step)

            v_transition = (1+sin(2*pi*(dt_randomWalk:dt_randomWalk:T_transition)'/2/(T_transition) - pi/2))/2;

            %fixed max speed
            % v_lap = v_max*[zeros(T_stop/dt_randomWalk,1); v_transition; ones(T_run/dt_randomWalk,1); v_transition(end:-1:1)];
            v_lap = v_max*ones(T_run/dt_randomWalk,1);
            v = repmat(v_lap,1,numStops); v = v(:);
                v = repmat(v_lap,1,numStops); 
                v = v(:);

            %variable max speed
            % v_lap = [zeros(T_stop/dt_randomWalk,1); v_transition; ones(T_run/dt_randomWalk,1); v_transition(end:-1:1)];
            % v = repmat(v_lap,1,numStops); v = v(:);
            %     v_max_repmat = repmat(abs(normrnd(v_max,0.1,1,numStops)),length(v_lap),1);
            %     v = repmat(v_lap,1,numStops).*v_max_repmat; 
            %     v = v(:);

            x = load_randomWalk(L,T,100*v,speed_sigma,theta_sigma,dt_randomWalk,[],spatialDim);
            x(:,2:3) = x(:,2:3)/100;
            v = [diff(x(:,2)), diff(x(:,3))]/mode(diff(x(:,1))); v = [v(1,:); v];
            speed = sqrt(v(:,1).^2 + v(:,2).^2);
            angle = atan2(v(:,2),v(:,1));
            x = [x speed angle]; x = [0 x(1,2:end); x];
            x = compute_dataInterpolation(x,(dt:dt:T)');
            v = x(:,4);
            x = x(:,2:3);

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
    end

    % plot(x)
    % hold on
    % plot(abs(diff(v)/dt))
    % hold off
    % keyboard


%frequency of sequence generation
    %constant theta freq 
        % g_theta = (sign(sin(2*pi*(dt:dt:T)'*(1/T_theta_max)-pi/8-pi/4)+0.4)+1)/2; 
        % g_theta_FF = (sign(sin(2*pi*(dt:dt:T)'*(1/T_theta_max)+pi/2-pi/4)-0.9)+1)/2; 

    % %smooth changes to theta frequence
    %     f_theta = v/v_max*(1/T_theta_max - 1) + T_replay_min; 
    %     f_theta_cumsum = cumsum(f_theta);
        % g_theta = (square(2*pi*dt*f_theta_cumsum,70)+1)/2; 
        % g_theta_FF = (square(2*pi*dt*f_theta_cumsum+0.5,15)+1)/2; 

    %ITI sampling
        %nonstationary distribution: sample theta ITI's based on running speed
        if spatialDim == 1
            time = dt:dt:T;
            sigma_theta = 0.5;
            f_theta = v/max(v)*(1/T_theta_max) + (1-v/max(v))*(1/T_replay_min);
            isi = 1./normrnd(f_theta(1),sigma_theta);
            while isi<0
                isi = 1./normrnd(f_theta(1),sigma_theta);
            end
            t = 1+isi;
            Theta = t;
            while t<T
                % isi = 1./normrnd(interp1(dt:dt:T,f_theta,t),sigma_theta);
                % while isi<0
                %     isi = 1./normrnd(interp1(dt:dt:T,f_theta,t),sigma_theta);
                % end

                ind = round(t)/dt-5/dt+1:round(t)/dt+5/dt;
                    ind(ind<=0) = [];
                    ind(ind>length(time)) = [];
                isi = 1./normrnd(interp1(time(ind),f_theta(ind),t),sigma_theta);
                while isi<0
                    isi = 1./normrnd(interp1(time(ind),f_theta(ind),t),sigma_theta);
                end

                t = t+isi;
                Theta = [Theta; t];
            end
            Theta(isnan(Theta)) = [];

        end

        %stationary distribution
        if spatialDim == 2
            sigma_theta = 1.2;%0.6;
            ITI = 1./normrnd((1/T_theta_max),sigma_theta,ceil(T/T_theta_max),1);
                ITI(ITI<=0) = 1/T_theta_max;
            Theta = cumsum(ITI);

            % histogram(ITI,linspace(0,0.25,20),'FaceColor','k','EdgeColor','none','Normalization','probability')
            % axis square, xlabel('theta period (sec)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
            % xlim([0 0.25])
            % keyboard
        end

        %theta modulation
            if spatialDim == 1
                theta = [Theta(1:end-1) Theta(1:end-1) + 0.6*diff(Theta)];
            else
                theta = [Theta(1:end-1) Theta(1:end-1) + 0.5*diff(Theta)];
            end
            theta_ff = [Theta(1:end-1) - 0.1*diff(Theta) Theta(1:end-1) + 0.1*diff(Theta)];
            % theta_ff = [Theta(1:end-1) - 0.02*diff(Theta) Theta(1:end-1) + 0.02*diff(Theta)];
            % theta_ff = [Theta(1:end-1) - 0.1*diff(Theta) Theta(1:end-1) + 0.2*diff(Theta)];

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
    if spatialDim==2
        w = mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec^2*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w;% - w_inh; 
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

        w = mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec^2*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w - w_inh; 
        plot(w(:,end/2),'k','linewidth',5)
        xlabel('neuron'), ylabel('weight'), axis square, set(gca,'fontsize',16), set(gcf,'color','w')
        
    else
        w = normpdf(X,(n/2+1)/n,sig_rec); w = w/max(w); w = w_rec*w;% - w_inh; 
        if use_periodic == 1
            w_ft=fft(fftshift(w));
        else
            w_ft=fft((w),big);
        end
    end

    % hold on, plot(w(:,end/2)), sum(w(:,end/2))
    % continue


%initialize vectors
    activityCOM = nan(T/dt,spatialDim);
    activitySpread = nan(T/dt,1);
    if spatialDim==2
        r = zeros(n,n); a = zeros(n,n);
    else
        r = zeros(n,1); a = zeros(n,1); f = zeros(n,1);
    end
    spk_count = zeros(N,1);
    if use_plot == 1, spk_mat = zeros(N,T/dt); end

    if use_plot==1 && spatialDim == 1
        spk_mat = nan(n,T/dt);
        r_mat = nan(n,T/dt);
        a_mat = nan(n,T/dt);
        f_mat = nan(n,T/dt);
    end

%simulation
    for t=1:T/dt
        % if mod(t,10*T_plot/dt)==0, t/(T/dt), end
        
        %combined input conductance
            if spatialDim==2
                if use_periodic == 1
                    g_rec = real(ifft2(fft2(r).*w_ft));
                else
                    g_rec = real(ifft2(fft2(r,big,big).*w_ft)); g_rec = g_rec(n/2+1:big-n/2,n/2+1:big-n/2);
                end
                g_FF = reshape(mvnpdf([X(:) Y(:)],x(t,:),sig_FF^2*eye(2)),n,n); g_FF = w_FF*g_FF/max(g_FF(:));
            else
                if use_periodic == 1
                    g_rec = real(ifft(fft(r).*w_ft));
                else
                    g_rec = real(ifft(fft(r,big).*w_ft)); g_rec = g_rec(n/2+1:big-n/2);
                end
                g_FF = normpdf(X,x(t,:),sig_FF); g_FF = w_FF*g_FF/max(g_FF(:));
            end
            if use_feedforward==0, g_FF = 0; end
            g_feedbackInh = -w_inh*sum(r(:));
            g_unifInput = beta_0;
            if spatialDim == 1
                g_adapt = -g_transition(t)*w_adapt*a;
            else
                g_adapt = -w_adapt*a;
            end

            if use_theta==1
                G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF;
                if use_facilitation == 1 && spatialDim == 1
                    g_facil = w_facil*f;
                    G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF + g_rest(t)*g_theta(t)*A.*g_facil;
                    % G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF + g_theta(t)*A.*g_facil;
                    % G = g_transition(t)*(g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh)) + g_theta_FF(t)*g_FF + g_rest(t)*g_theta(t)*A.*g_facil;
                end
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
                    if spatialDim==2
                        hist_data = nanmean(r,i);
                    else
                        hist_data = r;
                    end
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
                if use_plot == 1, spk_mat(:,t) = spk(:); end
                
                %update firing rates
                r = r + spk - r*dt/tau_s;
                
                %update adaptation dynamics
                a = a + spk - a*dt/tau_adapt;

                if use_facilitation == 1 && spatialDim == 1
                    f = f + spk*(1-g_rest(t)) - f*dt/tau_facil; 
                    % f = f + spk - f*dt/tau_facil; 
                end
                
            else
                %update firing rates
                r = r + F*dt - r*dt/tau_s;
                
                %update adaptation dynamics
                a = a + F*dt - a*dt/tau_adapt;
            end

            if use_plot==1 && spatialDim == 1
                spk_mat(:,t) = spk;
                r_mat(:,t) = r;
                a_mat(:,t) = a;
                f_mat(:,t) = f;
            end

            % if use_plot==1 && mod(t,T_plot/dt)==0% && t>1/dt  
            %     neuron = 25;
            % 
            %     set(gcf,'color','w')
            %     subplot(6,1,[1 2]),   
            %         activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread) | [compute_sequenceJumps(activityCOM_sub); 0]>1) = nan;
            %         plot(dt:dt:t*dt,activityCOM_sub(1:t),'k-','linewidth',2)
            %         hold on, plot(dt:dt:t*dt,n*x(1:t)-0.5,'r.','linewidth',4), hold off
            %         hold on, hline(neuron,'g:'), hold off
            %         ylabel('neuron'), 
            %         set(gca,'fontsize',14,'xticklabel',[])
            % 
            %     subplot(613), 
            %         tt = dt:dt:t*dt;
            %         vline_efficient(tt(spk_mat(neuron,1:t)>0),[0 1],[])
            %         ylabel('spikes')
            %         ylabel('\sigma^{spk}(t)')
            %         set(gca,'xticklabel',[],'fontsize',14)
            % 
            %     subplot(614), 
            %         plot(dt:dt:t*dt,r_mat(neuron,1:t),'k','linewidth',2)
            %         ylabel('s(t)')
            %         set(gca,'xticklabel',[],'fontsize',14)
            % 
            %     subplot(615), 
            %         plot(dt:dt:t*dt,a_mat(neuron,1:t),'k','linewidth',2)
            %         hold on, plot(dt:dt:t*dt,nanmean(a_mat(1:n/2,1:t)),'m','linewidth',2), hold off
            %         hold on, plot(dt:dt:t*dt,nanmean(a_mat(n/2+1:n,1:t)),'c','linewidth',2), hold off
            %         ylabel('a(t)')
            %         set(gca,'xticklabel',[],'fontsize',14)
            % 
            %     subplot(616), 
            %         plot(dt:dt:t*dt,f_mat(neuron,1:t),'k','linewidth',2)
            %         hold on, plot(dt:dt:t*dt,nanmean(f_mat(1:n/2,1:t)),'m','linewidth',2), hold off
            %         hold on, plot(dt:dt:t*dt,nanmean(f_mat(n/2+1:n,1:t)),'c','linewidth',2), hold off
            %         ylabel('f(t)'), xlabel('time (sec)')
            %         set(gca,'fontsize',14)
            % 
            %     set(findall(gcf,'Type','axes'),'xlim',[0 t*dt])
            %     % set(findall(gcf,'Type','axes'),'xlim',[44 104])
            %     drawnow
            % end
            % continue
        
        %Plot
        if use_plot==1 &&  mod(t,T_plot/dt)==0% && t>1/dt    
            if spatialDim == 1
                % subplot(331), plot(a)
                subplot(331), plot(r,'k','linewidth',2), ylim([0 3]), vline(n*x(t),'r'), xlabel('neuron'), ylabel('firing rate')
                % subplot(334), plot(a,'k','linewidth',2), vline(n*x(t),'r'), xlabel('neuron'), ylabel('adaptation'), %ylim([0 30])
                %     if w_facil ~=0, hold on, plot(f,'b','linewidth',2), hold off, end
                subplot(337), plot(dt:dt:T,v,'k','linewidth',2), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
                subplot(1,3,[2 3])
                    activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread) | [compute_sequenceJumps(activityCOM_sub); 0]>1) = nan;
                    plot(dt*(1:t),activityCOM_sub(1:t),'k-','linewidth',2)
                    hold on, plot(t*dt,n*x(t)-0.5,'ro','linewidth',2), hold off
                    hold on, plot((1:t)*dt,n*x(1:t)-0.5,'r.','markersize',4), hold off
                    ylim([0 n]), xlabel('time (sec)'), ylabel('neuron')
                set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
            else
                % subplot(331), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('firing rate'), clim([0 3])
                % subplot(334), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('adaptation (-)'), clim([0 30])
                % subplot(337), plot(dt:dt:T,v), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
                subplot(231), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('Firing rate'), clim([0 3])
                subplot(234), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('Adaptation'), clim([0 30])
                subplot(1,3,[2 3]), 
                    activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread(:,1)) | [compute_sequenceJumps(activityCOM_sub); 0]>1,:) = nan;
                    % plot((activityCOM_sub(1:t,1)),(activityCOM_sub(1:t,2)),'.','color',0.8*ones(1,3),'markersize',8)
                    % hold on, cplot((activityCOM_sub(t-1/dt+1:t,1)),(activityCOM_sub(t-1/dt+1:t,2)),1:1/dt,'.','markersize',12), hold off
                    cplot((activityCOM_sub(1:t,1)),(activityCOM_sub(1:t,2)),1:t,'.','markersize',12)
                    hold on, plot((n*x(t,1)-0.5),(n*x(t,2)-0.5),'ko','linewidth',2), hold off
                    hold on, plot((n*x(1:t,1)-0.5),(n*x(1:t,2)-0.5),'k','linewidth',2), hold off
                    axis([0 n 0 n]), set(gca,'xtick',[],'ytick',[]), 
                    axis square, 
                    xlabel('neuron'), ylabel('neuron')
                    %title('Center-of-mass trajectory')
                    % camroll(-90)
                set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')

                % set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
                % xl = xlim; yl = ylim; text(xl(1)+0.01*diff(xl),yl(1),strcat('elapsed time:',{' '},num2str(t*dt,'%0.2f'),' sec'),'HorizontalAlignment','left','VerticalAlignment','top','fontsize',14);
                % writeVideo(V,getframe(gcf));
            end
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
            t = dt*(1:T/dt)';
            if spatialDim == 1, angle = [sign(diff(x)); 0];
            else, angle = [atan2(diff(x(:,2)),diff(x(:,1))); nan];
            end
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
            if spatialDim == 1, id_angleColumns = []; else, id_angleColumns = 5; end
            theta_rat = compute_dataInterpolation(rat,theta_time,id_angleColumns);

            [pks,locs] = findpeaks(abs(diff(g_theta_FF)));
            theta_bounds_FF = [locs(1:2:end-1)+1 locs(2:2:end)+1]; 
            theta_timeBounds_FF = t(theta_bounds_FF);
            theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),id_angleColumns);

        %Rat running periods
            v_thr = 0.03;
            d = 0*(v>v_thr) + 1*(v<v_thr);
            [~,locs] = findpeaks(abs(diff(d)));
            rat_runningTime_start = t(locs(1:2:end));
            rat_runningTime_end = t(locs(2:2:end));
            rat_running_timeBounds = [rat_runningTime_start(1:min(length(rat_runningTime_start),length(rat_runningTime_end))) rat_runningTime_end(1:min(length(rat_runningTime_start),length(rat_runningTime_end)))];
            numRunningPeriods = size(rat_running_timeBounds,1);

            %time since rat running
            theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
            theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);
        
        %Rat stopping periods
            v_thr = 0.01;
            d = 0*(v>v_thr) + 1*(v<v_thr);
            [~,locs] = findpeaks(abs(diff(d)));
            rat_stoppingTime_start = t(locs(2:2:end-1));
            rat_stoppingTime_end = t(locs(3:2:end));
            rat_stopping_timeBounds = [rat_stoppingTime_start rat_stoppingTime_end];
            numStoppingPeriods = size(rat_stopping_timeBounds,1);

            %time since rat stopping
            theta_timeSinceStopping = theta_time - repmat(rat_stopping_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceStopping(theta_timeSinceStopping<0) = nan;
            theta_timeSinceStopping = min(theta_timeSinceStopping,[],2);

            % plot(rat(:,1),rat(:,2))
            % vline(rat_running_timeBounds(:,1),'m--')
            % vline(rat_running_timeBounds(:,2),'b--')            
            % vline(rat_stopping_timeBounds(:,1),'r')
            % vline(rat_stopping_timeBounds(:,2),'g')
            % keyboard

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if use_plot == 0
            data_params(kk,kkk).T = T;
            data_params(kk,kkk).dt = dt;
            data_params(kk,kkk).spatialDim = spatialDim;
            data_params(kk,kkk).n = n;
            data_params(kk,kkk).x = x;
            data_params(kk,kkk).theta_bounds = theta_bounds;
            data_params(kk,kkk).theta_bounds_FF = theta_bounds_FF;
            data_params(kk,kkk).rat_stopping_timeBounds = rat_stopping_timeBounds;
            data_params(kk,kkk).rat_running_timeBounds = rat_running_timeBounds;
            data_params(kk,kkk).activityCOM = activityCOM;
            data_params(kk,kkk).T_transition = T_transition;
            data_params(kk,kkk).T_run = T_run;
            data_params(kk,kkk).T_stop = T_stop;
            data_params(kk,kkk).v_max = v_max;
            data_params(kk,kkk).T_theta_max = T_theta_max;
            data_params(kk,kkk).w_adapt = w_adapt;
            data_params(kk,kkk).tau_adapt = tau_adapt;
            if spatialDim == 1
                data_params(kk,kkk).w_facil = w_facil; 
                data_params(kk,kkk).tau_facil = tau_facil; 
            end

            if spatialDim ==1
                save('data_params_1D_sweep_shortRange.mat','data_params','w_facil_vec','tau_adapt_vec','-v7.3')
            else
                save('data_params_2D_bigEnvironment','data_params','v_max_vec','T_theta_max_vec','w_adapt_vec','tau_adapt_vec','sig_rec_vec','sig_FF_vec','-v7.3')
            end
        end

end
end

% save('data_params_2D_bigEnvironment','data_params','v_max_vec','T_theta_max_vec','w_adapt_vec','tau_adapt_vec','-v7.3')

keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%2D analysis
    
    load data_params_2D_bigEnvironment.mat

    % load data_params_2D_bigEnvironment_2.mat
    % load data_params_2D_bigEnvironment_ratSpeed
    % load data_params_2D_bigEnvironment_ratSpeed_2
    % load data_params_2D_bigEnvironment_variableSpeed

%%

    data_relativeAngleToRat_cell = cell(6,1);
    data_relativeAngleToPrev_cell = cell(6,1);
    data_pairwise_opposite_cell = cell(6,1);
    data_distance_cell = cell(6,1);
    data_relativeAngleToRat_mean = [];
    data_relativeAngleToPrev_mean = [];
    for kk = 1:6 %v_max, T_theta, w_adapt, tau_adapt, w_rec, w_ff
        for kkk = 3:10
        [kk,kkk]
        
        T = data_params(kk,kkk).T;
        dt = data_params(kk,kkk).dt;
        spatialDim = data_params(kk,kkk).spatialDim;
        n = data_params(kk,kkk).n;
        x = data_params(kk,kkk).x;
        theta_bounds = data_params(kk,kkk).theta_bounds;
        theta_bounds_FF = data_params(kk,kkk).theta_bounds_FF;
        rat_stopping_timeBounds = data_params(kk,kkk).rat_stopping_timeBounds;
        rat_running_timeBounds = data_params(kk,kkk).rat_running_timeBounds;
        activityCOM = data_params(kk,kkk).activityCOM;
        T_transition = data_params(kk,kkk).T_transition;
        T_run = data_params(kk,kkk).T_run;
        T_stop = data_params(kk,kkk).T_stop;
        v_max = data_params(kk,kkk).v_max;
        T_theta_max = data_params(kk,kkk).T_theta_max;
        w_adapt = data_params(kk,kkk).w_adapt;
        tau_adapt = data_params(kk,kkk).tau_adapt;

        rat_running_timeBounds = [rat_running_timeBounds(:,1)-0.3 rat_running_timeBounds(:,2)];

        %rat position
            t = dt*(1:T/dt)';
            if spatialDim == 1, angle = sign(diff(x));
            else, angle = [atan2(diff(x(:,2)),diff(x(:,1))); nan];
            end
            v = [diff(x(:,1)), diff(x(:,2))]/mode(diff(t)); v = [v(1,:); v];
            speed = sqrt(v(:,1).^2 + v(:,2).^2);
            rat = [t x speed angle];

        %Theta cycles
            theta_timeBounds = t(theta_bounds);
            theta_time = theta_timeBounds(:,1);%nanmean(theta_timeBounds,2);
            numThetaCycles = size(theta_bounds,1);

        %rat VTE
            [FilterA,FilterB]=butter(2,0.02); 
            HD = unwrap(rat(:,5)); 
            HD_diff = abs(diff(HD)); ind_NaN = find(isnan(HD_diff)); HD_diff(ind_NaN) = HD_diff(ind_NaN-1);
            HD_diff_smooth = filtfilt(FilterA,FilterB,HD_diff);
            zldphi = compute_zscore(HD_diff_smooth);
            % HD_diff_smooth_running = HD_diff_smooth(100*rat(:,4)>20);
            % zldphi = (HD_diff_smooth - nanmean(HD_diff_smooth_running))/nanstd(HD_diff_smooth_running);
            zldphi = [zldphi; zldphi(end)];
            theta_VTE = compute_dataInterpolation([rat(:,1) zldphi],theta_time);
            theta_VTE = theta_VTE(:,2);
            
        %rat properties
            if spatialDim == 1, id_angleColumns = []; else, id_angleColumns = 5; end
            theta_rat = compute_dataInterpolation(rat,theta_time,5);

            theta_timeBounds_FF = t(theta_bounds_FF);
            theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),id_angleColumns);

        %Rat running and stopping periods
            numRunningPeriods = size(rat_running_timeBounds,1);
            numStoppingPeriods = size(rat_stopping_timeBounds,1);

            %time since rat running
            theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
            theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);

            %time since rat stopping
            theta_timeSinceStopping = theta_time - repmat(rat_stopping_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceStopping(theta_timeSinceStopping<0) = nan;
            theta_timeSinceStopping = min(theta_timeSinceStopping,[],2);

            %Straight running bouts with no stopping periods
            [pks,locs] = findpeaks(abs(diff(unwrap(angle))),'minPeakHeight',0.8);
            rat_runningStraight_timeBounds = [t(locs(1:end-1)) t(locs(2:end))];
            vec = ones(size(rat_runningStraight_timeBounds,1),1);
            for i = 1:size(rat_stopping_timeBounds,1)
                ind = find(compute_isNumWithinRange(nanmean(rat_stopping_timeBounds(i,:)),rat_runningStraight_timeBounds));
                vec(ind) = 0;
            end
            rat_runningStraight_timeBounds = rat_runningStraight_timeBounds(vec==1,:);
            rat_runningStraight_timeBounds(rat_runningStraight_timeBounds(:,2)-rat_runningStraight_timeBounds(:,1)<0,:) = [];

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

            ind_thetaCycles = find(100*theta_rat(:,4)>20);
            % ind_thetaCycles = find(~isnan(theta_rat_running(:,2)));
            % [~,theta_rat_running] = compute_dataTemporalConcatenation([theta_rat(:,1) theta_rat(:,1)],rat_runningStraight_timeBounds);

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

            %angular regularity
                edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
                data_relativeAngleToRat = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges_angles,'normalization','probability');
                data_relativeAngleToPrev = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges_angles,'normalization','probability');
                
                data_relativeAngleToRat = smooth(data_relativeAngleToRat,3)'; 
                data_relativeAngleToPrev = smooth(data_relativeAngleToPrev,3)';

            %sequence distance
                % edges_distance = linspace(0,40,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
                edges_distance = 0:2:45; centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
                data_distance = histcounts(100*theta_seqDist(ind_thetaCycles,1),edges_distance,'normalization','probability');

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
            
            %save
                data_relativeAngleToRat_cell{kk} = [data_relativeAngleToRat_cell{kk}; data_relativeAngleToRat];
                data_relativeAngleToPrev_cell{kk} = [data_relativeAngleToPrev_cell{kk}; data_relativeAngleToPrev];
                data_pairwise_opposite_cell{kk} = [data_pairwise_opposite_cell{kk}; data_pairwise_opposite(1,:)];
                data_distance_cell{kk} = [data_distance_cell{kk}; data_distance];


            % set(gcf,'color','w')
            % subplot(121)
            % plot(centers_time,100*data_pairwise_opposite(1,:),'k','linewidth',4)
            % xlabel('time between events (sec)'), ylabel('% alternation')
            % xlim([0 lags(13)*thetaPeriod])
            % ylim([0 100])
            % hline(50,'k--')
            % vline(lags(2:13)*thetaPeriod)
            % axis square, box on
            % 
            % subplot(122)
            % thetaPeriod = 0.126;
            % plot(lags,-xcorr_values,'k-','linewidth',4)
            % xlabel('sweep #'), ylabel('autocorrelation')
            % xlim([0 lags(13)])
            % ylim(0.9*[-1 1])
            % hline(0,'k--'), 
            % vline(lags(2:13))
            % axis square, box on
            % set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16); 
            % 
            % 
            % 
            % tt = (min(theta_time):0.001:max(theta_time))';
            % [~,tt_NaN] = compute_dataTemporalConcatenation([tt tt],theta_timeBounds(theta_seqDir(:,2)>0,:));
            % tt_left = zeros(size(tt));
            % tt_right = zeros(size(tt));
            % tt_left(~isnan(tt_NaN(:,2))) = 1;
            % tt_right(isnan(tt_NaN(:,2))) = 1;
            % 
            % [x,lags] = xcorr(tt_left,tt_right,10000,'unbiased');
            % plot(lags,x)
            % 
            % keyboard

        %plots

            % T = 200;
            % data_ratSpeed_evolution = nan(numRunningPeriods,T);
            % for i = 1:numRunningPeriods
            %     i
            %     positions_running = compute_dataTemporalConcatenation([rat(:,1),rat(:,4)],rat_running_timeBounds(i,:));
            %     TT = min(T,size(positions_running,1));
            %     data_ratSpeed_evolution(i,1:TT) = positions_running(1:TT,2);
            % end
            % plot(1/posSampRate*(1:T),nanmean(data_ratSpeed_evolution))


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

        keyboard



             % %relative angle vs rat speed
        %     [edges_speed,centers_speed] = load_timeBins([0 90],1,10);
        %     edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
        %     a = cool(length(edges_speed));
        %     H1 = []; H2 = [];
        %     H1_data = nan(length(centers_speed),2);
        %     H2_data = nan(length(centers_speed),2);
        %     for i = 1:length(centers_speed)
        %         ind_thetaCycles = find(100*theta_rat(:,4)>edges_speed(i,1) & 100*theta_rat(:,4)<=edges_speed(i,2));
        %         % if length(ind_thetaCycles)<15, continue, end
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
            % subplot(131)
            %     shadedErrorBar(centers_time,100*data_pairwise_opposite(1,:),100*data_pairwise_opposite(2,:),'lineprops',{'k','linewidth',2})
            %     xlabel('time between events (sec)'), ylabel('% alternation'), set(gca,'fontsize',14)
            %     hline(50,'k--'), axis square, box on, axis tight
            %     ylim([0 100])
            % 
            % subplot(132)
            %     yyaxis left, plot(180/pi*centers_angles,data_relativeAngleToPrev,'k','linewidth',4)
            %         ylabel('fraction')
            %     yyaxis right, plot(180/pi*centers_angles,data_relativeAngleToRat,'r','linewidth',4)
            %     xlabel('angle (deg)'), axis square, set(gca,'fontsize',14,'xtick',[0 90 180]), set(gcf,'color','w'), xlim(180*[0 1])
            %     ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; 
            %     legend('prev. sweep','rat heading')
            %     ax = gca;
            %     ax.YAxis(1).Color = 'k';
            %     ax.YAxis(2).Color = 'r';
            %     set(gca,'xtick',0:60:180)
            % 
            % subplot(133)
            %     histogram(100*theta_seqDist(ind_thetaCycles,1),edges_distance,'facecolor','k','edgecolor','none','Normalization','probability');
            %     axis square, xlabel('distance (cm)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
            % 
            %     keyboard



            %plot
                % subplot(131)
                % edges_angles = linspace(0,pi,2*21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
                % H1 = []; H2 = [];
                % for i = 1:size(thetaSequenceNum_evolution,1)
                %     data_sub = data_relativeAngleToRat_evolution_mat(:,thetaSequenceNum_evolution(i,1):thetaSequenceNum_evolution(i,2));
                %     H1 = [H1, histcounts((data_sub(:)),edges_angles)'];
                % 
                %     data_sub = data_relativeAngleToPrev_evolution_mat(:,thetaSequenceNum_evolution(i,1):thetaSequenceNum_evolution(i,2));
                %     H2 = [H2, histcounts((data_sub(:)),edges_angles)'];
                % end
                % h = pcolor(centers_thetaSequenceNum_evolution,180/pi*centers_angles,H1); set(h,'edgecolor','none')
                % xlabel('sweep number'), ylabel({'relative angle to','rat heading (deg)'})
                % axis square, set(gca,'fontsize',14), box on, set(gcf,'color','w'), axis tight
                % % hline(0,'k--')
                % cb = colorbar(); ylabel(cb,'count','FontSize',14,'Rotation',270)
                % 
                % subplot(132)
                % yyaxis left, shadedErrorBar(centers_thetaSequenceNum_evolution,180/pi*nanmean(data_relativeAngleToRat_evolution),180/pi*nanstd(data_relativeAngleToRat_evolution)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
                % ylabel({'relative angle to','rat heading (deg)'})
                % yyaxis right, plot(centers_thetaSequenceNum_evolution,nanmean(data_ratSpeed_evolution),'b','linewidth',2)
                % xlabel('sweep number'), ylabel('rat speed')
                % axis square, set(gca,'fontsize',14), box on, set(gcf,'color','w'),
                % ax = gca;
                % ax.YAxis(1).Color = 'k';
                % ax.YAxis(2).Color = 'b';
                % % ylim([0.15 0.45])
                % 
                % subplot(133)
                % yyaxis left, shadedErrorBar(centers_thetaSequenceNum_evolution,100*nanmean(data_pairwise_opposite_evolution),100*nanstd(data_pairwise_opposite_evolution)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
                % ylabel('% alternation'), 
                % yyaxis right, plot(centers_thetaSequenceNum_evolution,nanmean(data_ratSpeed_evolution),'b','linewidth',2)
                % xlabel('sweep number'), ylabel('rat speed')
                % axis square, set(gca,'fontsize',14), box on, set(gcf,'color','w'),
                % ax = gca;
                % ax.YAxis(1).Color = 'k';
                % ax.YAxis(2).Color = 'b';
                % % ylim([0.15 0.45])
                % 
                % keyboard

            % for i = 1:size(rat_runningStraight_timeBounds,1)%numRunningPeriods%
            % 
            %     Times_snippet = rat_runningStraight_timeBounds(i,:);
            %     % Times_snippet = [3.5 4.6];
            %     % Times_snippet = [23.1 23.8];
            %     % Times_snippet = Times_snippet + [0.5 -0.5]
            % 
            %     [theta_time_snippet,theta_time_snippet_nan] = compute_dataTemporalConcatenation([theta_time theta_time],Times_snippet);
            %         ind_thetaCycles = find(~isnan(theta_time_snippet_nan(:,2)));
            %         theta_time_snippet = theta_time_snippet(:,2); 
            %     x_NaN_snippet_origin_rot_full_sub = [reshape(activityCOM_rot_full(ind_thetaCycles,:,1),[],1) reshape(activityCOM_rot_full(ind_thetaCycles,:,2),[],1)];
            % 
            %     positions_snippet = compute_dataTemporalConcatenation(rat,Times_snippet);
            %     x_NaN_snippet = compute_dataTemporalConcatenation([rat(:,1) activityCOM],Times_snippet);
            % 
            %     subplot(1,3,[1 2])
            %     plot(n*positions_snippet(:,2)-0.5,n*positions_snippet(:,3)-0.5,'k:','linewidth',2)
            %     hold on, cplot(x_NaN_snippet(:,2),x_NaN_snippet(:,3),1:size(x_NaN_snippet,1),'.','markersize',20), hold off
            %     axis square, axis equal, set(gca,'fontsize',14)
            % 
            %     for j = 1:length(ind_thetaCycles)
            %         ii = ind_thetaCycles(j);
            %         hold on, quiver(n*theta_rat(ii,2)-0.5,n*theta_rat(ii,3)-0.5,cos(theta_seqDir(ii,1)),sin(theta_seqDir(ii,1)),10,'r','linewidth',2,'maxHeadSize',0.6), hold off
            %     end
            % 
            %     subplot(233)
            %     plot(theta_time_snippet-Times_snippet(1),180/pi*theta_seqDir(ind_thetaCycles,2),'r','linewidth',3)
            %     hold on, plot(theta_time_snippet-Times_snippet(1),180/pi*theta_seqDir(ind_thetaCycles,3),'k','linewidth',3), hold off
            %     ylabel({'relative angle'}), xlabel('time (sec)'), set(gca,'fontsize',14),
            %     hline(0,'k--'), xlim([0 diff(Times_snippet)]), axis square, %pbaspect([4 2 1]), 
            %     legend('rat heading','prev sweep')
            % 
            %     subplot(236)
            %     plot(-x_NaN_snippet_origin_rot_full_sub(:,2),x_NaN_snippet_origin_rot_full_sub(:,1),'k.','markersize',16)
            %     axis square, set(gca,'fontsize',14), set(gcf,'color','w')
            %     axis([-10 10 -5 10])
            %     vline(0,'k--'), hline(0,'k--')
            %     xlabel('x (cm)'), ylabel('y (cm)')
            % 
            %     keyboard, continue
            % 
            % end


       end
    end

    %%

    save('param_sweep_2D_summary.mat','data_pairwise_opposite_cell','data_relativeAngleToPrev_cell','data_relativeAngleToRat_cell','data_distance_cell')

    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%

    %%

    load('param_sweep_2D_summary.mat')

    v_max_vec = 0.1:0.1:1;
    T_theta_max_vec = 0.08:0.02:0.26;
    w_adapt_vec = 10:10:100;
    tau_adapt_vec = 0.4:0.4:4;
    sig_rec = 1;%0.14;
    sig_rec_vec = sig_rec*(0.4:0.2:2.2);
    sig_FF = 1;%0.03; 
    sig_FF_vec = sig_FF*(0.4:0.2:2.2);

    edges_time = linspace(0,2,60); centers_time = edges_time(1:end-1) + mean(diff(edges_time))/2;
    edges_angles = linspace(0,pi,21); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
    edges_distance = linspace(0,40,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;

    for k = 1:4
        if k == 1, Data = data_distance_cell; xticks = (centers_distance-mean(diff(centers_distance)/2)); clabel = 'fraction';    
        elseif k == 2, Data = data_pairwise_opposite_cell; xticks = centers_time; clabel = '% alternation';
        elseif k == 3, Data = data_relativeAngleToRat_cell; xticks = 180/pi*(centers_angles-mean(diff(centers_angles)/2)); clabel = 'fraction';
        elseif k == 4, Data = data_relativeAngleToPrev_cell; xticks = 180/pi*(centers_angles-mean(diff(centers_angles)/2)); clabel = 'fraction';
        end
     
        kk_vec = [3 4 5 6];
        kk_vec = [5 6];
        KK = 2;
        for kk = 1:KK
            data = Data{kk_vec(kk)};
            data_max = repmat(max(data,[],2),1,size(data,2));
            % data = data./data_max;
            if k==2, data = 100*data; end

            if kk_vec(kk) == 1, yticks = v_max_vec; 
            elseif kk_vec(kk) == 2, yticks = T_theta_max_vec; 
            elseif kk_vec(kk) == 3, yticks = w_adapt_vec;
            elseif kk_vec(kk) == 4, yticks = tau_adapt_vec;
            elseif kk_vec(kk) == 5, yticks = sig_rec_vec; 
            elseif kk_vec(kk) == 6, yticks = sig_FF_vec; 
            end 

            subplot(KK,4,k + 4*(kk-1))
            % subplot(2,4,k + 4*(kk-1-4))
            % subplot(1,4,k)
            h = pcolor(xticks,yticks,data); set(h,'edgecolor','none')
            axis square, set(gcf,'Color','w'), set(gca,'fontsize',14,'ytick',yticks(2:2:end) + mean(diff(yticks))/2,'yticklabel',yticks(2:2:end))
            colorbar

            if k == 1 
                caxis([0 0.4])
            elseif k == 2
                caxis([30 70])
            elseif k == 3
                caxis([0 0.18]), 
                set(gca,'xtick',0:60:180)
                vline(30,'k:')
            elseif k == 4
                caxis([0 0.2]), 
                set(gca,'xtick',0:60:180)
                vline(60,'k:')
            end   
                
            if kk_vec(kk) == 1, hline(0.3 + mean(diff(yticks))/2,'r--')
            elseif kk_vec(kk) == 2, hline(0.12 + mean(diff(yticks))/2,'r--')
            elseif kk_vec(kk) == 3, hline(40 + mean(diff(yticks))/2,'r--')
            elseif kk_vec(kk) == 4, hline(0.8 + mean(diff(yticks))/2,'r--')
            elseif kk_vec(kk) == 5, hline(sig_rec + mean(diff(yticks))/2,'r--')
            elseif kk_vec(kk) == 6, hline(sig_FF + mean(diff(yticks))/2,'r--')
            end 

            if k == 1
                if kk_vec(kk) == 1, ylabel('v_{max}')
                elseif kk_vec(kk) == 2, ylabel('T_{\Theta}')
                elseif kk_vec(kk) == 3, ylabel('w_{adapt}')
                elseif kk_vec(kk) == 4, ylabel('\tau_{adapt}')
                elseif kk_vec(kk) == 5, ylabel('\sigma_{rec}')
                elseif kk_vec(kk) == 6, ylabel('\sigma_{ff}')
                end 
            end

            if kk == KK
                if k == 1 
                    xlabel('distance (cm)')
                elseif k == 2
                    xlabel({'time between events (sec)'})
                elseif k == 3
                    xlabel({'relative angle','to rat (deg)'})
                elseif k == 4
                    xlabel({'relative angle','to prev (deg)'})
                end   
            end

            cb = colorbar(); ylabel(cb,clabel,'FontSize',14,'Rotation',270)
        end
    end

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16); 





%%
    % %distance as a function of adaptation strength
    %     kk = 3; %v_max, T_theta, w_adapt, tau_adapt
    %     data = nan(2,10);
    %     for kkk = 1:10
    %     [kk,kkk]
    % 
    %     T = data_params(kk,kkk).T;
    %     dt = data_params(kk,kkk).dt;
    %     spatialDim = data_params(kk,kkk).spatialDim;
    %     n = data_params(kk,kkk).n;
    %     x = data_params(kk,kkk).x;
    %     theta_bounds = data_params(kk,kkk).theta_bounds;
    %     rat_stopping_timeBounds = data_params(kk,kkk).rat_stopping_timeBounds;
    %     rat_running_timeBounds = data_params(kk,kkk).rat_running_timeBounds;
    %     activityCOM = data_params(kk,kkk).activityCOM;
    % 
    %     %rat position
    %         t = dt*(1:T/dt)';
    %         if spatialDim == 1, angle = sign(diff(x));
    %         else, angle = [atan2(diff(x(:,2)),diff(x(:,1))); nan];
    %         end
    %         v = [diff(x(:,1)), diff(x(:,2))]/mode(diff(t)); v = [v(1,:); v];
    %         speed = sqrt(v(:,1).^2 + v(:,2).^2);
    %         rat = [t x speed angle];
    % 
    %     %Theta cycles
    %         theta_timeBounds = t(theta_bounds);
    %         theta_time = nanmean(theta_timeBounds,2);
    %         numThetaCycles = size(theta_bounds,1);
    % 
    %     %rat properties
    %         if spatialDim == 1, id_angleColumns = []; else, id_angleColumns = 5; end
    %         theta_rat = compute_dataInterpolation(rat,theta_time,5);
    % 
    %         %Straight running bouts with no stopping periods
    %         [pks,locs] = findpeaks(abs(diff(unwrap(angle))),'minPeakHeight',0.8);
    %         rat_runningStraight_timeBounds = [t(locs(1:end-1)) t(locs(2:end))];
    %         vec = ones(size(rat_runningStraight_timeBounds,1),1);
    %         for i = 1:size(rat_stopping_timeBounds,1)
    %             ind = find(compute_isNumWithinRange(nanmean(rat_stopping_timeBounds(i,:)),rat_runningStraight_timeBounds));
    %             vec(ind) = 0;
    %         end
    %         rat_runningStraight_timeBounds = rat_runningStraight_timeBounds(vec==1,:);
    %         rat_runningStraight_timeBounds(rat_runningStraight_timeBounds(:,2)-rat_runningStraight_timeBounds(:,1)<2,:) = [];
    % 
    %     %sequence direction
    %         theta_seqDist = nan(numThetaCycles,1);
    %         for i = 2:numThetaCycles
    %             theta_sub = activityCOM(theta_bounds(i,1):theta_bounds(i,2),:);
    %             theta_seqDist(i) = compute_sequenceDistance(unwrap(theta_sub))/n;
    %         end
    % 
    %         [~,theta_rat_running] = compute_dataTemporalConcatenation([theta_rat(:,1) theta_rat(:,1)],rat_runningStraight_timeBounds);
    %         ind_thetaCycles = find(~isnan(theta_rat_running(:,2)));
    % 
    %         data(1,kkk) = nanmean(100*theta_seqDist(ind_thetaCycles));
    %         data(2,kkk) = nanstd(100*theta_seqDist(ind_thetaCycles))./sqrt(length(ind_thetaCycles));
    %     end
    % 
    %     w_adapt_vec = 10:10:100;
    %     plot(w_adapt_vec,data(1,:),'k','linewidth',2)
    %     axis square, axis tight, xlabel('w_{adapt}'), ylabel('distance (cm)'), set(gca,'fontsize',14), set(gcf,'color','w'), box on

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sudo find /media/john/WorkingData/ -type f -exec chmod 666 {} \;


%1D analysis

    % load data_params_1D
    load data_params_1D_shortAdaptationTimeConstant
    % load data_params_1D_alwaysOn

    load data_params_1D_sweep.mat
    
    n = 64;
    data_pairwise_opposite_cell = cell(2,10);
    data_difference_cell = cell(2,10);
    for kk = 2:2
        L = length(data_params(kk,:));
        for kkk = 1:L
            [kk,kkk]
    
            if isempty(data_params(kk,kkk).T), continue, end
    
            T = data_params(kk,kkk).T;
            dt = data_params(kk,kkk).dt;
            spatialDim = data_params(kk,kkk).spatialDim;
            n = data_params(kk,kkk).n;
            x = data_params(kk,kkk).x;
            theta_bounds = data_params(kk,kkk).theta_bounds;
            theta_bounds_FF = data_params(kk,kkk).theta_bounds_FF;
            rat_stopping_timeBounds = data_params(kk,kkk).rat_stopping_timeBounds;
            rat_running_timeBounds = data_params(kk,kkk).rat_running_timeBounds;
            activityCOM = data_params(kk,kkk).activityCOM;
            T_transition = data_params(kk,kkk).T_transition;
            T_run = data_params(kk,kkk).T_run;
            T_stop = data_params(kk,kkk).T_stop;
            v_max = data_params(kk,kkk).v_max;
            T_theta_max = data_params(kk,kkk).T_theta_max;
            w_adapt = data_params(kk,kkk).w_adapt;
            tau_adapt = data_params(kk,kkk).tau_adapt;
            w_facil = data_params(kk,kkk).w_facil;
            tau_facil = data_params(kk,kkk).tau_facil;
        
            %rat position
                t = dt*(1:T/dt)';
                angle = sign(diff(x)); angle = [angle(1,:); angle];
                speed = diff(x(:,1))/mode(diff(t)); speed = [speed(1,:); speed];
                rat = [t x speed angle];
        
            %Theta cycles
                theta_timeBounds = t(theta_bounds);
                theta_time = theta_timeBounds(:,1);%nanmean(theta_timeBounds,2);
                numThetaCycles = size(theta_bounds,1);
        
            %rat properties
                theta_rat = compute_dataInterpolation(rat,theta_time,[]);
        
                theta_timeBounds_FF = t(theta_bounds_FF);
                theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),[]);
        
            %Rat running and stopping periods
                numRunningPeriods = size(rat_running_timeBounds,1);
                numStoppingPeriods = size(rat_stopping_timeBounds,1);
        
                %time since rat running
                theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
                theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
                theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);
        
                %time since rat stopping
                theta_timeSinceStopping = theta_time - repmat(rat_stopping_timeBounds(:,1)',numThetaCycles,1);
                theta_timeSinceStopping(theta_timeSinceStopping<0) = nan;
                theta_timeSinceStopping = min(theta_timeSinceStopping,[],2);
        
            %sequence direction
                theta_seqDir = nan(numThetaCycles,1);
                theta_seqDist = nan(numThetaCycles,1);
                for i = 2:numThetaCycles
                    theta_sub = unwrap(activityCOM(theta_bounds(i,1):theta_bounds(i,2),:));
                    
                    delta = diff(theta_sub,[],1);
                    delta(prod(delta,2)==0,:) = [];
                    theta_seqDir(i,1) = sign(nanmean(delta));
                    theta_seqDist(i) = abs(theta_sub(end) - theta_sub(1));
        
                    hold on, plot(t(theta_bounds(i,1):theta_bounds(i,2)),theta_sub-(n*rat(theta_bounds(i,1),2)-0.5))
                    % title(theta_seqDist(i))
                    title(theta_seqDir(i,1))
                    drawnow, %keyboard
                end

            keyboard
        
            theta_seqDist_thr = 3;
        
            %forward vs reverse rate
                [edges_time,centers_time] = load_timeBins([0 12],0.5,0.5);
                edges_distance = linspace(0.15,0.46,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
                data_distance = nan(numStoppingPeriods,length(edges_distance)-1);
                data_forward = nan(numStoppingPeriods,length(edges_time)-1);
                data_reverse = nan(numStoppingPeriods,length(edges_time)-1);
                data_difference = nan(numStoppingPeriods,length(edges_time)-1);
                data_eventRate = nan(numStoppingPeriods,length(edges_time)-1);
                data_length = nan(numStoppingPeriods,length(edges_time)-1);
                for i = 1:numStoppingPeriods
                    data_sub = compute_dataTemporalConcatenation([theta_time theta_timeSinceStopping theta_seqDir theta_rat(:,4) theta_seqDist],rat_stopping_timeBounds(i,:));
                    data_sub(data_sub(:,end)<theta_seqDist_thr,:) = [];
                    if size(data_sub,1)<2, continue, end
        
                    %distance
                        data_distance(i,:) = histcounts(data_sub(:,5),edges_distance,'normalization','pdf');
        
                    for j = 1:length(centers_time)
                        %forward vs. reverse rate
                        ind = find(data_sub(:,2)>edges_time(j,1) & data_sub(:,2)<=edges_time(j,2));
                        if length(ind)>0
                            data_forward(i,j) = sum(data_sub(ind,3)==1)/mean(diff(edges_time,[],2));
                            data_reverse(i,j) = sum(data_sub(ind,3)==-1)/mean(diff(edges_time,[],2));
                            data_difference(i,j) = data_forward(i,j) - data_reverse(i,j);
                            data_length(i,j) = nanmean(data_sub(ind,5));
                            data_eventRate(i,j) = length(ind);
                        else
                            data_forward(i,j) = 0;
                            data_reverse(i,j) = 0;
                            data_eventRate(i,j) = 0;
                            data_difference(i,j) = 0;
                        end
                    end
        
                    data_forward(i,:) = conv(data_forward(i,:),setUp_gaussFilt([1 10],1),'same');
                    data_reverse(i,:) = conv(data_reverse(i,:),setUp_gaussFilt([1 10],1),'same');
                    data_difference(i,:) = conv(data_difference(i,:),setUp_gaussFilt([1 10],1),'same');
        
                end
        
                % p_forVsRev = nan(1,length(centers_time));
                % for j = 1:length(centers_time)
                %     [~,p_forVsRev(j)] = ttest(data_forward(:,j),data_reverse(:,j));
                % end
        
            %boostraps and shuffles
                [~,theta_time_nan] = compute_dataTemporalConcatenation([theta_time theta_time],rat_stopping_timeBounds);
                ind_thetaCycles = find(theta_seqDist>theta_seqDist_thr & ~isnan(theta_time_nan(:,2)));
                replayEvents = [theta_time(ind_thetaCycles) theta_seqDir(ind_thetaCycles)];
                [edges_time_pairwise,centers_time_pairwise] = load_timeBins([0 10],0.1,0.4);
        
                %mean
                delT = squareform(pdist(replayEvents(:,1))); delT = delT(:);
                delX = squareform(pdist(replayEvents(:,2))); delX = delX(:)/2;
                data_opposite = nan(1,length(centers_time_pairwise));
        
                %bootstraps
                numBoostraps = 500;
                data_bootstraps = nan(numBoostraps,length(centers_time_pairwise));
                for j = 1:length(centers_time_pairwise)
                    ind = find(delT>edges_time_pairwise(j,1) & delT<=edges_time_pairwise(j,2));
                    data_opposite(j) = nanmean(delX(ind));
        
                    for i = 1:numBoostraps
                        data_bootstraps(i,j) = nanmean(datasample(delX(ind),length(ind)));
                    end
                end
                bootstraps_margins = quantile(data_bootstraps,[0.0250 0.975]);
        
                % %shuffles
                % numShuffles = 500;
                % data_shuffles = nan(numShuffles,length(centers_time_pairwise));
                % for i = 1:numShuffles
                %     i
                %     replayEvents_shuffle = replayEvents;
                %     replayEvents_shuffle(:,1) = replayEvents(randperm(size(replayEvents,1))',1);
                %     delT = squareform(pdist(replayEvents_shuffle(:,1))); delT = delT(:);
                %     delX = squareform(pdist(replayEvents_shuffle(:,2))); delX = delX(:)/2;
                % 
                %     delX(delT>max(edges_time_pairwise(:))) = [];
                %     delT(delT>max(edges_time_pairwise(:))) = [];
                %     for j = 1:length(centers_time_pairwise)
                %         ind = find(delT>edges_time_pairwise(j,1) & delT<=edges_time_pairwise(j,2));
                %         data_shuffles(i,j) = nanmean(delX(ind));
                %     end
                % end
                % shuffles_margins = quantile(data_shuffles,[0.0250 0.975]);
        
            % %plot
            %     figure
            %     colors = [.4660 0.6740 0.1880;0.4940 0.1840 0.5560];
            %     subplot(131)
            %     shadedErrorBar(centers_time,nanmean(data_forward),nanstd(data_forward)./sqrt(size(data_forward,1)),'lineprops',{'color',colors(1,:),'linewidth',2})
            %     hold on, shadedErrorBar(centers_time,nanmean(data_reverse),nanstd(data_reverse)./sqrt(size(data_reverse,1)),'lineprops',{'color',colors(2,:),'linewidth',2}), hold off
            %     xlabel('Time since arrival (s)'), ylabel('Events/sec')
            %     axis square, set(gca,'fontsize',14), set(gcf,'color','w'), box on, axis tight
            %     hold on, plot(centers_time(p_forVsRev<0.05),0.7,'k.','markersize',10), hold off, legend('forward','reverse')
            %     xlim([0 10])
            %     ylim([0 0.8])
            % 
            %     subplot(132)
            %     shadedErrorBar(centers_time,nanmean(data_difference),nanstd(data_difference)./sqrt(size(data_difference,1)),'lineprops',{'k','linewidth',2})
            %     xlabel('Time since arrival (s)'), ylabel('Events/sec'),
            %     axis square, set(gca,'fontsize',14), set(gcf,'color','w'), box on, axis tight
            %     hline(0,'k')
            %     hold on, plot(centers_time(p_forVsRev<0.05),0.4,'k.','markersize',10), hold off
            %     xlim([0 10])
            %     ylim([-0.5 0.5])
            % 
            %     subplot(133)
            %     shadedErrorBar(centers_time_pairwise,100*data_opposite,100*abs(bootstraps_margins-data_opposite),'lineprops',{'k','linewidth',2})
            %     % hold on, plot(centers_time_pairwise,100*shuffles_margins,'k--'), hold off
            %     xlabel('time between events (sec)'), ylabel('% opposite'), set(gca,'fontsize',14)
            %     axis square, box on, axis tight
            %     axis([0 10 0 100])
            %     hold on, plot(centers_time_pairwise(data_mean>shuffles_margins(2,:) | data_mean<shuffles_margins(1,:)),90,'k.','markersize',10), hold off
    

            data_pairwise_opposite_cell{kk,kkk} = data_opposite;
            data_difference_cell{kk,kkk} = nanmean(data_difference);
        end
    end

    tau_adapt_vec = 0.4:0.4:4; %[0.8,3];
    w_adapt_vec = 2:2:16;

    for kk = 1:2
        data_pairwise_opposite = [];
        data_difference = [];
        colors = cool(10);

        if kk == 1 
            param_vec = string(tau_adapt_vec);
            param_vec(1) = strcat('tau_a =',param_vec(1));
        else
            param_vec = string(w_adapt_vec);
            param_vec(1) = strcat('w_a =',param_vec(1));
        end

        for kkk = 1:10
            if isempty(data_difference_cell{kk,kkk}), continue, end

            subplot(2,2,2*(kk-1)+1)
            if kkk==1, plot(centers_time,data_difference_cell{kk,kkk},'color',colors(kkk,:),'linewidth',2);
            else, hold on, plot(centers_time,data_difference_cell{kk,kkk},'color',colors(kkk,:),'linewidth',2); hold off
            end
            xlabel('time since stopping (sec)'), ylabel('rate difference (event/sec)'), set(gca,'fontsize',14), set(gcf,'color','w')
            legend(param_vec)
            
            subplot(2,2,2*(kk-1)+2)
            if kkk==1, plot(centers_time_pairwise,data_pairwise_opposite_cell{kk,kkk},'color',colors(kkk,:),'linewidth',2);
            else, hold on, plot(centers_time_pairwise,data_pairwise_opposite_cell{kk,kkk},'color',colors(kkk,:),'linewidth',2); hold off
            end
            xlabel('time difference (sec)'), ylabel('% opposite'), set(gca,'fontsize',14), 
            
        end
    end














    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % 
    % cwd = pwd;
    % cd /media/john/WorkingData/MATLAB/replay_attractorDynamics/thetaSequences/
    % load data_params_2D_bigEnvironment.mat
    % cd(cwd)
    % 
    % data_relativeAngleToRat_cell = cell(4,10);
    % data_relativeAngleToPrev_cell = cell(4,10);
    % data_pairwise_opposite_cell = cell(4,10);
    % data_relativeAngleToRat_evolution_cell = cell(4,10);
    % data_relativeAngleToPrev_evolution_cell = cell(4,10);
    % data_pairwise_opposite_evolution_cell = cell(4,10);
    % data_distance_cell = cell(4,10);
    % data_relativeAngleToRat_meanVariance_oneSided_cell = cell(4,10);
    % data_relativeAngleToPrev_meanVariance_oneSided_cell = cell(4,10);
    % for kk = 1:4 %v_max, T_theta, w_adapt, tau_adapt
    %     for kkk = 1:10
    %     [kk,kkk]
    % 
    %     T = data_params(kk,kkk).T;
    %     dt = data_params(kk,kkk).dt;
    %     spatialDim = data_params(kk,kkk).spatialDim;
    %     n = data_params(kk,kkk).n;
    %     x = data_params(kk,kkk).x;
    %     theta_bounds = data_params(kk,kkk).theta_bounds;
    %     theta_bounds_FF = data_params(kk,kkk).theta_bounds_FF;
    %     rat_stopping_timeBounds = data_params(kk,kkk).rat_stopping_timeBounds;
    %     rat_running_timeBounds = data_params(kk,kkk).rat_running_timeBounds;
    %     activityCOM = data_params(kk,kkk).activityCOM;
    %     T_transition = data_params(kk,kkk).T_transition;
    %     T_run = data_params(kk,kkk).T_run;
    %     T_stop = data_params(kk,kkk).T_stop;
    %     v_max = data_params(kk,kkk).v_max;
    %     T_theta_max = data_params(kk,kkk).T_theta_max;
    %     w_adapt = data_params(kk,kkk).w_adapt;
    %     tau_adapt = data_params(kk,kkk).tau_adapt;
    % 
    %     %rat position
    %         t = dt*(1:T/dt)';
    %         if spatialDim == 1, angle = sign(diff(x));
    %         else, angle = [atan2(diff(x(:,2)),diff(x(:,1))); nan];
    %         end
    %         v = [diff(x(:,1)), diff(x(:,2))]/mode(diff(t)); v = [v(1,:); v];
    %         speed = sqrt(v(:,1).^2 + v(:,2).^2);
    %         rat = [t x speed angle];
    % 
    %     %Theta cycles
    %         theta_timeBounds = t(theta_bounds);
    %         theta_time = nanmean(theta_timeBounds,2);
    %         numThetaCycles = size(theta_bounds,1);
    % 
    %     %rat properties
    %         if spatialDim == 1, id_angleColumns = []; else, id_angleColumns = 5; end
    %         theta_rat = compute_dataInterpolation(rat,theta_time,5);
    % 
    %         theta_timeBounds_FF = t(theta_bounds_FF);
    %         theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),id_angleColumns);
    % 
    %     %Rat running and stopping periods
    %         numRunningPeriods = size(rat_running_timeBounds,1);
    %         numStoppingPeriods = size(rat_stopping_timeBounds,1);
    % 
    %         %time since rat running
    %         theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
    %         theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
    %         theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);
    % 
    %         %time since rat stopping
    %         theta_timeSinceStopping = theta_time - repmat(rat_stopping_timeBounds(:,1)',numThetaCycles,1);
    %         theta_timeSinceStopping(theta_timeSinceStopping<0) = nan;
    %         theta_timeSinceStopping = min(theta_timeSinceStopping,[],2);
    % 
    %         %Straight running bouts with no stopping periods
    %         [pks,locs] = findpeaks(abs(diff(unwrap(angle))),'minPeakHeight',0.8);
    %         rat_runningStraight_timeBounds = [t(locs(1:end-1)) t(locs(2:end))];
    %         vec = ones(size(rat_runningStraight_timeBounds,1),1);
    %         for i = 1:size(rat_stopping_timeBounds,1)
    %             ind = find(compute_isNumWithinRange(nanmean(rat_stopping_timeBounds(i,:)),rat_runningStraight_timeBounds));
    %             vec(ind) = 0;
    %         end
    %         rat_runningStraight_timeBounds = rat_runningStraight_timeBounds(vec==1,:);
    %         rat_runningStraight_timeBounds(rat_runningStraight_timeBounds(:,2)-rat_runningStraight_timeBounds(:,1)<1.5,:) = [];
    % 
    %     %sequence direction
    %         theta_seqDir = nan(numThetaCycles,3);
    %         theta_seqDist = nan(numThetaCycles,1);
    %         activityCOM_rot_full = nan(numThetaCycles,400,2);
    %         activityCOM_rot = nan(size(activityCOM));
    %         for i = 2:numThetaCycles
    %             theta_sub = activityCOM(theta_bounds(i,1):theta_bounds(i,2),:);
    % 
    %             distanceToRat_snippet = sqrt((theta_sub(:,1)-n*theta_rat(i,2)-0.5).^2 + (theta_sub(:,2)-n*theta_rat(i,3)-0.5).^2 );
    % 
    %             %method 1
    %                 % delta = diff(theta_sub,[],1);
    %                 %     delta(prod(delta,2)==0,:) = [];
    %                 % theta_seqDir(i,1) = circ_mean(atan2(delta(:,2),delta(:,1)));
    %                 % theta_seqDir(i,2) = circ_dist(theta_rat(i,5),theta_seqDir(i,1));
    %                 % theta_seqDir(i,3) = circ_dist(theta_seqDir(i-1,1),theta_seqDir(i,1));
    % 
    %             %method 2
    %                 ind_max = find(distanceToRat_snippet==max(distanceToRat_snippet)); ind_max = ind_max(1);
    %                 ind_min = find(distanceToRat_snippet==min(distanceToRat_snippet)); ind_min = ind_min(1);
    %                 delta = diff([n*theta_rat(i,2:3)-0.5; theta_sub(ind_max,:)],[],1);
    %                 % delta = diff([positions_snippet_cut(ind_max,2:3); x_NaN_snippet_cut(ind_max,2:3)],[],1);
    %                 % delta = diff([x_NaN_snippet_cut(ind_min,2:3); x_NaN_snippet_cut(ind_max,2:3)],[],1);
    %                 theta_seqDir(i,1) = atan2(delta(:,2),delta(:,1));
    %                 theta_seqDir(i,2) = circ_dist(theta_rat(i,5),theta_seqDir(i,1));
    %                 theta_seqDir(i,3) = circ_dist(theta_seqDir(i-1,1),theta_seqDir(i,1));
    % 
    %             %sequence distance
    %                 % theta_seqDist(i) = compute_sequenceDistance(unwrap(theta_sub))/n;
    %                 theta_sub_unwrapped = unwrap(theta_sub);
    %                 theta_seqDist(i) = sqrt((theta_sub_unwrapped(1,1)-theta_sub_unwrapped(end,1)).^2 + (theta_sub_unwrapped(1,2)-theta_sub_unwrapped(end,2)).^2);
    % 
    %             %rotate theta sequence by rat heading direction
    %                 theta_sub_origin = theta_sub - repmat(n*theta_rat(i,2:3)-0.5,size(theta_sub,1),1);
    %                 M = [cos(-theta_rat(i,5)) -sin(-theta_rat(i,5)); sin(-theta_rat(i,5)) cos(-theta_rat(i,5))];
    %                 theta_sub_origin_rot = (M*theta_sub_origin')';
    % 
    %                 activityCOM_rot(theta_bounds(i,1):theta_bounds(i,2),:) = theta_sub_origin_rot;
    %                 activityCOM_rot_full(i,:,1) = compute_vecBuffered(theta_sub_origin_rot(:,1),400);
    %                 activityCOM_rot_full(i,:,2) = compute_vecBuffered(theta_sub_origin_rot(:,2),400);
    %         end
    % 
    %     %angle of sequences relative to rat heading direction across run periods
    %         edges_angles = linspace(0,pi,20); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
    %         edges_time_pairwise = linspace(0,2,40); centers_time_pairwise = edges_time_pairwise(1:end-1) + mean(diff(edges_time_pairwise))/2;
    %         edges_distance = linspace(0,20,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
    %         thetaSequenceNum_evolution = load_timeBins([1 30],1,3); 
    %         data_relativeAngleToRat = nan(numRunningPeriods,length(edges_angles)-1);
    %         data_relativeAngleToPrev = nan(numRunningPeriods,length(edges_angles)-1);
    %         data_distance = nan(numRunningPeriods,length(edges_distance)-1);
    %         data_pairwise_opposite = nan(numRunningPeriods,length(edges_time_pairwise)-1);
    %         data_relativeAngleToRat_evolution = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_relativeAngleToRat_evolution_mat = nan(numRunningPeriods,50);
    %         data_relativeAngleToPrev_evolution = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_relativeAngleToPrev_evolution_mat = nan(numRunningPeriods,50);
    %         data_ratSpeed_evolution = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_pairwise_opposite_evolution  = nan(numRunningPeriods,length(thetaSequenceNum_evolution));
    %         data_relativeAngleToRat_meanVariance_oneSided = nan(numRunningPeriods,2);
    %         data_relativeAngleToPrev_meanVariance_oneSided = nan(numRunningPeriods,2);
    %         for i = 1:size(rat_runningStraight_timeBounds,1)%numRunningPeriods%
    % 
    %             Times_snippet = rat_runningStraight_timeBounds(i,:);
    %             % Times_snippet = [3.5 4.6];
    % 
    %             [theta_time_snippet,theta_time_snippet_nan] = compute_dataTemporalConcatenation([theta_time theta_time],Times_snippet);
    %                 ind_thetaCycles = find(~isnan(theta_time_snippet_nan(:,2)));
    %                 theta_time_snippet = theta_time_snippet(:,2); 
    %             % x_NaN_theta_snippet = compute_dataTemporalConcatenation(x_NaN_theta,Times_snippet);
    %             % x_NaN_slow_theta_snippet = compute_dataTemporalConcatenation(x_NaN_slow_theta,Times_snippet);
    %             x_NaN_snippet_origin_rot_full_sub = [reshape(activityCOM_rot_full(ind_thetaCycles,:,1),[],1) reshape(activityCOM_rot_full(ind_thetaCycles,:,2),[],1)];
    % 
    %                 %plot
    %                     % positions_snippet = compute_dataTemporalConcatenation(rat,Times_snippet);
    %                     % x_NaN_snippet = compute_dataTemporalConcatenation([rat(:,1) activityCOM],Times_snippet);
    % 
    %                     % subplot(1,3,[1 2])
    %                     % plot(n*positions_snippet(:,2)-0.5,n*positions_snippet(:,3)-0.5,'k:','linewidth',2)
    %                     % hold on, cplot(x_NaN_snippet(:,2),x_NaN_snippet(:,3),1:size(x_NaN_snippet,1),'.','markersize',20), hold off
    %                     % axis square, axis equal, set(gca,'fontsize',14)
    %                     % 
    %                     % for j = 1:length(ind_thetaCycles)
    %                     %     ii = ind_thetaCycles(j);
    %                     %     hold on, quiver(n*theta_rat(ii,2)-0.5,n*theta_rat(ii,3)-0.5,cos(theta_seqDir(ii,4)),sin(theta_seqDir(ii,4)),10,'r','linewidth',2,'maxHeadSize',0.6), hold off
    %                     % end
    %                     % 
    %                     % 
    %                     % subplot(233)
    %                     % plot(theta_time_snippet-Times_snippet(1),180/pi*theta_seqDir(ind_thetaCycles,2),'r','linewidth',3)
    %                     % hold on, plot(theta_time_snippet-Times_snippet(1),180/pi*theta_seqDir(ind_thetaCycles,3),'k','linewidth',3), hold off
    %                     % ylabel({'relative angle'}), xlabel('time (sec)'), set(gca,'fontsize',14),
    %                     % hline(0,'k--'), xlim([0 diff(Times_snippet)]), axis square, %pbaspect([4 2 1]), 
    %                     % legend('rat heading','prev sweep')
    %                     % 
    %                     % subplot(236)
    %                     % plot(-x_NaN_snippet_origin_rot_full_sub(:,2),x_NaN_snippet_origin_rot_full_sub(:,1),'k.','markersize',16)
    %                     % axis square, set(gca,'fontsize',14), set(gcf,'color','w')
    %                     % vline(0,'k--'), hline(0,'k--')
    %                     % xlabel('x (cm)'), ylabel('y (cm)')
    %                     % keyboard, continue
    % 
    % 
    %             data_sub = compute_dataTemporalConcatenation([theta_time theta_timeSinceRunning theta_seqDir theta_rat(:,4) theta_seqDist],Times_snippet);
    %             if size(data_sub,1)<2, continue, end
    % 
    %             %Angle relative to rat
    %                 data_relativeAngleToRat(i,:) = histcounts(abs(data_sub(:,4)),edges_angles,'normalization','probability');
    %                 data_relativeAngleToPrev(i,:) = histcounts(abs(data_sub(:,5)),edges_angles,'normalization','probability');
    % 
    %             %distance
    %                 data_distance(i,:) = histcounts(data_sub(:,end),edges_distance,'normalization','pdf');
    % 
    %             %percent opposite
    %                 delT = squareform(pdist(data_sub(:,1))); delT = delT(:);
    %                 delX = squareform(pdist(sign(data_sub(:,5)))); delX = delX(:)/2;
    %                 for j = 1:length(edges_time_pairwise)-1
    %                     ind = find(delT>edges_time_pairwise(j) & delT<=edges_time_pairwise(j+1));
    %                     data_pairwise_opposite(i,j) = nanmean(delX(ind));
    %                 end
    % 
    %             %Evolution across theta sequences
    %                 data_sub = compute_vecBuffered(data_sub,50);
    %                 for j =  1:length(thetaSequenceNum_evolution)
    %                     data_relativeAngleToRat_evolution(i,j) = nanmean(abs(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),4)));
    %                     data_relativeAngleToPrev_evolution(i,j) = nanmean(abs(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),5)));
    % 
    %                     data_ratSpeed_evolution(i,j) = nanmean(abs(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),6)));
    % 
    %                     data_sub_sub = sign(data_sub(thetaSequenceNum_evolution(j,1):thetaSequenceNum_evolution(j,2),5));
    %                     data_sub_sub = abs(data_sub_sub(1:end-1)+data_sub_sub(2:end))==0;
    % 
    %                     data_pairwise_opposite_evolution(i,j) = nanmean(data_sub_sub);
    %                 end
    %                 data_relativeAngleToRat_evolution_mat(i,:) = compute_vecBuffered(data_sub(:,4),50)';
    %                 data_relativeAngleToPrev_evolution_mat(i,:) = compute_vecBuffered(data_sub(:,5),50)';
    % 
    %             %one-sided angle mean and variance
    %                 data_relativeAngleToRat_meanVariance_oneSided(i,1) = nanmean(abs(data_sub(:,4)));
    %                 data_relativeAngleToRat_meanVariance_oneSided(i,2) = nanvar(abs(data_sub(:,4)));
    %                 data_relativeAngleToPrev_meanVariance_oneSided(i,1) = nanmean(abs(data_sub(:,5)));
    %                 data_relativeAngleToPrev_meanVariance_oneSided(i,2) = nanvar(abs(data_sub(:,5)));
    % 
    %         end
    % 
    %         [~,theta_rat_running] = compute_dataTemporalConcatenation([theta_rat(:,1) theta_rat(:,1)],rat_runningStraight_timeBounds);
    %         ind_thetaCycles = find(~isnan(theta_rat_running(:,2)));
    % 
    %         % subplot(211)
    %         % delT = squareform(pdist(theta_rat(ind_thetaCycles,1))); delT = delT(:);
    %         % delX = squareform(pdist(sign(theta_seqDir(ind_thetaCycles,6)))); delX = delX(:)/2;
    %         % edges = linspace(0,2,40); centers = edges(1:end-1) + mean(diff(edges))/2;
    %         % data = nan(3,length(edges)-1);
    %         % for i = 1:length(edges)-1
    %         %     ind = find(delT>=edges(i) & delT<edges(i+1));
    %         %     data(1,i) = nanmean(delX(ind));
    %         %     data(2,i) = nanstd(delX(ind))./sqrt(length(ind));
    %         %     data(3,i) = length(ind);
    %         % end
    %         % shadedErrorBar(centers,100*data(1,:),100*data(2,:),'lineprops',{'k','linewidth',2})
    %         % xlabel('time between events (sec)'), ylabel('% alternation'), set(gca,'fontsize',14)
    %         % hline(50,'k--'), axis square, box on, axis tight
    %         % 
    %         % subplot(212)
    %         % edges = linspace(0,pi,16); centers = edges(1:end-1)+mean(diff(edges))/2;
    %         % h1 = histcounts(abs(theta_seqDir(ind_thetaCycles,2)),edges,'Normalization','pdf');
    %         % h2 = histcounts(abs(theta_seqDir(ind_thetaCycles,3)),edges,'Normalization','pdf');
    %         % yyaxis left, plot(180/pi*centers,h2,'k','linewidth',4)
    %         %     ylabel('fraction')
    %         % yyaxis right, plot(180/pi*centers,h1,'r','linewidth',4)
    %         %     ylabel('fraction')
    %         % xlabel('angle (deg)'), axis square, set(gca,'fontsize',14,'xtick',[0 90 180]), set(gcf,'color','w'), xlim(180*[0 1])
    %         % ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; 
    %         % h = gca; set(h.YLabel, 'Rotation', -90) % Rotate ylabel
    %         % legend('rel. to prev. sweep','rel. to rat heading')
    %         % ax = gca;
    %         % ax.YAxis(1).Color = 'k';
    %         % ax.YAxis(2).Color = 'r';
    %         % 
    %         % keyboard
    % 
    % 
    %         % subplot(211)
    %         % yyaxis left, plot(180/pi*centers_angles,nanmean(data_relativeAngleToPrev),'k','linewidth',4)
    %         %     ylabel({'fraction','(relative to prev sweep)'})
    %         % yyaxis right, plot(180/pi*centers_angles,nanmean(data_relativeAngleToRat),'k:','linewidth',4)
    %         %     ylabel({'fraction','(relative to rat heading)'})
    %         % xlabel('angle (deg)'), axis square, set(gca,'fontsize',14,'xtick',[0 90 180]), set(gcf,'color','w'), xlim(180*[0 1])
    %         % ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; 
    %         % h = gca; set(h.YLabel, 'Rotation', -90) % Rotate ylabel
    %         % legend('rel. to prev. sweep','rel. to rat heading')
    %         % 
    %         % subplot(212)
    %         % shadedErrorBar(centers_time_pairwise,100*nanmean(data_pairwise_opposite,1),100*nanstd(data_pairwise_opposite)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
    %         % axis square, xlabel('time between events (sec)'), ylabel('% alternation'), set(gca,'fontsize',14), axis tight
    %         % hline(50,'k--')
    %         % 
    %         % keyboard
    % 
    % 
    % 
    % 
    % 
    %         % subplot(121)
    %         % ind_thetaCycles = find(theta_rat(:,4)>0.4);
    %         % activityCOM_rot_full_sub = [reshape(activityCOM_rot_full(ind_thetaCycles,:,1),[],1) reshape(activityCOM_rot_full(ind_thetaCycles,:,2),[],1)];
    %         % edges_x = linspace(-5,5,50); centers_x = edges_x(1:end-1) + mean(diff(edges_x))/2;
    %         % edges_y = linspace(-3,7,50); centers_y = edges_y(1:end-1) + mean(diff(edges_y))/2;
    %         % h = histcounts2(activityCOM_rot_full_sub(:,1),activityCOM_rot_full_sub(:,2),edges_y,edges_x,'Normalization','pdf'); 
    %         % h = pcolor(binSize*centers_x,binSize*centers_y,h); set(h,'edgecolor','none')
    %         % caxis([0 0.012]), 
    %         % cb = colorbar(); ylabel(cb,'fraction','FontSize',14,'Rotation',270)
    %         % xlabel('x (cm)'), ylabel('y (cm)')
    %         % axis square, set(gca,'fontsize',14), set(gcf,'color','w')
    %         % vline(0,'k--'), hline(0,'k--')
    % 
    %         % subplot(223)
    %         % yyaxis left, shadedErrorBar(centers_thetaSequenceNum_evolution,180/pi*nanmean(data_relativeAngleToRat_evolution),180/pi*nanstd(data_relativeAngleToRat_evolution)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
    %         % ylabel({'angle relative to','rat heading (deg)'})
    %         % yyaxis right, plot(nanmean(timeBins_evolution,2),meanRatSpeed(1,:),'b','linewidth',2)
    %         % xlabel('theta sequence #'), ylabel('rat speed')
    %         % axis square, set(gca,'fontsize',14), box on, set(gcf,'color','w'),
    %         % ax = gca;
    %         % ax.YAxis(1).Color = 'k';
    %         % ax.YAxis(2).Color = 'b';
    %         % ylim([0.15 0.45])
    %         % 
    %         % subplot(224)
    %         % yyaxis left, shadedErrorBar(centers_thetaSequenceNum_evolution,100*nanmean(data_pairwise_opposite_evolution),100*nanstd(data_pairwise_opposite_evolution)./sqrt(numRunningPeriods),'lineprops',{'k','linewidth',2})
    %         % ylabel('% alternation'), 
    %         % yyaxis right, plot(nanmean(timeBins_evolution,2),meanRatSpeed(1,:),'b','linewidth',2)
    %         % xlabel('theta sequence #'), ylabel('rat speed')
    %         % axis square, set(gca,'fontsize',14), box on, set(gcf,'color','w'),
    %         % ax = gca;
    %         % ax.YAxis(1).Color = 'k';
    %         % ax.YAxis(2).Color = 'b';
    %         % ylim([0.15 0.45])
    % 
    % 
    %         data_relativeAngleToRat_cell{kk,kkk} = data_relativeAngleToRat;
    %         data_relativeAngleToPrev_cell{kk,kkk} = data_relativeAngleToPrev;
    %         data_pairwise_opposite_cell{kk,kkk} = data_pairwise_opposite;
    %         data_distance_cell{kk,kkk} = data_distance;
    %         data_relativeAngleToRat_evolution_cell{kk,kkk} = data_relativeAngleToRat_evolution;
    %         data_relativeAngleToPrev_evolution_cell{kk,kkk} = data_relativeAngleToPrev_evolution;
    %         data_pairwise_opposite_evolution_cell{kk,kkk} = data_pairwise_opposite_evolution;
    %         data_relativeAngleToRat_meanVariance_oneSided_cell{kk,kkk} = data_relativeAngleToRat_meanVariance_oneSided;
    %         data_relativeAngleToPrev_meanVariance_oneSided_cell{kk,kkk} = data_relativeAngleToPrev_meanVariance_oneSided;
    % 
    %    end
    % end
    % 
    % %%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%
    % 
    % for k = 1:3
    %     if k == 1, Data = data_relativeAngleToPrev_cell; xticks = 180/pi*(centers_angles-mean(diff(centers_angles)/2)); clabel = 'fraction';
    %     elseif k == 2, Data = data_pairwise_opposite_cell; xticks = centers_time_pairwise; clabel = '% alternation';
    %     elseif k == 3, Data = data_distance_cell; xticks = centers_distance; clabel = 'fraction';
    %     elseif k == 4, Data = data_relativeAngleToRat_evolution_cell; xticks = nanmean(thetaSequenceNum_evolution,2); clabel = 'rel. angle to rat (deg)';
    %     elseif k == 5, Data = data_relativeAngleToPrev_evolution_cell; xticks = nanmean(thetaSequenceNum_evolution,2); clabel = 'rel. angle to prev. (deg)';
    %     elseif k == 6, Data = data_pairwise_opposite_evolution_cell; xticks = nanmean(thetaSequenceNum_evolution,2); clabel = '% alternation';
    %     end
    % 
    %     for kk = 1:4
    %         data = [];
    % 
    %         for kkk = 1:size(data_relativeAngleToRat_cell,2)
    %             data_sub = Data{kk,kkk};
    %                 if k==2, data_sub = 100*data_sub; end
    %                 if k==4 || k==5, data_sub = 180/pi*data_sub; end
    %                 if k==6, data_sub = 100*data_sub; end
    %             data = [data; nanmean(data_sub,1)];
    %         end
    %         if kk == 1, yticks = v_max_vec;
    %         elseif kk == 2, yticks = T_theta_max_vec;
    %         elseif kk == 3, yticks = w_adapt_vec;
    %         elseif kk == 4, yticks = tau_adapt_vec;
    %         end 
    % 
    %         figure(k)
    %         subplot(2,2,kk)
    %         h = pcolor(xticks,yticks,data); set(h,'edgecolor','none')
    %         axis square, set(gcf,'Color','w'), set(gca,'fontsize',14)
    %         colorbar
    % 
    %         if kk == 1, ylabel('v_{max}'), hline(0.3,'r--')
    %         elseif kk == 2, ylabel('T_{\Theta}'), hline(0.12,'r--')
    %         elseif kk == 3, ylabel('w_{adapt}'), hline(40,'r--')
    %         elseif kk == 4, ylabel('\tau_{adapt}'), hline(0.8,'r--')
    %         end 
    % 
    %         if k == 1 
    %             xlabel({'relative angle (deg)'})
    %             % caxis([0 1]), vline(0,'k--')
    %         elseif k == 2
    %             xlabel({'time between events (sec)'})
    %             % caxis([30 80])
    %         elseif k == 3
    %             xlabel({'distance (cm)'})
    %         elseif k >=4
    %             % caxis([25 50])
    %             xlabel('theta sequence #')
    %         end    
    % 
    %         cb = colorbar(); ylabel(cb,clabel,'FontSize',14,'Rotation',270)
    %     end
    % 
    % end
    % 
    % %%%%%%%%%%%%%%%%%%%
    % 
    % colors = {'k','r','b','g'};
    % for k = 1:2
    %     subplot(1,2,k)
    %     for kk = 1:4
    %         if kk == 1, xticks = v_max_vec/max(v_max_vec);
    %             elseif kk == 2, xticks = T_theta_max_vec/max(T_theta_max_vec);
    %             elseif kk == 3, xticks = w_adapt_vec/max(w_adapt_vec);
    %             elseif kk == 4, xticks = tau_adapt_vec/max(tau_adapt_vec);
    %         end 
    % 
    %         Data = nan(2,10);
    %         for kkk = 1:10
    %             Data(1,kkk) = nanmean(data_relativeAngleToPrev_meanVariance_oneSided_cell{kk,kkk}(:,k));
    %             Data(2,kkk) = nanstd(data_relativeAngleToPrev_meanVariance_oneSided_cell{kk,kkk}(:,k))./sqrt(length(data_relativeAngleToRat_meanVariance_oneSided_cell{kk,kkk}(:,k)));
    %         end
    % 
    %         shadedErrorBar(xticks,180/pi*Data(1,:),180/pi*Data(2,:),'lineprops',{'color',colors{kk},'linewidth',2})
    %         xlabel('% change in paramter'), set(gca,'fontsize',14), set(gcf,'color','w'), axis square, box on
    %     end
    %     if k==1, ylabel('angle mean (deg)'), else, ylabel('angle variance (deg)'), end
    %     if k==1, legend('v_{max}','T_{\Theta}','w_{adapt}','\tau_{adapt}'), end
    % end




       
        % %plot example trajectory with theta sweeps
        %     % data = [theta_rat(:,1) theta_timeBounds theta_seqDir]; 
        %     data = compute_dataTemporalConcatenation([theta_rat(:,1) theta_timeBounds theta_seqDir],rat_running_timeBounds(4,:));
        %     rat_run =  compute_dataTemporalConcatenation(rat,data(:,2:3));
        %     [~,thetaSequences_run_left] = compute_dataTemporalConcatenation([t activityCOM],data(data(:,5)<0,2:3));
        %     [~,thetaSequences_run_right] = compute_dataTemporalConcatenation([t activityCOM],data(data(:,5)>=0,2:3));
        % 
        %     plot(unwrap(n*rat_run(:,2)-0.5),n*rat_run(:,3)-0.5,'k:','linewidth',2)
        %     hold on, plot(thetaSequences_run_left(:,2),thetaSequences_run_left(:,3),'r.','markersize',12), hold off
        %     hold on, plot(thetaSequences_run_right(:,2),thetaSequences_run_right(:,3),'b.','markersize',12), hold off
        %     % hold on, plot(n*theta_rat_start(:,2)-0.5,n*theta_rat_start(:,3)-0.5,'ko'), hold off
        %     axis([0 n 0 n]), set(gca,'xtick',[],'ytick',[]), axis square, set(gcf,'color','w')
        % 
        %     % plot(180/pi*theta_seqDir(:,2),'k','linewidth',2)
        %     % hline(0,'k--')
        %     % xlabel('theta sequence #'), ylabel({'angle relative to','rat heading (deg)'})
        %     % set(gca,'fontsize',14), set(gcf,'color','w'), axis square, axis tight

        % %evolution of angle relative to rat
        %     edges_angles = linspace(-pi,pi,40); centers_angles = edges_angles(1:end-1) + mean(diff(edges_angles))/2;
        %     timeBins_evolution = load_timeBins([1 30],1,3);
        %     H = [];
        %     for i = 1:size(timeBins_evolution,1)
        %         data_sub = data_relativeAngleToRat_evolution_mat(:,timeBins_evolution(i,1):timeBins_evolution(i,2));
        %         H = [H, histcounts(data_sub(:),edges_angles,'normalization','probability')'];
        %     end
        % 
        %     h = pcolor(nanmean(timeBins_evolution,2),180/pi*centers_angles,H); set(h,'edgecolor','none')
        %     xlabel('theta sequence #'), ylabel({'angle relative to','rat heading (deg)'})
        %     axis square, set(gca,'fontsize',14), box on, set(gcf,'color','w'), axis tight
        %     hline(0,'k--')
        %     caxis([0 0.1])
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    
