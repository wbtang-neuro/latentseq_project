use_spiking = 1;
use_feedforward = 1;
use_envelope = 0;
use_periodic = 1;

spatialDim = 2;
use_plot = 1;

%parameters
    %simulation duration
    T_run = 8;
    numStops = 1; %number of stops
    T = numStops*(T_run); %total simulation time

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
            w_adapt = 10; %amplitude of adaptation
            tau_adapt = 0.8; %adaptation time constant !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            w_FF = 3000; %amplitude of feedforward exc.
            sig_FF = 0.1; %width of feedforward exc.
            v_max = 0.1;
            T_theta_max = 0.4; 
            T_replay_min = 0.8;
        else
            tau_s = 20/1000; %synaptic time constant
            beta_0 = 600; %uniform excitation    
            w_rec = 24; %strength of recurrent exc.
            sig_rec = 0.14;%0.02; %width of recurrent exc.
            w_inh = w_rec; %global inhibition
            w_adapt = 40; %amplitude of adaptation 
            tau_adapt = 0.8; %adaptation time constant 
            w_FF = 2200; %amplitude of feedforward exc.
            sig_FF = 0.08; %0.001; %width of feedforward exc.
            v_max = 0.3;
            T_theta_max = 0.12;
        end

        T_plot = T_theta_max;%/12;

        %place field centers
        if spatialDim == 2
            [X,Y] = meshgrid((1:n)/n,(1:n)/n); 
        else
            X = (1:n)'/n;
        end

        %activity envelope
        if spatialDim == 1
            A = ones(n,1);
        else
            A = ones(n,n);
        end

%rat trajectory
    if spatialDim == 1
        %rat max speed
        % v_max = 0.1;

        v_lap = v_max*[ones(T_run/dt,1)];
        v = repmat(v_lap,1,numStops); v = v(:);
            if use_plot == 1, v = circshift(v,-(T_run/4)/dt); end %uncomment this if sim starts with rat running full speed
        x = mod(0.1+cumsum(v)*dt,1);
        acc = [abs(diff(v)); 0];

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

            %fixed max speed
            v_lap = v_max*[ones(T_run/dt_randomWalk,1)];
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
                g_FF = normpdf([X-1; X; X+1],x(t,:),sig_FF); g_FF = w_FF*g_FF/max(g_FF(:));
                g_FF = g_FF(1:n) + g_FF(n+1:2*n) + g_FF(2*n+1:3*n);
            end
            if use_feedforward==0, g_FF = 0; end
            g_feedbackInh = -w_inh*sum(r(:));
            g_unifInput = beta_0;
            g_adapt = -w_adapt*a;

            G = A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_FF;
            % G = g_FF;

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
        if use_plot==1 &&  mod(t,T_plot/dt)==0% && t>1/dt    
            if spatialDim == 1
                % subplot(331), plot(r/max(r)), hold on, plot(g_FF/max(g_FF)),hold off
                % subplot(331), plot(r,'k','linewidth',2), ylim([0 3]), vline(n*x(t),'r'), xlabel('neuron'), ylabel('firing rate')
                % subplot(334), plot(a,'k','linewidth',2), vline(n*x(t),'r'), xlabel('neuron'), ylabel('adaptation'), %ylim([0 30])
                %     if w_facil ~=0, hold on, plot(f,'b','linewidth',2), hold off, end
                % subplot(337), plot(dt:dt:T,v,'k','linewidth',2), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
                % subplot(1,3,[2 3])

                % subplot(121)
                activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread) | [compute_sequenceJumps(activityCOM_sub); 0]>1) = nan;
                % plot(dt*(1:t),activityCOM_sub(1:t),'r.','markersize',8)
                cplot(dt*(1:t),activityCOM_sub(1:t),1:t,'.','markersize',12)
                hold on, plot(t*dt,n*x(t)-0.5,'ko','linewidth',2), hold off
                hold on, plot((1:t)*dt,n*x(1:t)-0.5,'k.','markersize',4), hold off
                ylim([0 n]), xlabel('time (sec)'), ylabel('neuron')
                axis square
                set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
            else
                % subplot(331), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('firing rate'), clim([0 3])
                % subplot(334), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('adaptation (-)'), clim([0 30])
                % subplot(337), plot(dt:dt:T,v), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
                % subplot(231), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('Firing rate'), clim([0 3])
                % subplot(234), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('Adaptation'), clim([0 30])
                % subplot(1,3,[2 3]), 

                % subplot(122)
                activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread(:,1)) | [compute_sequenceJumps(activityCOM_sub); 0]>1,:) = nan;
                % plot((activityCOM_sub(1:t,1)),(activityCOM_sub(1:t,2)),'r.','markersize',8)
                % hold on, cplot((activityCOM_sub(t-1/dt+1:t,1)),(activityCOM_sub(t-1/dt+1:t,2)),1:1/dt,'.','markersize',12), hold off
                cplot((activityCOM_sub(1:t,1)),(activityCOM_sub(1:t,2)),1:t,'.','markersize',12)
                hold on, plot((n*x(t,1)-0.5),(n*x(t,2)-0.5),'ko','linewidth',2), hold off
                hold on, plot((n*x(1:t,1)-0.5),(n*x(1:t,2)-0.5),'k','linewidth',2), hold off
                axis([0 n 0 n]), %set(gca,'xtick',[],'ytick',[]), 
                axis square, %title('Center-of-mass trajectory')
                xlabel('neuron'), ylabel('neuron')
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