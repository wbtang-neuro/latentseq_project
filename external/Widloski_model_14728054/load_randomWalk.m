function [x] = load_randomWalk(boundaries,T,s_mean,s_var,theta_var,dt,barriers,spatialDim)%,thr_wallLingering)

%sample duration of random walk
%     T = poissrnd(T_mean);
%     T = round(exprnd(T_mean));

%initial conditions
    % x_prev = round(diff(boundaries)*rand(1,spatialDim)+min(boundaries));
    x_prev = 100*[0.5 0.1]; %!!!!!!!!!!!!!!
    s_prev = s_mean;

    if spatialDim==1 
        theta_prev = binornd(1,0.5)*2-1;
    else
        theta_prev = pi/2; % !!!!!!!!!!!!!!!!!!!!!
        % theta_prev = rand*2*pi;
    end

%empty storage matrix for trajectory
    x = nan(round(T/dt),spatialDim); 
    x(1,:) = x_prev;

if size(s_mean,1)==1, s_mean = s_mean*ones(T/dt,1); end
    
%random walk dynamics
    for j = 2:T/dt
        wall_crossings = 1;
        barrier_crossings = ones(length(barriers),1);
        iter = 1;
        while wall_crossings+sum(barrier_crossings)>0
            
            s_new = abs(normrnd(s_mean(j),s_var));
            if spatialDim==1
                % theta_new = normrnd(theta_prev,iter*theta_var);
                % delx = s_new*cos(theta_new)*dt;

                b = binornd(1,iter*theta_var);
                if b==0 
                    theta_new = theta_prev; 
                else
                    theta_new = -theta_prev;
                end
                delx = s_new*theta_new*dt;
            elseif spatialDim == 2
                theta_new = normrnd(theta_prev,iter*theta_var);
                delx = s_new*[cos(theta_new),sin(theta_new)]*dt;
            end
            x_new = x_prev + delx;
            
            %check for wall crossings
            if spatialDim==1
                if x_new(1)>boundaries(2) || x_new(1)<boundaries(1)
                    wall_crossings = 1;
                else
                    wall_crossings = 0;
                end
            elseif spatialDim==2
                if x_new(1)>boundaries(2) || x_new(1)<boundaries(1) || x_new(2)>boundaries(2) || x_new(2)<boundaries(1)
                    wall_crossings = 1;
                else
                    wall_crossings = 0;
                end
            end
            
            %check for barrier crossings
            if spatialDim==2
                for k = 1:length(barriers)
                    barrier_crossings(k) = ~isempty(InterX([x_prev;x_new]',barriers{k}'));
                end
            end
                
            iter = iter+1; 
                        
        end
        
        x(j,:) = x_new;
        x_prev = x_new;
        s_prev = s_new;
        theta_prev = theta_new;

        % if spatialDim==1
        %     plot(x(1:j))
        % else
        %     plot(x(1:j,1),x(1:j,2),'.-')
        %     axis([0-0.1 1+0.1 0-0.1 1+0.1])
        % end
        % drawnow
    end
    
    x = [(dt:dt:T)' x];