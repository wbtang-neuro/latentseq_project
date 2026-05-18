% function interactiveEllipse()
%     close all
% 
%     % Create a figure for the sliders and plot
%     fig = figure('Name', 'Interactive Ellipse', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);
%     set(gcf,'color','w')
% 
%     % Create axes for the plot
%     ax = axes('Parent', fig, 'Position', [0.4, 0.4, 0.4, 0.4]);
%     ax_2 = axes('Parent', fig, 'Position', [0.4, 0.1, 0.4, 0.2]);
% 
%     axis equal;
%     grid on;
%     hold on;
%     xlabel('x');
%     ylabel('y');
% 
%     % Default ellipse parameters
%     defaultL = 0.12/2;
%     defaultW = 0.06/2;
%     defaultR = 0.5;
%     defaultSpeed = 1;
%     defaultTheta0 = 39*(pi/180);
% 
% 
%     % Slider for semi-major axis (l)
%     uicontrol('Style', 'text', 'Position', [20, 420, 100, 20], 'String', 'Semi-Major (l)');
%     sliderL = uicontrol('Style', 'slider', 'Min', 0, 'Max', 0.5, 'Value', defaultL, ...
%                         'Position', [20, 400, 150, 20], 'Callback', @(src, ~) updatePlot());
% 
%     % Slider for semi-minor axis (w)
%     uicontrol('Style', 'text', 'Position', [20, 360, 100, 20], 'String', 'Semi-Minor (w)');
%     sliderW = uicontrol('Style', 'slider', 'Min', 0, 'Max', get(sliderL, 'Value'), 'Value', defaultW, ...
%                         'Position', [20, 340, 150, 20], 'Callback', @(src, ~) updatePlot());
% 
%     % Slider for rat starting location
%     uicontrol('Style', 'text', 'Position', [20, 300, 100, 20], 'String', 'displacement (r)');
%     sliderR = uicontrol('Style', 'slider', 'Min', 0, 'Max', 1, 'Value', defaultR, ...
%                             'Position', [20, 280, 150, 20], 'Callback', @(src, ~) updatePlot());
% 
%     % Slider for rat starting location
%     uicontrol('Style', 'text', 'Position', [20, 240, 100, 20], 'String', 'speed (v)');
%     sliderSpeed = uicontrol('Style', 'slider', 'Min', 0, 'Max', 2, 'Value', defaultSpeed, ...
%                             'Position', [20, 220, 150, 20], 'Callback', @(src, ~) updatePlot());
% 
%     % Slider for rat starting location
%     uicontrol('Style', 'text', 'Position', [20, 180, 100, 20], 'String', 'angle (theta0)');
%     sliderTheta0 = uicontrol('Style', 'slider', 'Min', 0, 'Max', pi, 'Value', defaultTheta0, 'SliderStep',[pi/1000 pi/100],...
%                             'Position', [20, 160, 150, 20], 'Callback', @(src, ~) updatePlot());
% 
% 
%     % Initial plot
%     updatePlot();
% 
%     % Nested function to update the plot
%     function updatePlot()
%         % Get slider values
%         l = get(sliderL, 'Value');
%         w = get(sliderW, 'Value');
%         r = l*get(sliderR, 'Value');
%         v = get(sliderSpeed, 'Value');
%         theta0 = get(sliderTheta0, 'Value');
%         delT = 1/8;
% 
%         % Generate ellipse points
%         t = linspace(0, 2*pi, 500);
%         x_ellipse = l * cos(t);
%         y_ellipse = w * sin(t);
% 
%         x_ellipse_grad = 2.8*l * cos(t);
%         y_ellipse_grad = 2.8*w * sin(t);
% 
%         x0 = -r; y0 = 0;
%         x = -r + v*delT*cos(theta0);
%         y = v*delT*sin(theta0);
%         theta = acos((x*cos(theta0)/l^2 + y*sin(theta0)/w^2)/sqrt(x^2/l^4 + y^2/w^4));
% 
%         theta0_vec = linspace(0,pi,40)';
%         x_vec = -r + v*delT*cos(theta0_vec);
%         y_vec = v*delT*sin(theta0_vec);
%         theta_vec = acos((x_vec.*cos(theta0_vec)/l^2 + y_vec.*sin(theta0_vec)/w^2)./sqrt(x_vec.^2/l^4 + y_vec.^2/w^4));
% 
%         % [X,Y] = meshgrid(linspace(-0.3,0.3,100),linspace(-0.3,0.3,100)); X = X(:); Y = Y(:);
%         % theta_mat = atan2(2*Y/l^2,2*X/l^2);
% 
%         % R = [cos(theta0+theta) -sin(theta0+theta); sin(theta0+theta) cos(theta0+theta)];
%         % X_ellipse_rot = R*[x_ellipse; y_ellipse];
%         % x_ellipse_rot = X_ellipse_rot(1,:) + r + v*delT*cos(theta0+theta);
%         % y_ellipse_rot = X_ellipse_rot(2,:) + v*delT*sin(theta0+theta);
% 
% 
%         % Clear and re-plot
%         cla(ax); % Clear axes
%         plot(ax, x_ellipse, y_ellipse,'k', 'LineWidth', 3);
%         hold(ax,'on'), plot(ax, x_ellipse_grad, y_ellipse_grad,'k:', 'LineWidth', 3); hold(ax,'off')
%         % hold(ax,'on'), plot(ax, x_ellipse_rot, y_ellipse_rot,'r', 'LineWidth', 2); hold(ax,'off')
%         hold(ax,'on'), plot(ax,x0,y0,'ko'), hold(ax,'off')
%         hold(ax,'on'), plot(ax,[x0,x],[y0,y],'g','linewidth',3), hold(ax,'off')
%         hold(ax,'on'), quiver(ax,x,y,cos(theta0+theta),sin(theta0+theta),0.1,'color','r','linewidth',2), hold(ax,'off')
%         xlim(ax, [-0.2 0.2]);
%         ylim(ax, [-0.1 0.2]);
%         title(ax, sprintf('l=%.2f, w=%.2f, r=%.2f, v=%.2f, theta0=%.2f, theta=%.2f', l, w, r, v, round(theta0*180/pi), round(theta*180/pi)));
%         axis(ax,'square'), box on
%         set(ax,'fontsize',14), axis(ax,'equal') 
%         hline(0,'k-'),vline(0,'k-')
% 
%         cla(ax_2); % Clear axes
%         plot(ax_2,180/pi*theta0_vec,180/pi*theta_vec,'k','linewidth',2),
%         yl = ylim;
%         hold(ax_2,'on'), plot(180/pi*theta0_vec,180/pi*theta0_vec,'k:','linewidth',3), hold(ax,'off')
%         hold(ax_2,'on'), plot(180/pi*[theta0,theta0],yl,'g-','linewidth',3), hold(ax_2,'off')
%         xlim(ax_2,180/pi*[min(theta0_vec) max(theta0_vec)])
%         ylim(ax_2,yl)
% 
%         box on, grid off
%         xlabel('\theta_1'), ylabel('\theta_2')
%         xlim([0 90])
%         set(ax_2,'xtick',[0:30:90],'fontsize',16)
%         set(ax_2,'ytick',[0:30:90],'fontsize',16)
% 
%         drawnow
%     end
% end

%%


    w = 0.06/2;
    r = 0.5;
    delT = 1/8;

    theta0_vec = linspace(0,pi,100)';
    l_vec = w*(1:0.25:3);
    v_vec = linspace(0.01,1,100)';
    theta_star = nan(10,100,2);

    for j = 1:length(l_vec)
        for i = 1:length(v_vec)
            l = l_vec(j);
            a = l/w;

            x_vec = -r*l + v_vec(i)*delT*cos(theta0_vec);
            y_vec = v_vec(i)*delT*sin(theta0_vec);

            % theta_vec = acos((x_vec.*cos(theta0_vec)/l^2 + y_vec.*sin(theta0_vec)/w^2)./sqrt(x_vec.^2/l^4 + y_vec.^2/w^4));
            theta_vec = acos((x_vec.*cos(theta0_vec)/a^2 + y_vec.*sin(theta0_vec))./sqrt(x_vec.^2/a^4 + y_vec.^2));

            [a,b] = InterX([theta0_vec theta0_vec],[theta0_vec theta_vec]);
            ind = b(end,1);
            if ind==1; continue, end

            theta_star(j,i,1) = real(a(end,1));
            if isnan(theta_star(j,i,1)), continue, end
            theta_star(j,i,2) = abs((theta_vec(ind+1) - theta_vec(ind))/(theta0_vec(ind+1) - theta0_vec(ind)));

            % plot(theta0_vec,theta_vec)
            % hold on, plot(theta0_vec,theta0_vec,'r:'), hold off
            % keyboard
        end
    end


    a = cool(length(l_vec));
    for j = 1:length(l_vec)
        l = l_vec(j);
        % r = sqrt(l^2-w^2)/l;

        % Generate ellipse points
        t = linspace(0, 2*pi, 500);
        x_ellipse = l * cos(t);
        y_ellipse = w * sin(t);

        % subplot(411)
        subplot(131)
        if j==1, plot(0,0,'ko','markersize',12,'linewidth',2), hold on, plot(x_ellipse+r*l,y_ellipse,'color',a(j,:),'linewidth',3), hold off
        else, hold on, plot(0,0,'ko','markersize',12,'linewidth',2), plot(x_ellipse+r*l,y_ellipse,'color',a(j,:),'linewidth',3), hold off
        end
        axis tight, axis equal, axis off

        % subplot(4,1,[2 4])
        subplot(132)
        if j==1, plot(v_vec,180/pi*theta_star(j,:,1),'color',a(j,:),'linewidth',3)
        else, hold on, plot(v_vec,180/pi*theta_star(j,:,1),'color',a(j,:),'linewidth',3), hold off
        end
        xlabel('rat speed (m/s)'), ylabel('|sweep angle| (deg)')
        axis square, 
        ylim([0 90]), set(gca,'ytick',[0:30:90])

        subplot(133)
        if j==1, plot(v_vec,theta_star(j,:,2),'color',a(j,:),'linewidth',3)
        else, hold on, plot(v_vec,theta_star(j,:,2),'color',a(j,:),'linewidth',3), hold off
        end
        xlabel('rat speed (m/s)'), ylabel('|slope|')
        axis square,
        hline(1,'k:')

    end

    set(gcf,'color','w')
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 18); 

    %%