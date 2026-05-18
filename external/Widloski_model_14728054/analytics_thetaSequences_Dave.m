function interactiveEllipse()
    % Create a figure for the sliders and plot
    fig = figure('Name', 'Interactive Ellipse', 'NumberTitle', 'off', 'Position', [100, 100, 600, 600]);

    % Create axes for the plot
    ax = axes('Parent', fig, 'Position', [0.3, 0.35, 0.6, 0.6]);
    ax_2 = axes('Parent', fig, 'Position', [0.3, 0.1, 0.6, 0.2]);

    axis equal;
    grid on;
    hold on;
    xlabel('x');
    ylabel('y');

    % Default ellipse parameters
    defaultL = 0.3/2;
    defaultW = 0.15/2;
    defaultR = 0.1;
    defaultSpeed = 2;
    defaultTheta0 = pi/4;

    % Slider for semi-major axis (l)
    uicontrol('Style', 'text', 'Position', [20, 420, 100, 20], 'String', 'Semi-Major (l)');
    sliderL = uicontrol('Style', 'slider', 'Min', 0, 'Max', 0.5, 'Value', defaultL, ...
                        'Position', [20, 400, 150, 20], 'Callback', @(src, ~) updatePlot());

    % Slider for semi-minor axis (w)
    uicontrol('Style', 'text', 'Position', [20, 360, 100, 20], 'String', 'Semi-Minor (w)');
    sliderW = uicontrol('Style', 'slider', 'Min', 0, 'Max', get(sliderL, 'Value'), 'Value', defaultW, ...
                        'Position', [20, 340, 150, 20], 'Callback', @(src, ~) updatePlot());

    % Slider for rat starting location
    uicontrol('Style', 'text', 'Position', [20, 300, 100, 20], 'String', 'displacement (r)');
    sliderR = uicontrol('Style', 'slider', 'Min', 0, 'Max', 1, 'Value', defaultR, ...
                            'Position', [20, 280, 150, 20], 'Callback', @(src, ~) updatePlot());

    % Slider for rat starting location
    uicontrol('Style', 'text', 'Position', [20, 240, 100, 20], 'String', 'speed (v)');
    sliderSpeed = uicontrol('Style', 'slider', 'Min', 0, 'Max', 10, 'Value', defaultSpeed, ...
                            'Position', [20, 220, 150, 20], 'Callback', @(src, ~) updatePlot());

    % Slider for rat starting location
    uicontrol('Style', 'text', 'Position', [20, 180, 100, 20], 'String', 'angle (theta0)');
    sliderTheta0 = uicontrol('Style', 'slider', 'Min', 0, 'Max', pi, 'Value', defaultTheta0, ...
                            'Position', [20, 160, 150, 20], 'Callback', @(src, ~) updatePlot());


    % Initial plot
    updatePlot();

    % Nested function to update the plot
    function updatePlot()
        % Get slider values
        l = get(sliderL, 'Value');
        w = get(sliderW, 'Value');
        r = l*get(sliderR, 'Value');
        v = get(sliderSpeed, 'Value');
        theta0 = get(sliderTheta0, 'Value');
        delT = 1/8;

        % Generate ellipse points
        t = linspace(0, 2*pi, 500);
        x_ellipse = l * cos(t);
        y_ellipse = w * sin(t);

        x0 = -r; y0 = 0;
        x = -r + v*delT*cos(theta0);
        y = v*delT*sin(theta0);
        theta = acos((x*cos(theta0)/l^2 + y*sin(theta0)/w^2)/sqrt(x^2/l^4 + y^2/w^4));
        
        % Del = [x^2/l^2,y^2/w^2];
        % theta_2 = atan2(Del(:,2),Del(:,1))-theta0;

        theta0_vec = linspace(0,pi,40)';
        x_vec = -r + v*delT*cos(theta0_vec);
        y_vec = v*delT*sin(theta0_vec);
        theta_vec = acos((x_vec.*cos(theta0_vec)/l^2 + y_vec.*sin(theta0_vec)/w^2)./sqrt(x_vec.^2/l^4 + y_vec.^2/w^4));

        % Del = [x_vec.^2/l^2,y_vec.^2/w^2];
        % theta_vec_2 = atan2(Del(:,2),Del(:,1))-theta0;

        % Clear and re-plot
        cla(ax); % Clear axes
        plot(ax, x_ellipse, y_ellipse,'k', 'LineWidth', 2);
        hold(ax,'on'), plot(ax,x0,y0,'ko'), hold(ax,'off')
        hold(ax,'on'), plot(ax,[x0,x],[y0,y],'g','linewidth',3), hold(ax,'off')
        hold(ax,'on'), quiver(ax,x,y,cos(theta0+theta),sin(theta0+theta),0.15,'color','r','linewidth',2), hold(ax,'off')
        % hold(ax,'on'), quiver(ax,x,y,cos(theta0+theta_2),sin(theta0+theta_2),0.15,'color','b','linewidth',2), hold(ax,'off')
        xlim(ax, [-0.25 0.5]);
        ylim(ax, [-0.25 0.5]);
        title(ax, sprintf('l=%.2f, w=%.2f, r=%.2f, v=%.2f, theta0=%.2f, theta=%.2f', l, w, r, v, theta0*180/pi, theta*180/pi));
        axis(ax,'square')

        cla(ax_2); % Clear axes
        plot(ax_2,theta0_vec,theta_vec,'k','linewidth',2),
        % hold(ax_2,'on'), plot(ax_2,theta0_vec,theta_vec_2,'b','linewidth',2), hold(ax,'off')
        yl = ylim;
        hold(ax_2,'on'), plot(theta0_vec,theta0_vec,'k:','linewidth',3), hold(ax,'off')
        hold(ax_2,'on'), plot([theta0,theta0],yl,'g-','linewidth',3), hold(ax_2,'off')
        xlim(ax_2,[min(theta0_vec) max(theta0_vec)])
        ylim(ax_2,yl)
        pbaspect(ax_2,[3 1 1])

        drawnow
    end
end


    l = 0.2/2;
    w = 0.06/2;
    r = 0.3;
    delT = 1/8;

    theta0_vec = linspace(0,pi,100)';
    l_vec = linspace(w,3*w,10);
    r_vec = linspace(0.01,1,10);
    v_vec = linspace(0.01,4,100)';
    theta_star = nan(10,100,2);

    for j = 1:10
        for i = 1:length(v_vec)
            % r = r_vec(j);
            l = l_vec(j);
            x_vec = -r + v_vec(i)*delT*cos(theta0_vec);
            y_vec = v_vec(i)*delT*sin(theta0_vec);

            theta_vec = acos((x_vec.*cos(theta0_vec)/l^2 + y_vec.*sin(theta0_vec)/w^2)./sqrt(x_vec.^2/l^4 + y_vec.^2/w^4));

            [a,b] = InterX([theta0_vec theta0_vec],[theta0_vec theta_vec]);
            ind = b(end,1);
            if ind==1; continue, end

            theta_star(j,i,1) = real(a(end,1));
            theta_star(j,i,2) = abs((theta_vec(ind+1) - theta_vec(ind))/(theta0_vec(ind+1) - theta0_vec(ind)));

            % plot(theta0_vec,theta_vec)
            % hold on, plot(theta0_vec,theta0_vec,'r:'), hold off
            % keyboard
        end
    end

    a = cool(10);
    for j = 1:10
        % r = r_vec(j);
        l = l_vec(j);

        % Generate ellipse points
        t = linspace(0, 2*pi, 500);
        x_ellipse = l * cos(t);
        y_ellipse = w * sin(t);

        subplot(131)
        % if j==1, plot(x_ellipse,y_ellipse,'k','linewidth',2), hold on, plot(-r*l,0,'o','color',a(j,:),'markersize',12,'linewidth',2), hold off 
        % else, hold on, plot(-r*l,0,'o','color',a(j,:),'markersize',12,'linewidth',2), hold off
        % end
        if j==1, plot(-r*l,0,'o','color',a(j,:),'markersize',12,'linewidth',2), hold on, plot(x_ellipse,y_ellipse,'color',a(j,:),'linewidth',2), hold off
        else, hold on, plot(-r*l,0,'o','color',a(j,:),'markersize',12,'linewidth',2), plot(x_ellipse,y_ellipse,'color',a(j,:),'linewidth',2), hold off
        end
        axis tight, axis equal, axis off

        subplot(132)
        if j==1, plot(v_vec,180/pi*theta_star(j,:,1),'color',a(j,:),'linewidth',2)
        else, hold on, plot(v_vec,180/pi*theta_star(j,:,1),'color',a(j,:),'linewidth',2), hold off
        end
        xlabel('speed'), ylabel('angle (deg)')
        axis square, set(gca,'fontsize',14)

        subplot(133)
        if j==1, plot(v_vec,theta_star(j,:,2),'color',a(j,:),'linewidth',2)
        else, hold on, plot(v_vec,theta_star(j,:,2),'color',a(j,:),'linewidth',2), hold off
        end
        xlabel('speed'), ylabel('|slope|')
        axis square, set(gca,'fontsize',14)
        hline(1,'k:')
    end
    set(gcf,'color','w')