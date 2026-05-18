function [outputArg1,outputArg2] = scatterEgodirAcorrs(p)
%SCATTEREGODIRACORRS Summary of this function goes here
%   Detailed explanation goes here
arguments
p.cols = [.2,.2,.2; 1,0,0];
p.acorrs = [];
p.ax = [];
p.tl = [];
p.plotHalf = 0;
p.shadedThetaCycles = 1;
end
acorrs = p.acorrs;
nlags = size(acorrs, 2);
lags = -floor(nlags/2):floor(nlags/2);
if p.plotHalf
    acorrs = acorrs(:, ceil(nlags/2):end);
    lags = lags(:, ceil(nlags/2):end);
    nlags = size(acorrs, 2);
end

cols = p.cols;
if isempty(p.ax)
    p.ax = nexttile;
end
ax = p.ax;
ax.YLim =[-.65,.65];
ax.YLim =[-1,1];
ax.XLim = minmax(lags)+[-.5,.5];
if p.shadedThetaCycles
    x = ([0,1,1,0,0]-.5);
    y = [vec(ax.YLim+[0;0]+.01);ax.YLim(1)]';
    cens = lags(2:2:end);
    for c = cens
       patch(ax, x+c, y, [1,1,1]*.9, 'faceAlpha', 1, 'edgeColor', 'none');
    end
end

clear means
for g = 1:size(acorrs, 2)
    y = [acorrs(:, g)];
    x = ones(size(y))*lags(g);
    swarmchart(x, y, 8, 'MarkerFaceColor', [.5,.5,.5], 'XJitterWidth',.7)
    plot([1,1]*lags(g), [-std(y), std(y)]+mean(y),'k', 'LineWidth',1);

    errorbar([1,1]*lags(g), mean(y), std(y),'k', 'LineWidth',1);
% errorbar(0, mean(coflickering), std(coflickering)./numel(coflickering), 'k')
    h(g) = scatter(lags(g), mean(y), 30, 'MarkerEdgeColor','k', 'MarkerFaceColor','r', 'LineWidth',1);
    means(g) = mean(y);
    %     text(x(1)+2, .9, sprintf("mean %.2f \n std %.2f", mean(y), std(y)))
end
plot(lags, means, 'k');
ylabel("Correlation (r)")
xlabel("Lags (theta cycles)")
yticks([-1,0,1])
yline(0)
ax = gca;
ax.FontSize = 10;
end

