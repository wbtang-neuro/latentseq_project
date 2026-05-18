function [fcn, nIn, nOut, prmNames, fsName, plotFcn] = basisFcnPreset(presetName)

switch lower(presetName)
    case 'vm'
        denom = @(kp) 2 * pi * besseli(0, kp);
        fcn = @(X, th, kp) exp(kp*cos(X-th)) ./ denom(kp);
        nIn = 1;
        nOut = 1;
        prmNames = {'kappa', 'thetaHat'};
        fsName = 'Von-Mises';
        plotFcn = @plot1dCirc;
    case 'rcos'
        fcn = @(X, mu, beta) rcos(X, mu, beta);
        nIn = 1;
        nOut = 1;
        prmNames = {'mu', 'beta'};
        fsName = 'raised-cosine';
        plotFcn = @plot1d;
    case 'rcoslog'
        fcn = @(X, mu, beta, s) rcoslog(X, mu, beta, s);
        nIn = 1;
        nOut = 1;
        prmNames = {'mu', 'beta', 's'};
        fsName = 'log-raised-cosine';
        plotFcn = @plot1d;
    case 'gaus'
        fcn = @(X, mu, sigma) normpdf(X, mu, sigma);
        nIn = 1;
        nOut = 1;
        prmNames = {'mu', 'sigma'};
        fsName = '1D-Gaussian';
        plotFcn = @plot1d;
    case 'gaus2'
        fcn = @(X, mu, sigma) mvnpdf(X, mu, sigma.^2); % mvnpdf takes *covariance*
        nIn = 2;
        nOut = 1;
        prmNames = {'mu', 'sigma'};
        fsName = '2D-Gaussian';
        plotFcn = @plot2d;
    case 'conv'
        fcn = @(X, v) conv(X, v, 'same');
        nIn = 1;
        nOut = 1;
        prmNames = {'kernel'};
        fsName = '1D-convolution';
        plotFcn = @plot1dConv;
    case 'conv-causal'
        fcn = @(X, v, nPad) convCausal(X, v, nPad);
        nIn = 1;
        nOut = 1;
        prmNames = {'kernel', 'npad'};
        fsName = '1D-causal-convolution';
        plotFcn = @plot1dConv;
end

end

function h = plot1d(ax, x, y)
    h = line(ax, x, y);
    ax.XLim = x([1, end]);
    xlabel(ax, 'X');
    ylabel(ax, 'Y');
end

function h = plot2d(ax, x, y, z)
    h = imagesc(ax, x, y, z);
    ax.YDir = "normal";
    ax.DataAspectRatio = [1, 1, 1];
%     axis(ax, 'xy', 'equal');
%     xlabel(ax, 'X_1');
%     ylabel(ax, 'X_2');
end

function h = plot1dCirc(ax, x, y)
    h = plot1d(ax, x, y);
%     ax.XTick = [0, pi, 2*pi];
    ax.XTick = [-pi, 0, pi];
    ax.XTickLabel = {'-\pi', '0', '\pi'};
    ax.XLim = [-pi, pi];
end

function h = plot1dConv(ax, x, y)
    h = line(ax, x, y);
    ax.XLim = x([1, end]);
    ylabel(ax, 'Y');
    xlabel(ax, 'index');
end

function y = rcoslog(x, mu, beta, s)
% log-transformed raised cosine
nlin = @(x)log(x+1e-20);
y = rcos(nlin(x+s), mu, beta);
end

function y = rcos(x, mu, beta)
% mu - center position of each cosine
% beta - x-scaling factor
b = (x-mu)*pi/beta/2;
b = min(pi, b);
b = max(-pi, b);
y = (cos(b)+1)/2;
end

function y = convCausal(x, h, nPad)
h = h(:);
nh = numel(h);
y = conv(x, h, 'full');
y = [zeros(nPad, 1); y(1:end-nh-nPad+1,:)];
end