% Examples of using the function-space tools

% FUNCTION SPACES
% A function space is a space defined by a set of functions, each of which
% acts as a basis vector. This means that a point in the function space
% represents a linear combination of the function values.
%
% A FunctionSpace object represents a collection of BasisFunction objects,
% which collectively define the function space.

%% Example 1a: create a high-D space with 1-D gaussian basis vectors.

% Create parameters for the gaussian basis functions.
nDims = 50;                     % Number of dimensions in function space
mu = linspace(0, 1, nDims)';    % Basis function mean (evenly spread)
sigma = 0.05 * ones(nDims, 1);  % Basis function sigma (same for all)

% Create the FunctionSpace and display its basis functions.
%
% The "createFcnSpace" function is a convenient way of creating a
% FunctionSpace, with some useful presets for particular types of basis
% function. Here we specify the function type "gaus", which indicates a 1-D
% gaussian function. We pass the gaussian mu and sigma parameters, and also
% a 1-D grid of coordinates which is used for generating plots.

plottingGrid = {linspace(0, 1, 200)'};
prms = num2cell([mu, sigma], 1);
F = createFcnSpace('gaus', prms, 'plottingCoords', plottingGrid);

% Display the basis function values across the default plotting grid
figure();
F.plotAll([], 'showLabels', true);
suptitle('FunctionSpace 1-D gaussian basis vectors');

%% Example 1b: Use the FunctionSpace to decompose a 1-D random variable.

% FunctionSpace objects can be used to decompose data in the basis function
% domain. The decomposition maps the data from the domain to the codomain, 
% for each of the basis functions.

% Create a random 1D variable
nPoints = 500;
sig = smooth(rand(nPoints, 1), 5);
sig = sig-min(sig);
sig = sig/max(sig);

% We can perform functional decomposition on the signal using the
% "evaluate" method, which passes signal "sig" to each of the basis
% functions and returns matrix Y, each column of which contains the values
% for one basis function.
YD = F.evaluate(sig);

% Plot original signal
figure();
subplot(2, 1, 1)
line(1:numel(sig), sig, 'color', 'k');
title('1D signal');
ylim([0, 1]);

% Plot decomposed signal
subplot(2, 1, 2);
imagesc(YD')
xlabel('Time');
ylabel('Basis function');
axis xy
title(sprintf('%u-D decomposition', nDims));

%% Example 1b: combine basis functions to reconstruct a tuning curve

% The complementary process to decomposing a 1D variable using basis
% functions is combining basis functions in particular proportions to
% reconstruct a signal.

% Let's simulate a variable that is tuned to the random signal we generated
% above. This function describes the variable's tuning curve.
tuningFcn = @(x) exp(5*-abs(x-0.5)) .* sin(2 * 2*pi*x);
yReg = tuningFcn(sig) + 0.1 * randn(nPoints, 1);

% Use regression to fit a tuning curve of yReg to the decomposed signal.
b = regress(yReg, [ones(nPoints, 1), YD]);
c = b(1);  % intercept
b(1) = [];

% Reconstructing the fitted tuning curve:
%
% Each element of vector "b" is a regression coefficient for one basis
% function, indicating the weighting of the basis function determined by
% the regression model. The weighted combination of basis functions 
% represented by vector "b" defines a point in the function space.
% The fitted tuning curve can be resconstructed by evaluating a grid of
% coordinates at point "b".
%
% We can reconstruct the tuning curve by mixing together the basis
% functions with weights given by the regression coefficients. We use the
% "evaluate" method again, using the grid of points over which we wish to
% construct the tuning curve, and giving "b" as the second argument to
% indicate the point in the function space (the weighting of the basis
% functions).
%
% Typically, this reconstruction shows signs of overfitting, due to the
% large number of basis functions compared to the number of observations.
% The dimensionality reduction in the next section provides a remedy
% against this.

tmp = F.evaluate(plottingGrid{1}, b');
tCurveFit = c + sum(tmp, 2);

figure();
subplot(2, 3, 1);
plot(sig, 'k');
title('Predictor signal');

subplot(2, 3, 2);
imagesc(YD');
title('Decomposed predictor');

subplot(2, 3, 3);
plot(plottingGrid{1}, tuningFcn(plottingGrid{1}), 'k');
title('Response tuning curve')
xlabel('Predictor value');
ylabel('Response value');

subplot(2, 3, 4);
plot(yReg, 'k');
xlabel('Time');
ylabel('Response');
title('Simulated response');

subplot(2, 3, 5);
bar(b);
title('Regression coefficients');
ylabel('b');

subplot(2, 3, 6);
plot(plottingGrid{1}, tCurveFit, 'k');
xlabel('Signal value');
ylabel('Fitted response');
title('Fitted tuning curve');
set(gca, 'ylim', [-1, 1]);

%% Example 2a: use a DRFunctionSpace to eliminate redundant dimensions

% The DRFunctionSpace is an extension of FunctionSpace which uses PCA to
% transform the basis function output, yielding a new orthogonal basis from
% which redundant dimensions can be discarded.
%
%             decomposition           PCA
%          X --------------> Y  --------------> Y'
%       DOMAIN            CODOMAIN           CODOMAIN
%
% As the above diagram shows, a DRFunctionSpace maps variable X in the 
% function domain, to variable Y' in the function codomain. Y' is a
% linearly transformed version of Y, yielded by multiplying Y by its
% eigenvectors.
%
% Despite the additional PCA step, a DRFunctionSpace can be used in the 
% same way as normal FunctionSpace for decomposing a variable or
% reconstructing a function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% We can create a DRFunctionSpace from an existing FunctionSpace object
% using the FunctionSpace "dimReduce" method. For this, we need to define
% the coordinates that we want to use for calculating the PCA coefficients.
% We can just use a uniform grid for this, which will ensure that all basis
% functions are equally prioritized by the PCA.
%
% The "dimReduce" method runs the PCA and returns a DRFunctionSpace object
% "FR".
xDecomp = linspace(0, 1, 1000)';
FR = F.dimReduce(xDecomp);

figure();
FR.plotAll();
suptitle('DRFunctionSpace gaussian-PCA basis vectors');

% FR will still contain all dimensions at this point. We can eliminate
% uninformative dimensions using the "discardThresh" method, which applies
% a threshold for how much variance must be explained, beyond which all
% further dimensions are discarded.
%
% Here we discard all dimensions above the 99th percentile
pcThresh = 0.99;
explained = FR.explained;
FR.discardThresh(0.99);
nKeep = FR.nDims;

figure()
ax = gca();
bar(explained * 100)
xlabel('PC');
ylabel('%');
title('Variance explained');
ax.YLim = ax.YLim;
h = line([1, 1]*nKeep+0.5, ax.YLim, 'color', 'r', 'lineStyle', '--');
legend(h, {sprintf('%.0f%% threshold', pcThresh*100)});


%% Example 2c: rescaling dimensions of a function space

% Before PCA, our gaussian decomposition had similar similar variances
% across all basis functions. However PCA drastically changes the 
% distribution of variance, and having large differences in variance
% can make fitting a regression more difficult.
%
% To help with this problem, the basis vectors of a FunctionSpace can be 
% rescaled easily, using the "rescale" method. This method takes a vector,
% with each element indicating the new scale value for the corresponding
% basis vector.

% Decompose the coordinate grid and calculate the standard deviation of the
% output columns

% 1) Use the original gaussian basis
clear Y stdY
YTmp = F.evaluate(xDecomp);
stdY{1} = std(YTmp);
Y{1} = YTmp;

% 2) Use the dim-reduced basis, without any scaling
FR = F.dimReduce(xDecomp);
FR.discardThresh(0.99);
YTmp = FR.evaluate(xDecomp);
stdY{2} = std(YTmp);
Y{2} = YTmp;

% 3) Use the dim-reduced basis, scaled to equalize variance across PCs
FR.rescale(1./stdY{2});
YTmp = FR.evaluate(xDecomp);
stdY{3} = std(YTmp);
Y{3} = YTmp;

figure(),
strs = {'Gaussian', 'Gaussian-PCA', 'Gaussian-PCA (scaled)'};

for n = 1:3
    
    subplot(2, 3, n);
    bar(stdY{n});
    ylabel('Std. dev');
    xlabel('x');
    title({'Std. dev.', strs{n}});
    
    subplot(2, 3, 3+n);
    imagesc(Y{n}');
    axis xy
    xlabel('x');
    ylabel('Component');
    cb = colorbar();
    ylabel(cb, 'Y');
    title({'Component values', strs{n}});
    set(gca, 'clim', prctile(Y{n}(:), [0.5, 99.5]));
end

%% Example 2b: repeat tuning curve fitting using DRFunctionSpace

% As before, we can use the "evaluate" method to decompose the predictor
% signal and reconstruct a tuning curve.
%
% We will do the fitting twice, first using the original scaling of the 
% PCs, and second after equalizing the variance.
%
% The results should show that the fitted tuning curve is unaffected by the
% scaling. This is because the FunctionSpace automatically applies the
% reverse scaling when reconstructing the tuning curve.
%
% The results should also show that the DR-fit for the tuning curve is less
% prone to overfitting problems, which are seen when using the original
% gaussian decomposition.

FR = F.dimReduce(xDecomp);
FR.discardThresh(0.99);

clear tCurveFitDR
for n = 1:2
    
    % Reset any previous scaling
    FR.resetScaling();
    
    % Equalize variance on the second round only
    if n==2
        FR.rescale(1./std(YD));
    end
    YD = FR.evaluate(sig);

    b = regress(yReg, [ones(nPoints, 1), YD]);
    c = b(1);  % intercept
    b(1) = [];

    % Combine the PC-basis vectors with the proportions given by the
    % regression fit.
    tmp = FR.evaluate(plottingGrid{1}, b');
    tCurveFitDR{n} = c + sum(tmp, 2);
end

figure()
clear h
plt = @(y, col, varargin) line(plottingGrid{1}, y, 'color', col, 'lineWidth', 2, varargin{:});

h(1) = plt(tuningFcn(plottingGrid{1}),'k');
h(2) = plt(tCurveFit,       'b');
h(3) = plt(tCurveFitDR{1},  'r', 'lineStyle', ':');
h(4) = plt(tCurveFitDR{2},  [0, 0.7, 0], 'lineStyle', '--');
legend(h, {'True', 'fit (all dims)', 'fit (DR)', 'fit (DR-scaled)'});
set(gca, 'ylim', [-2, 2]);