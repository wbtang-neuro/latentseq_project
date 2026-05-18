classdef (Abstract) BaseFunctionSpace < handle
    %BASEFUNCTIONSPACE interface for function-space classes
    
    properties
        nDims
        dimScales
        name string = '<no name>'
    end
    
    properties (Abstract)
        plotFcn
        defaultPlotX
        defaultPlotAxesGrid
    end
    
    properties (Dependent)
        fcnNDimsIn
        fcnNDimsOut
    end
    
    methods (Abstract, Access=protected)
        % Evaluate a M*N matrix of function-domain data at one point in the
        % function space, where M is the number of function-domain data
        % points and N is the number of function input dimensions.
        % Parameter "p" can be a 1*P vector with each element representing
        % the value of one function basis vector, where P is the number of
        % basis functions. Alternatively, "p" can be a M*P matrix
        % with each row representing one point in the function space.
        %
        % Must return = M*P
        Y = evaluateImpl(self, X, p)
        
        % Evaluate function-domain data points for one basis function.
        y = evaluate1Impl(self, X, idx)
        
        n = fcnNDimsInImpl(self)
        n = fcnNDimsOutImpl(self)
        
        D = decompImpl(self, X)
    end
    
    methods
        
        function self = BaseFunctionSpace(nDims)
            if nargin
                self.nDims = nDims;
                self.dimScales = ones(1, nDims);
            end
        end
        
        function Y = evaluate(self, X, p)
            % EVALUATE - evaluate basis functions
            %
            % Y = S.EVALUATE(X) evaluates the basis functions at the
            % function-domain coordinates specified in N-by-D matrix X, 
            % where N is the number of observations and D is the 
            % dimensionality of the function domain. The returned matrix Y 
            % contains the value of each function for each point in X. Y 
            % has dimensions [N, B], where B is the number of basis 
            % functions.
            %
            % Y = S.EVALUATE(X, P) evaluates the basis functions with
            % weightings specified by vector P, which defines a point in
            % the function space. P can be a scalar, in which case all
            % basis functions will be weighted equally, a row vector
            % of length B (to specify a single), or a N-by-B matrix
            
            nd = self.nDims;
            if nargin < 3 || isempty(p), p = ones(1, 'like', X); end
            if isscalar(p)
                p = p .* ones(1, nd);
            end
            Y = self.evaluateImpl(X, p);
            Y = self.scaleForward(Y);
        end
        
        function y = evaluate1(self, X, idx)
            % Evaluate one function at one point in the space
            y = self.evaluate1Impl(X, idx);
            y = self.scaleForward(y, idx);
        end
        
        function D = decompose(self, X)
            % DECOMPOSE perform functional decomposition
            % 
            % D = F.DECOMPOSE(X) evaluates the value of each of the basis
            % functions on input data matrix X, returning the decomposition
            % in FDecomp object D.
            D = self.decompImpl(X);
        end
        
        function fcn = createEvalFunction(self, p)
            fcn = @(X) self.evaluate(X, p);
        end
        
        function rescale(self, scales, relative)
            % RESCALE apply scaling to basis vectors
            %
            % F.RESCALE(S), where S is a row vector containing one element
            % for each basis vectors, sets the scales of the Nth basis
            % function to S(N). S may also be a scalar, in which case S is
            % used to set the scale for all basis vectors.
            %
            % F.RESCALE(S, REL) uses logical flag REL to indicate whether
            % the scaling is relative or not. If REL is true, the new basis
            % vector scales will be the product of the previous scale and
            % the scaling factor specified by S.
            if nargin < 3 || isempty(relative), relative = false; end
            if isscalar(scales), scales = scales .* ones(1, self.nDims); end
            assert(isequal(size(scales), [1, self.nDims]));
            if relative
                scales = self.scaleForward(scales);
            end
            self.dimScales = scales;
        end
        
%         function wrap(self, bFcn)
%             % WRAPINNER wrap the FunctionSpace's BFs around another BF
%             for b = 1:self.nDims
%                 self.basisFcns(b).wrap(bFcn);
%             end
%         end
        
        function resetScaling(self)
            self.rescale(1);
        end
        
        function tf = isScaled(self, inds)
            if nargin < 2, inds = 1:self.nDims; end
            tf = any(self.dimScales(inds) ~= 1);
        end
        
        function val = get.fcnNDimsIn(self)
            val = self.fcnNDimsInImpl;
        end
        
        function val = get.fcnNDimsOut(self)
            val = self.fcnNDimsOutImpl;
        end
        
        function [h, ax] = plot(self, varargin)
            % PLOT(S, P, COORDS, ARGS)
            % Plot evaluation for a point in func-space
            [h, ax] = self.plotSub('one', varargin{:});
            if ~nargout, clear h; end
        end
        
        function [h, ax] = plotAll(self, varargin)
            % PLOTALL(S, P, COORDS, ARGS)
            % Plot individual function evaluations for a point in func-space
            [h, ax] = self.plotSub('all', 1, varargin{:});
            if ~nargout, clear h; end
        end
        
    end
    
    methods (Hidden)
        
        function yS = scaleForward(self, y, inds)
            if nargin < 3, inds = 1:self.nDims; end
            if self.isScaled(inds)
                yS = y .* self.dimScales(inds);
            else
                yS = y;
            end
        end
        
        function yS = scaleReverse(self, y, inds)
            if nargin < 3, inds = 1:self.nDims; end
            if self.isScaled(inds)
                yS = y ./ self.dimScales(inds);
            else
                yS = y;
            end
        end
        
        function [h, axs] = plotSub(self, plotMode, p, varargin)
            % Plot individual function evaluations for a point in func-space
            if nargin < 3 || isempty(p), p = 1; end
            
            inp = inputParser();
            inp.addParameter('figure', []);
            inp.addParameter('showLabels', false);
            inp.addParameter('showNames', false);
            inp.addParameter('axes', []);
            inp.addParameter('plottingCoords', self.defaultPlotX)
            inp.addParameter('plottingAxesGrid', self.defaultPlotAxesGrid);
            inp.addParameter('transformFcn', []);
            inp.parse(varargin{:});
            P = inp.Results;
            
            if isempty(P.figure), P.figure = gcf(); end
            
            % Get the grid of function input values
            [~, gridMat, gridSz] = self.parseCoordsGrid(P.plottingCoords{:});
            
            % Get the grid for axis labelling
            if isempty(P.plottingAxesGrid)
                P.plottingAxesGrid = P.plottingCoords;
            end
            gridCoordsAx = self.parseCoordsGrid(P.plottingAxesGrid{:});
            
            singlePlot = strcmp(plotMode, 'one');
            
            if singlePlot
                z = sum(self.evaluate(gridMat, p), 2);
                z = reshape(z, gridSz);
                nPlots = 1;
                if isempty(P.axes)
                    ax = gca();
                else
                    ax = P.axes;
                end
            elseif strcmp(plotMode, 'all')
                nx = ceil(sqrt(self.nDims));
                ny = ceil(self.nDims/nx);
                layout = tiledlayout(ny, nx, "TileSpacing", "tight", "parent", P.figure);
                zAll = self.evaluate(gridMat, p);
                nPlots = self.nDims;
            end
            
            for d = 1:nPlots
                if ~singlePlot
                    ax = nexttile(layout);
%                     ax = subplot(ny, nx, d, 'parent', P.figure);
                    z = zAll(:, d);
                    z = reshape(z, gridSz);
                end
                if ~isempty(P.transformFcn)
                    z = feval(P.transformFcn, z);
                end
                h(d, 1) = self.plotFcn(ax, gridCoordsAx{:}, z);
                axs(d, 1) = ax;
                if ~singlePlot && P.showNames
                    title(ax, "BF #" + string(d));
                    % title(ax, self.basisFcns(d).name); % not valid for DRFunctionSpace
                end
            end
            if ~P.showLabels, axis(axs, 'off'); end

            if ~nargout, clear h; end
            
        end
        
        function [grid, fullGridMat, gridSz] = parseCoordsGrid(self, varargin)
            grid = varargin;
            nd = self.fcnNDimsIn;
            if numel(grid) == 1 && nd > 1
                grid = repmat(grid, 1, nd);
            end
            if nd > 1
                [tmp{1:nd}] = meshgrid(grid{:});
                gridSz= size(tmp{1});
                tmp = cellfun(@(x) {x(:)}, tmp);
                fullGridMat = [tmp{:}];
            else
                fullGridMat = grid{1};
                gridSz = size(grid{1});
            end
            
            
        end
        
    end
    
    
end