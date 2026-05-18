classdef FDecomp < matlab.mixin.Copyable
    %FDECOMP functional decomposition
    
    properties (SetAccess = protected)
        X           % input data    (points*dims)
        Y           % output data   (points*basis-functions)
        nPoints     % number of data points
        gridded     % 
    end
    
    properties (Access = protected)
        % Keep the fcnSpace restricted, to prevent direct manipulation of
        % parameters, which could conflict with a FDecomp object's own 
        % properties.
        
        fcnSpace    % FunctionSpace obj used for decomposition
    end
    
    properties (Hidden)
        XGrid       % Stores input coordinate grid, in case where gridded = true
    end
    
    properties (Dependent)
        nDimsX      % Dimensionality of input
        nDimsY      % Dimensionality of output
        scalesY     % Scaling of function basis vectors
    end
    
    methods
       
        function self = FDecomp(S, varargin)
            if nargin
                self.fcnSpace = S.copy();
                X = varargin;
                nX = numel(X);
                
                % Input data "X" can be supplied either as a single M*N
                % matrix (M points, in N-D function domain), or as N
                % vectors, specifying a coordinate grid.
                %
                % Determine whether the input data represents a grid
                self.gridded = nX > 1 || size(X{1}, 1) == 1;
                if self.gridded
                    self.XGrid = X;
                    [Xg{1:nX}] = ndgrid(X{:});
                    X = cat(nX+1, Xg{:});
                    sz = size(X);
                    X = reshape(X, prod(sz(1:nX)), nX);
                else
                    X = X{1};
                end
                self.nPoints = size(X, 1);
                self.X = X;
                self.Y = self.fcnSpace.evaluate(X);
            end
        end
        
        function rescale(self, scale, relative)
            % RESCALE change basis vector scaling
            % 
            % Only relative scaling will be applied, so convert absolute
            % scale factors to relative ones.
            if relative
                scaleRel = scale;
            else
                scaleRel = scale ./ self.scalesY;
            end
            self.fcnSpace.rescale(scaleRel, true);
            self.Y = self.Y .* scaleRel;
        end
        
        function scaleFactor = rescaleBy(self, scaleFcn, scaleConst, scaleMode, relative)
            % Rescale by result of a specified scaling function
            if nargin < 5 || isempty(relative), relative = false; end
            if nargin < 4 || isempty(scaleMode), scaleMode = 'inverse'; end
            if nargin < 3 || isempty(scaleConst), scaleConst = 1; end
            
            fscale = scaleFcn(self.Y);
            if strcmpi(scaleMode, 'inverse')
                fscale = 1./fscale;
            elseif ~strcmpi(scaleMode, 'linear')
                error('Argument ''scaleMode'', must be either ''normal'' or ''inverse''');
            end
            scaleFactor = scaleConst .* fscale;
            self.rescale(scaleFactor, relative);
        end
        
        function resetScaling(self)
            scale0 = self.scalesY;
            self.Y = self.Y ./ scale0;
            self.fcnSpace.rescale(1);
        end
        
        function Yi = interpolate(self, Xi, method, extrapVal)
            %INTERPOLATE interpolate a decomposition
            %
            % YI = FD.INTERPOLATE(XI) linearly interpolates the basis
            % function values at the coordinates specified in XI. XI must
            % be a N*D matrix, where N is the number of points and D is the
            % dimenisonality of the function domain.
            
            if nargin < 3 || isempty(method), method = 'linear'; end
            if nargin < 4, extrapVal = NaN; end
            nInterp = size(Xi, 1);
            Yi = zeros(nInterp, self.nDimsY);
            if self.gridded
                grid = self.XGrid;
                gridSz = cellfun(@numel, grid);
                if strcmpi(method, 'nearest')
                    validInterp = true(nInterp, 1);
                    for d = 1:self.nDimsFcnIn
                        [subs{d}, v] = gridnn(grid{d}, Xi(:, d));
                        validInterp = validInterp & v;
                    end
                    inds = sub2ind(gridSz, subs{:});
                    Yi = self.XR(inds, :);
                    if ~all(validInterp)
                        Yi(~validInterp, :) = extrapVal;
                    end
                else
                    [Xg{1:self.nDimsX}] = ndgrid(grid{:});
                    Xi = num2cell(Xi, 1);
                    gridSz = size(Xg{1});
                    for d = 1:self.nDimsY
                        Y =  reshape(self.Y(:, d), gridSz);
                        Yi(:, d) = interpn(Xg{:}, Y, Xi{:}, method, extrapVal);
                    end
                end
            else
                for d = 1:self.nDims
                    scint = scatteredInterpolant(self.X, self.Y(:, d), method);
                    Yi(:, d) = scint(Xi);
                end
            end
        end
        
        function tf = isScaled(self)
            tf = self.fcnSpace.isScaled();
        end
        
        function val = get.nDimsX(self)
            val = self.fcnSpace.fcnNDimsIn;
        end
        
        function val = get.nDimsY(self)
            val = self.fcnSpace.nDims;
        end
        
        function val = get.scalesY(self)
            val = self.fcnSpace.dimScales;
        end
        
    end
    
    methods (Access = protected)
        function cpself = copyElement(self)
            cpself = self.copyElement@matlab.mixin.Copyable();
            cpself.fcnSpace = self.fcnSpace.copy();
        end
    end
    
end