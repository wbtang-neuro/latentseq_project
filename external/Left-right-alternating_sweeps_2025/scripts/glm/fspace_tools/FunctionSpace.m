classdef FunctionSpace < BaseFunctionSpace & matlab.mixin.Copyable
    %FUNCTIONSPACE a space defined by a set of function basis vectors
   
    properties
        basisFcns
        plotFcn
        defaultPlotX = {};
        defaultPlotAxesGrid = {};
    end
    
    methods
        
        function self = FunctionSpace(varargin)
            % FUNCTIONSPACE constructor
            %
            % F = FUNCTIONSPACE() creates an empty FunctionSpace object F.
            %
            % F = FUNCTIONSPACE(BFUNCS) creates FunctionSpace object F,
            % using the array of BasisFunction objects BFUNCS as basis
            % vectors.
            if nargin
                bFuncs = varargin{1};
                nDims = numel(bFuncs);
                args = {nDims};
            else
                args = {};
            end
            self@BaseFunctionSpace(args{:});
            self.basisFcns = bFuncs;
        end
        
        function [DRFS, D] = dimReduce(self, X, alg)
            % DIMREDUCE generate dimension-reduced version of FunctionSpace
            % 
            % FR = F.DIMREDUCE(X) uses data matrix X to generate
            % DRFunctionSpace object FR, which represents a PCA-transformed
            % version of F.
            
            if nargin < 3 || isempty(alg), alg = PCA(); end
            DRFS = DRFunctionSpace(self, alg);
            DRFS.name = sprintf('%s -> %s', self.name, alg.name);
            D = DRFS.runDR(X);
            DRFS.plotFcn = self.plotFcn;
            DRFS.defaultPlotX = self.defaultPlotX;
            DRFS.defaultPlotAxesGrid = self.defaultPlotAxesGrid;
        end
        
        function str = toString(self)
            str = sprintf('%s %s, %u dimensions', class(self), self.name, self.nDims);
        end
        
    end
    
    methods (Access = protected)
        
        function cpSelf = copyElement(self)
            cpSelf = copyElement@matlab.mixin.Copyable(self);
            cpSelf.basisFcns = self.basisFcns.copy();
        end
        
        function val = nDimsImpl(self)
            val = numel(self.basisFcns);
        end
        
        function val = fcnNDimsInImpl(self)
            val = self.basisFcns(1).nDimsIn;
        end
        
        function val = fcnNDimsOutImpl(self)
            val = self.basisFcns(1).nDimsOut;
        end
        
        function Y = evaluateImpl(self, X, p)
            % Evaluate function-domain data points
            % dims are (points, nbas)
            Y = zeros(size(X, 1), self.nDims, 'like', X);
            for d = 1:self.nDims
                ytmp = p(:, d) .* self.evaluate1Impl(X, d);
                Y(:, d) = cast(ytmp, 'like', X); % older matlab versions don't allow implicit GPU -> CPU cast
            end
        end
        
        function y = evaluate1Impl(self, X, idx)
            % Evaluate one element
            y = self.basisFcns(idx).evaluate(X);
        end
        
        function D = decompImpl(self, X)
            D = FDecomp(self, X);
        end
        
    end
    
end