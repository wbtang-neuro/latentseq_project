classdef (Abstract) DRAlgorithm < matlab.mixin.Copyable
    
    properties
        explained
        nDims
        nDimsR
        dimIndsR
        initialized = false
        name
    end
    
    properties (Constant)
        MIN_EXPLAINED = 1e-8;
    end
    
    methods
        function self = DRAlgorithm(X)
            if nargin
                self.run(X);
            end
        end
        
        function run(self, X)
            self.nDims = size(X, 2);
            [self.nDimsR, self.explained] = self.onRun(X);
            self.dimIndsR = 1:self.nDimsR;
            self.initialized = true;
        end
        
        function X = transformFwd(self, X, dims)
            if nargin < 3 || isempty(dims), dims = (1:self.nDimsR)'; end
            X = self.onTransformFwd(X, dims);
        end
        
        function X = transformRev(self, XR, dims)
            if nargin < 4 || isempty(dims), dims = (1:self.nDimsR)'; end
            X = self.onTransformRev(XR, dims);
        end
        
        function str = toString(self)
            str = sprintf('%s, initialized=%u, nDims: original=%u, reduced=%u', ...
                class(self), self.initialized, self.nDims, self.nDimsR);
        end
        
        function indsKeep = discard(self, nKeep, mode, mult)
            if nargin < 3 || isempty(mode), mode = 'above'; end
            if nargin < 4 || isempty(mult), mult = 1; end
            nKeep = ceil(nKeep / mult) * mult;
            if strcmpi(mode, 'above') % discard dims above specified index
                indsKeep = 1:nKeep;
            elseif strcmpi(mode, 'below') % discard dims below
                nd = self.nDimsR;
                indsKeep = nd-nKeep+1 : nd;
            end
            self.discardDims(indsKeep);
        end
        
        function indsKeep = discardThresh(self, thresh, mode, mult)
            if nargin < 4, mult = []; end
            if nargin < 3, mode = 'above'; end
            cumExp = cumsum(self.explained);
            if strcmpi(mode, 'below')
                cumExp = flip(cumExp);
            end
            tol = self.MIN_EXPLAINED;
            n = find(cumExp+tol >= thresh, 1, 'first');
            indsKeep = self.discard(n, mode, mult);
        end
            
        function discardDims(self, indsKeep)
            self.onDiscardDimsR(indsKeep);
            self.explained = self.explained(indsKeep);
            self.nDimsR = numel(indsKeep);
            self.dimIndsR = self.dimIndsR(indsKeep);
        end
        
    end
    
    methods (Abstract)
        % Calculates the parameters of the dim-reduction. Must return the
        % number of dimensions in the reduced data matrix
        nDimsR = onRun(self, X);
        
        % Apply the dim-reduction to observations
        XR = onTransformFwd(self, X, dims);
        
        % Restore reduced observations to original basis
        X = onTransformRev(self, XR, dims);
        
        onDiscardDimsR(self, indsKeep);
    end
    
end