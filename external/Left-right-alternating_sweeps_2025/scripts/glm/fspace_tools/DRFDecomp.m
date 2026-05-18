classdef DRFDecomp < FDecomp
    %DRFDECOMP dimension-reduced functional decomposition
    %
    % This adapatation of FDecomp adds functionality for discarding
    % dimensions.
    
    properties (SetAccess = protected)
        explained
    end
    
    properties (Dependent)
        dimInds
    end
    
    methods
        
        function self = DRFDecomp(varargin)
            self@FDecomp(varargin{:});
        end
        
        function indsKeep = discard(self, nKeep, mode, mult)
            if nargin < 4,  mult = []; end
            if nargin < 3, mode = ''; end
            indsKeep = self.fcnSpace.discard(nKeep, mode, mult);
            self.discardDims(indsKeep);
        end
        
        function indsKeep = discardThresh(self, thresh, mode, mult)
            if nargin < 4, mult = []; end
            indsKeep = self.fcnSpace.discardThresh(thresh, mode, mult);
            self.discardDims(indsKeep);
        end
        
        function val = get.explained(self)
            val = self.fcnSpace.explained;
        end
        
        function val = get.dimInds(self)
            val = self.fcnSpace.dimInds;
        end
        
    end
    
    methods (Access = protected)
        function discardDims(self, indsKeep)
            self.Y = self.Y(:, indsKeep);
        end
    end
    
end