classdef (Abstract) LinearDRAlgorithm < DRAlgorithm
    
    % Base class for dimension reduction that can be expressed as a linear 
    % transformation
    
    properties
        useCentering = true;
        coeffs
        mu
    end
    
    methods (Abstract)
        C = onComputeCoeffs(Xc);
    end
    
    methods
        
        function [nDimsR, explained] = onRun(self, X)
            if self.useCentering
                self.mu = mean(X);
                X = bsxfun(@minus, X, self.mu);
            end
            [self.coeffs, explained] = self.onComputeCoeffs(X);
            nDimsR = size(self.coeffs, 2);
        end
        
        function X = onTransformFwd(self, X, dims)
            if self.useCentering
                X = bsxfun(@minus, X, self.mu);
            end
            X = X * self.coeffs(:, dims);
        end
        
        function X = onTransformRev(self, XR, dims)
            X = XR / self.coeffs(:, dims);
            if self.useCentering
                X = bsxfun(@plus, X, self.mu);
            end
        end
        
        function onDiscardDimsR(self, indsKeep)
            self.coeffs = self.coeffs(:, indsKeep);
        end
    end
    
end