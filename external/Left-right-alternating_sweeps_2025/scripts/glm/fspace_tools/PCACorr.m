classdef PCACorr < LinearDRAlgorithm
    % PCA based on correlation instead of covariance
    
    methods
    
        function self = PCACorr(varargin)
            self@LinearDRAlgorithm(varargin{:});
            self.name = 'PCACorr';
        end
        
        function [coeffs, explained] = onComputeCoeffs(self, Xc)
            n = vecnorm(Xc);
            R = (Xc' * Xc) ./ (n'*n);
            [V, D] = eig(R);
            [D, iSort] = sort(diag(D), 'descend');
            V = V(:, iSort);
            coeffs = V;
            explained = D/sum(D);
        end
        
    end
    
end