classdef PCA < LinearDRAlgorithm
    
    methods
        
        function self = PCA(varargin)
            self@LinearDRAlgorithm(varargin{:});
            self.name = 'PCA';
        end
        
        function [coeffs, explained] = onComputeCoeffs(self, X)
            n = size(X, 1);
            coeffs = X' * X ./ (n-1);
            [V, D] = eig(coeffs);
            [D, iSort] = sort(diag(D), 'descend');
            V = V(:, iSort);
            coeffs = V;
            explained = D/sum(D);
        end
    end
    
end