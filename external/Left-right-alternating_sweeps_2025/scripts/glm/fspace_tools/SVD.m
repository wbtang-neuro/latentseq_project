classdef SVD < LinearDRAlgorithm
    
    % X = U*S*V'
    %
    % This maps X <---> U
    
    properties (SetAccess = protected, Transient)
        svdRaw
    end
    
    methods
        
        function self = SVD(varargin)
            self@LinearDRAlgorithm(varargin{:});
            self.name = 'SVD';
        end
        
        function [coeffs, explained] = onComputeCoeffs(self, X)
            [U, S, V] = svd(X, 'econ');
            coeffs = pinv(S*V');
            DOF = size(X, 1) - 1;
            latent = diag(S).^2./DOF;
            explained = latent ./ sum(latent);
            self.svdRaw = struct('U', U, 'S', S, 'V', V);
        end
        
    end
    
end