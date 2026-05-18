classdef PoissonGLM < TuningModel

    properties
        X                           % GLM design matrix (variables x units)
        useIntercept = true         % treat column as an intercept term; this means it won't be penalized
        optimizeIntercept = true    % if true, intercept will be learned during fitting (
        lambdaL1 = 0;               % L1-regularization term
    end
    
    methods
        
        function self = PoissonGLM(varargin)
            self@TuningModel(varargin{:});
            if nargin
                X = varargin{2};
                self.X = X;
                self.nf = size(X, 2);
            end
        end
        
        function onInit(self)
            assert(~issparse(self.Y));
            self.F = randn(self.nf, self.nunits, "like", self.Y) * 1e-5;
            if self.useIntercept
                self.F(1, :) = log(mean(self.Y));
            end
        end

        function yh = onPredictLogY(self, iunit)
            % Predict log-yhat for all or a subset of units
            if nargin < 2 || strcmpi(iunit, "all")
                iunit = 1:self.nunits;
            end
            b = self.F(:, iunit);
            yh = self.X*b;
        end
        
        function [L, Fnew] = onOptimizeF(self, F, inpR)
            X = self.X;
            Y = self.Y(:);
            inpR = inpR(:);
            
            lam1 = self.parseLambda(self.lambdaL1);
            
            % If using intercept and it's enabled
            if self.useIntercept
                fi = self.F(1, :);
            else
                fi = [];
            end
            
            if self.optimizeIntercept
                % this case is simple, because all betas are optimized
                fcn = @(f) poissonLlh(X, Y, f, -1, inpR, lam1);
            else
                % this case is a bit more complex, because only a subset of
                % betas are optimized
                fcn = @(f) costfcn(X, Y, f, fi, inpR, lam1);
            end
            [Fnew, L] = minFunc(fcn, F(:), self.minFuncOptions);

            % Fit optimized subset of F into the complete F matrix
            Fnew = reshape(Fnew, size(self.F));
        end        
        
        function lam = parseLambda(self, lam0)
            lam = lam0;
            if isscalar(lam)
                lam = lam + zeros(self.nf, self.nunits);
            end 
            if self.useIntercept
                % don't penalize the intercept
                lam(1, :) = 0;
            end
        end

    end
    
end

function varargout = costfcn(X, Y, b, bi, inpR, lam1)
% X    - GLM design matrix
% Y    - nt*nunits spike count matrix
% b    - beta values to be optimized
% bi   - beta values for intercept (not optimized; optional)
% lam1 - L1 lambda values (

useIntercept = ~isempty(bi);
if useIntercept
    nunits = numel(bi);
    b = reshape(b, [], nunits);
    nc = size(b, 1);
    b = [bi; b];
end
b = b(:);

[varargout{1:nargout}] = poissonLlh(X, Y, b, -1, inpR, lam1);

if nargout == 2 && useIntercept
    % We'll get gradients for the intercept terms too, which we'll discard
    b = reshape(varargout{2}, 2, []);
    iint = 1 : (nc+1) : numel(b);
    b(iint) = [];
    b = b(:);
    varargout{2} = b;
end


end