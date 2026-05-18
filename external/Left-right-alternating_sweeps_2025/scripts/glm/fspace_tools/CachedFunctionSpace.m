classdef CachedFunctionSpace < handle
    % Simplified FunctionSpace wrapper with cache, useful where basis
    % functions are expensive to evaluate.
    
    properties (SetAccess = protected)
        fsobj
        cache
    end
    
    properties (Dependent)
       ndims 
    end
    
    methods
        function self = CachedFunctionSpace(fsobj)
            self.fsobj = fsobj;
            self.cache = SimpleCache();
        end
        
        function Y = evaluate(self, X, p)
            if nargin < 3 || isempty(p), p = 1; end
            [Y, match] = self.cache.get(X);
            if ~match
                Y = self.fsobj.evaluate(X);
                self.cache.put(X, Y);
            end
            Y = Y .* p;
        end
        
        function D = decompose(self, X)
            [D, match] = self.cache.get(X);
            if match
                D = D.copy();
            else
                D = self.fsobj.decompose(X);
                self.cache.put(X, D);
            end
        end
        
        function val = get.ndims(self)
            val = self.fsobj.ndims;
        end
    end
end