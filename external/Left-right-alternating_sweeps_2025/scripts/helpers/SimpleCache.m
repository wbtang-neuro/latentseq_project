classdef SimpleCache < matlab.mixin.Copyable
   
    properties
        compareFcn = @isequal
        compareAllFcn = []
        keyNotFoundFcn = []
    end
    
    properties (SetAccess = protected)
        keys    (:,1) cell = {}
        values  (:,1) cell = {}
    end
    
    properties (Dependent)
        nkeys
    end
    
    methods
        
        function clear(self)
            self.keys = {};
            self.values = {};
        end
        
        function put(self, key, val, replace)
            if nargin < 4 || isempty(replace)
                replace = true;
            end
            [match, idx] = self.iskey(key);
            if ~match
                idx = self.nkeys + 1;
            end
            self.keys{idx, 1} = key;
            self.values{idx, 1} = val;
        end
        
        function [val, match] = get(self, key)
            
            if self.nkeys
                if isempty(self.compareAllFcn)
                    [match, idx] = self.iskey(key);
                else
                    [match, idx] = self.compareAllFcn(key, self.keys);
                end
            else
                match = false;
            end
            
            if match
                val = self.values{idx};
            else
                val = self.onKeyNotFound();
            end
            
        end
        
        function [tf, idx] = iskey(self, key)
            tf = false;
            idx = [];
            fcn = self.compareFcn;
            k = 0;
            while k < self.nkeys
                k = k+1;
                try
                    if fcn(key, self.keys{k})
                        tf = true;
                        idx = k;
                        break;
                    end
                catch e
                    warning("Failed while comparing key #%u: '%s'. Removing this key.", e.message);
                    self.keys(k) = [];
                    self.values(k) = [];
                end
            end
        end
        
        function val = get.nkeys(self)
            val = numel(self.keys);
        end
    end
    
    methods (Access = protected)
        function val = onKeyNotFound(self)
            fcn = self.keyNotFoundFcn;
            if isempty(fcn)
                val = [];
            else
                val = feval(fcn);
            end
        end
    end
    
    methods (Static)
        function obj = getCache(cacheName)
            persistent cachecache
            if isempty(cachecache)
                cachecache = SimpleCache();
                cachecache.keyNotFoundFcn = @SimpleCache;
            end
            obj = cachecache.get(cacheName);
            cachecache.put(cacheName, obj);
        end
    end
    
end