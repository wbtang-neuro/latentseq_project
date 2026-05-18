classdef KFoldCV
    %KFOLDCV k-fold crossvalidation set with interleaving
    
    properties
        nFolds
        nPoints
        div = 1    % number of subdivisions of train/test blocks
        foldInds
    end
    
    methods
        function self = KFoldCV(nPoints, nFolds, div)
            self.nPoints = nPoints;
            self.nFolds = nFolds;
            if nargin == 3
                if numel(div)==1
                    div = ceil(linspace(0, nPoints, div+1));
                end
                self.div = div;
            end
            self = self.repartition();
        end
        
        function self = repartition(self)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for n = 1:self.nDiv()
                divInds = self.div(n)+1 : (self.div(n+1));
                divLen = numel(divInds);
                divFoldInds = ceil((1:divLen)/divLen*self.nFolds);
                self.foldInds(divInds) = divFoldInds;
            end
        end
        
        function v = train(self, foldIdx)
            v = ~self.test(foldIdx);
        end
        
        function v = test(self, foldIdx)
            v = false(self.nPoints, 1);
            v(self.foldInds == foldIdx) = true;
        end
        
        function n = nDiv(self)
            n = numel(self.div)-1;
        end
    end
end

