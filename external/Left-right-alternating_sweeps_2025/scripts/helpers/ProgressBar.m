classdef ProgressBar < handle
    %PROGRESSBAR simple command-line progress bar
    
    properties (SetAccess = protected)
        stepSize
        trim
        lastStep = -inf
        lastStepStr = ''
        finished = false
    end
    
    properties
        format = '%.0f'
    end
    
    methods
        
        function self = ProgressBar(stepSize, trim)
            if nargin < 1 || isempty(stepSize), stepSize = 0.01; end
            if nargin < 2 || isempty(trim), trim = true; end
            self.stepSize = stepSize;
            self.trim = trim;
        end
        
        function update(self, status)
            
            isEndPoint = status==0 || status==1;
            isNewStep = status >= self.lastStep+self.stepSize;
            
            if isEndPoint || isNewStep
                percent = status * 100;
                str = sprintf([self.format, '%% '], percent);
                if self.trim
                    bstr = repmat(sprintf('\b'), 1, length(self.lastStepStr));
                    fprintf('%s', bstr);
                end
                fprintf('%s', str);
                self.lastStep = status;
                self.lastStepStr = str;
            end
            
            if status==1
                self.finish();
            end
        end
        
        function reset(self)
            self.lastStep = 0;
        end
        
        function finish(self)
            if ~self.finished
                fprintf('\n')
                self.finished = true;
            end
        end
    end
    
end

