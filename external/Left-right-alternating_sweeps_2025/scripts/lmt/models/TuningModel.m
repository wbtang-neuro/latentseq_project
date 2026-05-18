classdef (Abstract) TuningModel < BaseModel
    % Abstract class adding tuning functionality

    properties
        nf                          % number of elements in tuning curve
        F                           % tuning curve parameters (dims nf x nunits)
        enableFStep = true          % enable optimization of F?
    end

    methods (Abstract)
        [L, Fnew] = onOptimizeF(self, F, inpLogR)
    end

    methods

        function self = TuningModel(varargin)
            self@BaseModel(varargin{:});
        end

        % Implement the step function
        function L = onStep(self, inpLogR)
            if self.enableFStep
                [L, self.F] = self.onOptimizeF(self.F, inpLogR);
            else
                L = self.L;
            end
        end

    end

end