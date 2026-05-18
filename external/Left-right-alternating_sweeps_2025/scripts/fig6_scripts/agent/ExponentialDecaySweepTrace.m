classdef ExponentialDecaySweepTrace < SweepTrace
    % A sweep trace that decays exponentially over time

    properties
        tau = 1
    end

    methods

        function self = ExponentialDecaySweepTrace(Finit, tau)
            self@SweepTrace(Finit);
            self.tau = tau;
        end

        
        function self = addFootprint(self, footprint)
            self.F = self.F + footprint;
        end

        function self = timeStep(self, dt)
            k = self.tau^dt;
            self.F = (self.F*k);
            self.F = self.F * self.tau;
        end
    end

end