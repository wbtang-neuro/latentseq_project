classdef SweepTrace

    properties
        F
        Fsize
    end

    methods (Abstract)
        % These functions must be implemented by the subclass, to specify
        % how new sweep footprints are integrated, and how the trace
        % evolves over time.
        self = addFootprint(self, footprint)
        self = timeStep(self, dt)
    end

    methods
        function self = SweepTrace(Finit)
            self.F = Finit;
            self.Fsize = size(Finit);
        end
    end

end