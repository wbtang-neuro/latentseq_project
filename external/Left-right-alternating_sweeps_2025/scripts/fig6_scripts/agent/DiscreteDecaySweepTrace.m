classdef DiscreteDecaySweepTrace < SweepTrace
    % A sweep trace whose memory lasts a discrete number of time steps
    properties
        memoryDuration
        footprintBuffer
    end

    methods

        function self = DiscreteDecaySweepTrace(Finit, memoryDuration)
            self@SweepTrace(Finit);
            self.memoryDuration = memoryDuration;
            self.footprintBuffer = zeros([self.Fsize, memoryDuration], "like", Finit);
        end

        function self = addFootprint(self, footprint)
            % self.F = self.F + footprint;

            % Circularly shift the footprint buffer back one step and place
            % the new footprint at the beginning.
            self.footprintBuffer = circshift(self.footprintBuffer, 1, 3);
            self.footprintBuffer(:, :, 1) = footprint;

            % The updated trace F now consists of the current sum over the
            % footprint buffer
            self.F = sum(self.footprintBuffer, 3);
        end

        function self = timeStep(self, dt)
            % the passage of time has no effect on the trace (we only need
            % to count time steps), so we don't do anything here.
        end

    end



end