classdef ModularSweepAgentSimulation < matlab.mixin.Copyable

    % To be decided: what to call the weights/profile of a single sweep?
    % currently "profile" is used, but the variables are named "ww".

    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % User-configurable params
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % These two parameters define the 2D position grid that covers the
        % agent's ambient space. The sweep profile and history trace will
        % be computed at each position in this grid.
        %
        % N.B. this is a dimensionless model, so the positions and
        % distances don't have units; it's only the relative scaling that
        % matters.
        posGridX = -2 : 0.01 : 2
        posGridY = -2 : 0.01 : 2

        % The inverse-distance weighting of 2D position bins means that as
        % distance approaches zero, the weighting approaches infinity. This
        % causes the model to become unstable. To mitigate this effect, we
        % ignore bins that are very near when integrating the sweep trace.
        % The property "minDist" defines the threshold radius below which
        % bins are ignored.
        minDist = 0.05

        % Grid module settings: here we define a set of spatial "modules"
        % analogous to grid-cell modules. Each module generates sweeps with 
        % a characteristic radius.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sweep footprint settings

        % sweepProfileType = "vm_invsqd" % von-mises / inverse-square distance
        sweepProfileType = "gaussian"

        % The sweep's angular profile is given by the Von-Mises
        % distribution, with a concentration parameter "kappa". Higher
        % values will produce narrower sweeps.
        kappa = 5

        % To run without separate modules, configure 'nModules' to 1 and
        % set the scale properties to NaN
        nModules = 1
        moduleScaleRatio = sqrt(2)
        firstModuleSweepLength = 0.10
        sigma = 0.3                     % TODO: use this more consistently

        moduleIndependentSweepDirections = false

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Trace forgetting: the sweep history trace decays over time, with
        % two different modes available ("exponential" or "discrete")

        sweepTraceDecayMode = "exponential"

        % When using the "exponential" mode, the trace decays by a fixed 
        % proportion with each unit time that passes. The rate of forgetting
        % is controlled by the parameter "tau", which ranges from 0 (
        % instant forgetting) to 1 (no forgetting).
        tau = 1

        % If using the "discrete" forgetting mode, the trace consists of
        % the equally weighted sum of the N most recent footprints. N is
        % set by the parameter "sweepTraceDiscreteN".
        sweepTraceDiscreteN = 1


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Other misc settings

        % When using the "run" function, this option determines whether the
        % final sweep in the sequence is added to the sweep trace.
        addFinalSweepToTrace = true

        ditherSweepDirections = true % generates non-discrete sweep directions by adding noise

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting
        plotType = "none"
        plotScale = 5
        plotSweepHistoryLength = 10
        contourLevels = 0
        plotAllContours = false
        traceNormalization = "trace"
        clipPercentile = 99.8
        arrowLength = 0

        moduleColorMaps

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Internal data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        initialized = false;

        % Spatial grids
        xx                  % vector of x locations in 2D position grid 
        yy                  % vector of y locations in 2D position grid 
        posBinWidth         % width of the square x,y bins
        dirGrid             % 

        % Sweep "profile" (the spatial weights of individual sweeps)
        lastSweepProfileM    % (matrix) profile of the previous sweep, separated by module

        profileTotalWeightM

        profileMaxWeightM    % maximum weight value of a sweep (for normalization; not used by agent)
        profileMaxWeight

        % Sweep "trace" (the sum of the profiles of past sweeps)
        sweepTrace                      % 
        sweepTraceM                     % sweep traces separated by module
        traceUpdateMode = "parallel"    % either "parallel" or "serial"

        castDataFcn

        % Plotting
        axes
        plotHandles

        useGpu
        floatClass = "single"

    end

    properties (Constant)
        DEFAULT_COLORMAP = bone(256)
    end

    properties (Dependent)
        doPlot
        posGridSize
        F
        FM
        dispersionParam
        dispersion
    end

    methods

        function self = ModularSweepAgentSimulation()
            % Determine whether machine has CUDA-capable GPU. If so, we'll
            % use it.
            try
                gpuDevice(1);
                self.useGpu = true;
            catch e
                self.useGpu = false;
            end
        end

        function configureModules(self, n, profileType, firstModuleSweepLength, fieldWidth, scaleRatio)
            assert(isnumeric(n)&&n>=0, "Argument 'n' should be a non-negative integer");
            if n==0
                n = 1;
                % the remaining arguments are ignored
                profileType = "vm_invsqd";
                firstModuleSweepLength = NaN;
                fieldWidth = NaN;
                scaleRatio = NaN;
            end
            self.sweepProfileType = profileType;
            self.nModules = n;
            self.moduleScaleRatio = scaleRatio;
            self.firstModuleSweepLength = firstModuleSweepLength;
            self.sigma = fieldWidth;
        end

        function initialize(self)
            % Call this function after configuring the model

            if self.useGpu
                self.castDataFcn = @(x) gpuArray(cast(x, self.floatClass));
            else
                self.castDataFcn = @(x) cast(x, self.floatClass);
            end
            self.initSpatialGrids();

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Modules
            
            % Calculate normalizing variables
            ww = 0;
            for m = 1:self.nModules
                wwM = self.calcModuleSweepProfile(0, 0, 0, m, false, true);

                % normalize module footprint so it sums to 1
                p = gather(sum(wwM));
                wwM = wwM/p;

                % add normalized footprint to total
                ww = ww+wwM;

                self.profileTotalWeightM(m) = p; % store normalizer
                self.profileMaxWeightM(m) = gather(max(wwM)); % store max value *after* normalization
            end

            self.profileMaxWeight = gather(max(ww)); % max of normalized total

            % Set up colormaps. Use 'bone' for a single module, and
            % independent R,G,B channels for up to 3 modules. More than 3
            % modules will require other solutions
            ncmap = size(self.DEFAULT_COLORMAP, 1);
            if self.nModules==1
                self.moduleColorMaps{1} = self.DEFAULT_COLORMAP;
            elseif self.nModules <= 3
                for m = 1:self.nModules

                    % Map modules to R,G,B channels
                    col = [0, 0, 0];
                    col(m) = 1;
                    cmap = col .* linspace(0, 1, ncmap)';

                    % use colormap from Settings
                    % cmap = S.col_grid_module(m) .* linspace(0, 1, ncmap)';
                    self.moduleColorMaps{m} = cmap;
                end
            end

            if self.doPlot
                self.plotHandles = self.initPlot();
            end

            self.initialized = true;
            
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Internal methods

        function initSpatialGrids(self)
            % Set up the binned environment
            gx = self.castData(self.posGridX);
            gy = self.castData(self.posGridY);
            [xx,yy] = meshgrid(gx, gy);
            self.xx = xx(:);
            self.yy = yy(:);
            bwx = gx(2)-gx(1);
            bwy = gy(2)-gy(1);
            assert(bwx==bwy, "posGridX and posGridY must have the same bin spacing");
            self.posBinWidth = bwx;

            % Set up trace arrays ("F" variables)
            F0 = zeros(prod(self.posGridSize), 1);
            F0 = self.castData(F0);
            F0M = zeros(prod(self.posGridSize), self.nModules);
            F0M = self.castData(F0M);

            if strcmpi(self.sweepTraceDecayMode, "exponential")
                createTraceFcn = @(f) ExponentialDecaySweepTrace(f, self.tau);
            elseif strcmpi(self.sweepTraceDecayMode, "discrete")
                createTraceFcn = @(f) DiscreteDecaySweepTrace(f, self.sweepTraceDiscreteN);
            else
                error("Parameter 'sweepTraceDecayMode' must be set to 'exponential' or 'discrete'");
            end

            % self.F = F0;
            self.sweepTrace = createTraceFcn(F0);

            % For analysis purposes, also store the individual module
            % last-sweep profile and whole history trace
            self.lastSweepProfileM = F0M;
            % self.FM = F0M;
            self.sweepTraceM = createTraceFcn(F0M);

            S = SweepsSettings();
            self.dirGrid = self.castData(S.gv.angular);
        end

        function createSquarePosGrid(self, step, range)
            % Convenience function for configuring a square position grid
            % centered on (0, 0). Call initialize() afterwards to complete
            % the setup.
            if isscalar(range)
                gv = -range : step : range;
            elseif numel(range)==2
                gv = range(1) : step : range(2);
            end
            self.posGridX = gv;
            self.posGridY = gv;
        end

        function A = castData(self, A)
            A = self.castDataFcn(A);
        end

        function [wwQ, dd, aa] = calcModuleSweepProfile(self, xp, yp, dirQ, moduleId, normalize, useMinDist)
            % INPUTS
            % xp, yp: The agent's current position (the origin of the
            %         sweeps)
            %
            % dirQ:     Vector of hypothetical sweep directions for which weights
            %         will be calculated
            % 
            % OUTPUTS
            % wwaq:   (Matrix) sweep profile weights, with rows 
            %         corresponding to the x,y-position grid and columns 
            %         corresponding to the input sweep directions dirQ.
            %
            % aa:     (Vector) offset direction of each position in the x,y
            %         grid with respect to the agent's position (xp, yp)
            % dd:     (Vector) 

            if nargin < 6 || isempty(normalize), normalize = true; end
            if nargin < 7 || isempty(useMinDist), useMinDist = true; end

            dirQ = dirQ(:)';

            ptype = self.sweepProfileType;

            sweepLength = self.firstModuleSweepLength * (self.moduleScaleRatio^(moduleId-1));

            if ptype=="vm_invsqd"
                [dd, aa] = self.calcPosGridOffsets(xp, yp);

                % Define the footprint density function. This is a product
                % of the two functions below:
                % (1) the Von-Mises function:           exp(k.*cos(a1-a2))
                % (2) the inverse-square distance:      1 ./ (d.^2)

                fcn = @(k, a1, a2, d) exp(k.*cos(a1-a2)) ./ (d.^2);

                if self.useGpu
                    % When using GPU, it's faster to call the density
                    % function elementwise via arrayfun
                    wwQ = arrayfun(fcn, self.kappa, aa, dirQ, dd);
                else
                    % On CPU, calling the function directly is faster
                    wwQ = fcn(self.kappa, aa, dirQ, dd);
                end

                % apply module radius thresholds
                f = (1 + self.sigma);
                rIn = sweepLength / f;
                rOut = sweepLength * f;
                if rIn>0 || rOut<Inf
                    invalid = dd<rIn | dd>=rOut;
                    wwQ(invalid, :) = 0;
                end
            elseif ptype=="gaussian"
                [u, v] = pol2cart(dirQ, sweepLength);
                dd = self.calcPosGridOffsets(xp+u, yp+v);
                aa = []; % don't bother to compute this
                f = self.sigma * sweepLength;
                wwQ = normpdf(dd, 0, f) ./ (sweepLength.^2);
            end

            if normalize
                assert(self.initialized, "Model is not initialized");
                wwQ = wwQ ./ self.profileTotalWeightM(moduleId);
            end

            if useMinDist && self.sweepProfileType=="vm_invsqd"
                wwQ = self.applyMinDist(dd, wwQ);
            end

        end

        function [ww, dd, aa] = calcTotalSweepProfile(self, xp, yp, dirQ, useMinDist)
            if nargin < 5 || isempty(useMinDist), useMinDist = true; end
            ww = 0;
            for m = 1:self.nModules
                [wwaM, dd, aa] = self.calcModuleSweepProfile(xp, yp, dirQ, m, true, useMinDist);
                ww = ww+wwaM;
            end
        end

        function [dd, aa] = calcPosGridOffsets(self, xp, yp)
            % Calculate the distance and angular offset between the agent and each
            % position bin in the environment
            dxx = self.xx-xp;
            dyy = self.yy-yp;
            dd = hypot(dxx, dyy);

            % don't allow distances less than one bin width
            dd = max(dd, self.posBinWidth);

            if nargout==2
               aa = atan2(dyy, dxx); % we only need aa with this footprint type
            end
        end

        function zza = applyMinDist(self, dd, zza)
            % this version needs DD to be single-column. ZZA can have
            % multiple columns.
            invalid = dd<self.minDist;
            if size(dd,2)==1
                zza(invalid, :) = 0;
            elseif size(dd,2)==size(zza,2)
                zza(invalid) = 0;
            else
                error("Inputs dd and zza have incompatible sizes");
            end
        end

        function [dirSweep, dirOptim, scoresAll] = chooseActualAndOptimalSweepDir(self, wwDirQ, dirSweep0)
            % Find optimal sweep direction, given an array of footprint
            % weights for all possible sweep directions
            scoresAll = sum(wwDirQ .* self.F);
            [~, iBestScore] = min(scoresAll);
            dirOptim = gather(self.dirGrid(iBestScore));
            dAng = gather(diff(self.dirGrid([1, 2])));
            if self.ditherSweepDirections
                dirOptim = dirOptim + (rand()-0.5)*dAng; % dither
            end
            [dirSweep, dirOptim] = parseSweepDir(dirSweep0, dirOptim);
        end

        function [dirOptim, scores, scoresM] = step(self, xp, yp, dt, dirSweep0, updateTrace)

            assert(self.initialized, "Model is not initialized");

            % Apply temporal updating to the sweep traces
            self.sweepTrace = self.sweepTrace.timeStep(dt);
            self.sweepTraceM = self.sweepTraceM.timeStep(dt);

            % % Exponential decay of existing trace
            % k = self.tau^dt;
            % self.F = (self.F*k);

            % Allocate output arrays. If "sweepDir" is "none" on this
            % step, these outputs will remain blank.
            blankM = nan(1, self.nModules);
            if self.moduleIndependentSweepDirections
                blank = blankM;
            else
                blank = nan;
            end
            dirActual = blank;
            dirOptim = blank;

            % Allocate scores array. This contains the overlap scores used
            % to choose the optimal sweep direction. If modules choose
            % their own sweep directions, it contains module-specific
            % scores. If a single direction is used for all modules, the
            % scores will be a sum across modules.

            na = numel(self.dirGrid);
            if self.moduleIndependentSweepDirections
                scores = nan([na, self.nModules]);
            else
                scores = nan([na, 1]);
            end
            scoresM = nan([na, self.nModules]);

            imposeSweepDirection = isnumeric(dirSweep0);

            if ~strcmpi(dirSweep0, "none")
                % Method 1: each module picks its own optimal sweep dir
                if self.moduleIndependentSweepDirections
                    for m = 1:self.nModules
                        [dirActual(m), dirOptim(m), scores(:, m)] = self.doSweepChoiceStep(xp, yp, dirSweep0, m, updateTrace);
                    end
                    % We already have the modulewise scores
                    scoresM = scores;
                else
                    % Method 2: choose sweep direction for total footprint
                    [dirActual, dirOptim, scores] = self.doSweepChoiceStep(xp, yp, dirSweep0, [], updateTrace);
                    % also calculate module-wise scores
                    for m = 1:self.nModules
                        wwQ = self.calcModuleSweepProfile(xp, yp, self.dirGrid, m);
                        scoresM(:, m) = sum(wwQ .* self.F);
                    end
                end

                if updateTrace
                    % update individual module traces
                    % self.FM = self.FM + self.lastSweepProfileM;
                    self.sweepTraceM = self.sweepTraceM.addFootprint(self.lastSweepProfileM);
                    if self.traceUpdateMode=="parallel"
                        % update main trace (modules in parallel)
                        self.sweepTrace = self.sweepTrace.addFootprint(sum(self.lastSweepProfileM, 2));
                        % self.F = self.F + sum(self.lastSweepProfileM, 2);
                    end
                end
            end

            if self.doPlot
                self.updatePlot(xp, yp, dirOptim, dirActual, imposeSweepDirection);
            end
        end

        function applyTraceForgetting(self)
            
        end


        function [dirSweep, dirOptim, scoresAll] = doSweepChoiceStep(self, xp, yp, dirSweep, moduleId, updateTrace)
            % Perform a sweep-choosing step
            useTotalFootprint = isempty(moduleId);

            if useTotalFootprint
                wwQ = self.calcTotalSweepProfile(xp, yp, self.dirGrid);
                imod = 1:self.nModules;
            else
                wwQ = self.calcModuleSweepProfile(xp, yp, self.dirGrid, moduleId);
                imod = moduleId;
            end

            % Use the supplied footprint weights to calculate the optimal
            % sweep direction
            [dirSweep, dirOptim, scoresAll] = self.chooseActualAndOptimalSweepDir(wwQ, dirSweep);

            % Compute single-module profiles of the selected sweep dir
            for m = imod
                self.lastSweepProfileM(:, m) = self.calcModuleSweepProfile(xp, yp, dirSweep, m, true, false);
                if updateTrace && self.traceUpdateMode=="serial"
                    % self.F = self.F + self.lastSweepProfileM(:, m);
                    self.sweepTrace = self.sweepTrace.addFootprint(self.lastSweepProfileM(:, m));
                end
            end

            % % Compute single-module profiles of the selected sweep dir
            % if useTotalFootprint
            %     for m = 1:self.nModules
            %         self.lastSweepProfileM(:, m) = self.calcModuleSweepProfile(xp, yp, dirSweep, m, true, false);
            %     end
            % else
            %     self.lastSweepProfileM(:, moduleId) = self.calcModuleSweepProfile(xp, yp, dirSweep, moduleId, true, false);
            % end

        end

        function [optimalSweepDirs, scores, scoresM, F, profileM] = run(self, xp, yp, dt, sweepInds, imposedSweepDirs)

            if nargin < 4 || isempty(dt), dt = 1; end % dt isn't used if tau==1
            if nargin < 5 || isempty(sweepInds), sweepInds = 1:numel(xp); end
            if nargin < 6, imposedSweepDirs = []; end

            % Check inputs
            nt = numel(xp);
            ndir = numel(self.dirGrid);

            if isscalar(dt), dt = repmat(dt, nt, 1); end
            doSweep = false(nt, 1);
            doSweep(sweepInds) = true;
            assert(all(sweepInds>=1));
            assert(all(sweepInds<=nt));
            csweep = 0;
            nsweeps = numel(sweepInds);

            imposeSweepDirs = isnumeric(imposedSweepDirs) && ~isempty(imposedSweepDirs);

            % allocate output arrays
            if self.moduleIndependentSweepDirections
                optimalSweepDirs = zeros(nt, self.nModules);
            else
                optimalSweepDirs = zeros(nt, 1);
            end

            Fsz = self.posGridSize;
            posGridSzM = [Fsz, self.nModules];

            if nargout >= 3
                sz = size(optimalSweepDirs);
                scores = zeros([sz(1), ndir, sz(2)], self.floatClass);
                scoresM = zeros(nt, ndir, self.nModules);
                F = zeros([nt, Fsz], self.floatClass);
                profileM = zeros([nt, posGridSzM], self.floatClass);
                outputAll = true;
            else
                outputAll = false;
            end

            % pb = ProgressBar();

            for i = 1:nt
                if doSweep(i)
                    csweep = csweep + 1;
                    if imposeSweepDirs
                        sweepDir = imposedSweepDirs(csweep);
                    else
                        if csweep == 1
                            % % when trace is blank, all directions are
                            % equally optimal
                            sweepDir = "first";
                        else
                            sweepDir = "optimal";
                        end
                    end
                    % pb.update(i/nt);
                else
                    sweepDir = "none";
                end

                % At the last sweep, we won't update the trace,
                % meaning that the model's final state will reflect the
                % moment just before the final sweep is executed.
                if self.addFinalSweepToTrace
                    addToTrace = true;
                else
                    addToTrace = csweep<nsweeps;
                end

                [   optimalSweepDirs(i,:), ...
                    scores(i, :, :), ...
                    scoresM(i, :, :) ] = self.step(xp(i), yp(i), dt(i), sweepDir, addToTrace);

                if outputAll
                    F(i, :, :, :) = reshape(self.F, [1, Fsz]);
                    profileM(i, :, :, :) = reshape(self.lastSweepProfileM, [1, posGridSzM]);
                end

            end
            
        end

        function h = initPlot(self)
            ax = self.axes;
            if isempty(ax)
                ax = gca;
            end
            h = struct();
            h.axes = ax;
            gx = self.posGridX;
            gy = self.posGridY;

            S = SweepsSettings();

            if self.plotType=="basic"
                h.trace = image(ax, gx, gy, nan(self.posGridSize));
            elseif self.plotType=="basic_mono"
                colormap(ax, self.DEFAULT_COLORMAP);
                h.trace = imagesc(ax, gx, gy, nan(self.posGridSize));
            end

            scale = self.plotScale;
            
            % Make contour plot for current sweep profile. We only do this
            % if the 'plotAllContours' option is disabled, because the
            % latter covers the current sweep anyway.
            if self.contourLevels && ~self.plotAllContours
                h.contour = self.makeFootprintContourPlot(ax);
            end
            h.contourAll = {};

            h.path = plot3(ax, nan, nan, nan, "color", S.col_pos_true, "lineWidth", scale/3);
            % h.pos = plot3(ax, nan, nan, nan, ".", "color", S.col_pos_true, "markerSize", scale*5);
            % h.posPast = plot3(ax, nan, nan, nan, ".", "color", S.col_pos_true, "markerSize", scale*5);

            if self.moduleIndependentSweepDirections
                n = self.nModules;
                cols = cellfun(@(x) {x(end, :)}, self.moduleColorMaps);
            else
                n = 1;
                cols = {S.col_covmodel};
            end

            for i = 1:n
                h.sweep.sim(i) = plot3(ax, nan, nan, nan, "color", cols{i}, 'lineWidth', scale/3);
            end
            h.sweep.real = plot3(ax, nan, nan, nan, "color", S.col_id, 'lineWidth', scale/3);

            blank = nan(1, self.plotSweepHistoryLength*3);
            for i = 1:n
                h.sweepPast.sim(i) = plot3(ax, blank, blank, blank, "color", cols{i}, 'lineWidth', scale/3);
            end
            h.sweepPast.real = plot3(ax, blank, blank, blank, "color", S.col_id, 'lineWidth', scale/3);

            axis(ax, "equal", "xy", "tight", "off");
            % ax.Colormap = 1-gray();
        end

        function h = makeFootprintContourPlot(self, ax, h)
            % helper function for creating individual-module or combined
            % footprint plot
            if nargin < 3, h = []; end
            patchArgs = {"faceAlpha", 0.3, "lineWidth", self.plotScale/3, "edgeAlpha", 1};

            hasObj = ~isempty(h);
            if ~hasObj
                clear h
            end

            if self.plotType=="basic"
                for m = 1 : self.nModules
                    if hasObj
                        htmp = h(m);
                    else
                        htmp = [];
                    end
                    col = self.moduleColorMaps{m}(end, :);
                    z = self.lastSweepProfileM(:, m);
                    z = z./self.profileMaxWeightM(m);
                    stackLevel = self.nModules - m + 1;
                    h(m) = self.contourPatchHelper(ax, htmp, z, stackLevel, patchArgs{:}, "edgeColor", col, "faceColor", col);
                end
            elseif self.plotType=="basic_mono"
                if hasObj
                    htmp = h;
                else
                    htmp = [];
                end
                z = sum(self.lastSweepProfileM, 2);
                z = z./self.profileMaxWeight;
                col = self.DEFAULT_COLORMAP(end, :);
                h = self.contourPatchHelper(ax, htmp, z, 1, patchArgs{:}, "edgeColor", col, "faceColor", col);
            end
        end

        function h = contourPatchHelper(self, ax, h, zdata, stackLevel, varargin)
            % Creates a new contourf plot, or updates an existing one if
            % the handle is supplied
            patchArgs = varargin;
            zdata = reshape(double(zdata), self.posGridSize);
            c = contourc(self.posGridX, self.posGridY, zdata, [1, 1]*self.contourLevels);
            x = c(1, 2:end)';
            y = c(2, 2:end)';
            z = ones(size(x))*stackLevel;
            if isempty(h)
                h = patch(ax, 'xdata', x, 'ydata', y, 'zdata', z, patchArgs{:});
            else
                % Patch already exists: just update its position
                h.XData = x;
                h.YData = y;
                h.ZData = z;
            end
            % zdata = reshape(zdata, self.posGridSize);
            % h.ZData = zdata;
        end


        % function h = contourfHelper(self, h, zdata, varargin)
        %     % Creates a new contourf plot, or updates an existing one if
        %     % the handle is supplied
        %     createArgs = varargin;
        %     if isempty(h)
        %         [c, h] = contourf(createArgs{:});
        %     end
        %     zdata = reshape(zdata, self.posGridSize);
        %     h.ZData = zdata;
        % end

        function F = getNormalizedTraceForPlotting(self, plotType)
            % Generate correctly normalized sweep trace, scaled to range
            % 0-1.

            traceProp = "F";
            normProp = "profileMaxWeight";
            if strcmpi(plotType, "basic")
                traceProp = traceProp + "M";
                normProp = normProp + "M";
            end

            F = gather(self.(traceProp));

            normPrm = self.traceNormalization;
            % Normalize by percentile of either the trace or the last
            % footprint

            if isnumeric(normPrm)
                F = F./normPrm;
            else
                if strcmpi(normPrm, "trace")
                    % use current trace
                    N = F;
                elseif strcmpi(normPrm, "footprint")
                    % use single-footprint density
                    N = self.(normProp);
                end
                p = prctile(N, self.clipPercentile);
                F = F./p;
            end

            % if strcmpi(normType, "off")
            %     % normalize using a precalculated constant and clip and
            %     % values that exceed 1
            %     p = self.(normProp);
            %     F = F ./ p;
            %     F(F>1) = 1;
            % else

                % F = F./prctile(N, self.clipPercentile);
            % end
        end

        function zrgb = getNormalizedTraceRgb(self, plotType)
            F = self.getNormalizedTraceForPlotting(plotType);
            % scale to color map and clip
            ncmap = size(self.DEFAULT_COLORMAP, 1);
            x = linspace(0, 1, ncmap);

            % % Plot the sweep traces (stored separately for each module)
            if strcmpi(plotType, "basic")
                zrgb = zeros(prod(self.posGridSize), 3);
                for m = 1:self.nModules
                    cmap = self.moduleColorMaps{m};
                    zq = F(:, m);
                    zrgb = zrgb + interp1(x, cmap, zq, "nearest", "extrap");
                end
            elseif strcmpi(plotType, "basic_mono")
                zq = F;
                zrgb = interp1(x, self.DEFAULT_COLORMAP, zq, "nearest", "extrap");
            end
        end

        function updatePlot(self, x, y, dirOptim, dirActual, imposeSweepDirection)

            h = self.plotHandles;
            ax = h.axes;

            h.pos.XData = x;
            h.pos.YData = y;

            if self.plotType == "basic"
                % with RGB image there's no axis clipping of color axis, so
                % we need to clip the data before plotting
                z = self.getNormalizedTraceRgb(self.plotType);
                z = reshape(z, [self.posGridSize, 3]);
                z(z<0) = 0;
                z(z>1) = 1;
            elseif self.plotType == "basic_mono"
                % no clipping needed when using imagesc
                z = self.getNormalizedTraceForPlotting(self.plotType);
                z = reshape(z, self.posGridSize);
            end

            h.trace.CData = double(sqrt(z));

            if any(self.contourLevels)
                if self.plotAllContours
                    % add new contour
                    if ~any(isnan(dirActual))
                        h.contourAll{end+1} = self.makeFootprintContourPlot(ax, []);
                    end
                else
                    % update last-sweep contour
                    h.contour = self.makeFootprintContourPlot(ax, h.contour);
                end

            end

            % Update path
            h.path.XData(end+1) = x;
            h.path.YData(end+1) = y;
            h.path.ZData(end+1) = 1;

            % Update sweep line plots
            if ~isnan(dirOptim)

                % % sweep position
                % h.posPast.XData(end+1) = x;
                % h.posPast.YData(end+1) = y;
                % h.posPast.ZData(end+1) = 1;

                % find maximum in sweep footprint
                z = self.lastSweepProfileM;
                if ~self.moduleIndependentSweepDirections
                    z = sum(z, 2);
                end
                % z = sum(self.lastSweepProfileM, 2);

                dd = self.calcPosGridOffsets(x, y);
                z = self.applyMinDist(dd, z);

                for m = 1:size(z, 2)

                    if strcmpi(self.arrowLength, "peak")
                        % arrows mark footprint peak positions
                        [~, imx] = max(z(:, m));
                        [i,j] = ind2sub(self.posGridSize, imx);
                        xarr = self.posGridX(j);
                        yarr = self.posGridY(i);
                    elseif self.arrowLength > 0
                        % Last part of string indicates arrow length, e.g.
                        % "direction1.5" means the arrow is 1.5 units long.
                        [u, v] = pol2cart(dirOptim(m), self.arrowLength);
                        xarr = x+u;
                        yarr = y+v;
                    else
                        xarr = NaN;
                        yarr = NaN;
                    end

                    xs = [x, xarr, nan];
                    ys = [y, yarr, nan];

                    zlevel = self.nModules - m + 1.5;

                    % current sweep
                    fd = "sim";
                    htmp = h.sweep.sim(m);
                    htmp.XData = xs;
                    htmp.YData = ys;
                    htmp.ZData = ones(size(xs)) * zlevel;

                    % sweep history
                    htmp = h.sweepPast.(fd)(m);
                    htmp.XData = [htmp.XData(4:end), xs];
                    htmp.YData = [htmp.YData(4:end), ys];
                    htmp.ZData = [htmp.ZData(4:end), ones(size(xs))*zlevel];
                end

            end

            self.plotHandles = h;
        end

        function val = get.doPlot(self)
            val = ~strcmpi(self.plotType, "none");
        end

        function val = get.posGridSize(self)
            val = [numel(self.posGridY), numel(self.posGridX)];
        end

        function val = get.F(self)
            val = self.sweepTrace.F;
        end
        
        function val = get.FM(self)
            val = self.sweepTraceM.F;
        end

    end
end

function [dirSweep, dirOptim] = parseSweepDir(dirSweep, dirOptim)
if strcmpi(dirSweep, "optimal")
    % self-driving mode: use the optimal direction
    dirSweep = dirOptim;
elseif strcmpi(dirSweep, "first")
    % self-driving mode, special case for first sweep:
    % pick a random direction
    dirSweep = 2*pi*rand();
    dirOptim = dirSweep;
else
    % driven mode: use supplied sweep direction
    assert(isnumeric(dirSweep), "Argument 'dirSweep' must be 'optimal', 'first' or a numeric value");
end
end