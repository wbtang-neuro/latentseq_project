classdef (Abstract) BaseModel < handle
    % Base class for a multi-variate Poisson spiking model

    properties

        name

        Y                       % spike trains     (t, dims)
        unitIds
        
        L = NaN
        
        % nf % TODO: reimplement as dependent prop; size(self.F, 1)
        
        floatClass
        
        % Iteration data
        iter  = 0                               % iteration counter
        niter = 500                             % iteration limit
        iterStringLast                          % command-line display string for previous iter
        
        status          = "uninitialized"
        
        minFuncOptions = defaultMinFuncOptions()
        
        % Fitting plots and command-line display
        display = true                          % enable/disable command-line display?
        catchPlottingErrors = false             % suppress errors in plotting functions?
        plotInterval = 10                       % refresh plots after this number of iterations

        version
    end

    properties (Constant)
        VERSION = 8
    end
    
    properties (Dependent)
        nunits
        nt
        nspikes
    end
    
    properties (Transient)
        plotFcns = {}
        figures             (:,1) matlab.ui.Figure
    end
    
    methods (Abstract)
        onInit(self)
        L = onStep(self)
        logy = onPredictLogY(self, icol, logy)
    end
    
    methods
        
        function self = BaseModel(varargin)
            if nargin
                Y = varargin{1};
                self.Y = Y;
                self.floatClass = underlyingType(Y);
            end
            self.version = self.VERSION;
        end
        
        function initialize(self)
            % Must be called before fitting, after all user-specified 
            % parameters have been set.
            assert(~isempty(self.Y), "Y cannot be empty");
            if self.status ~= "initialized"
                self.floatClass = class(gather(self.Y(1)));
            end
            
            % Initialize iteration variables
            self.onInit();
            self.iter = 0;
            self.status = "initialized";
        end
        
        function fit(self, niters)
            assert(self.status ~= "uninitialized")
            if nargin == 2
                % Optional second arg specifies iteration limit
                self.niter = self.iter + niters;
            end
            self.status = "fitting";
            while self.iter < self.niter
                self.step();
            end
            self.status = "finished fitting";
        end
        
        function L = step(self, inpLogR)
            if nargin < 2, inpLogR = []; end

            self.iter = self.iter + 1;
            self.onIterInit();

            L = self.onStep(inpLogR);
            self.L = L;

            self.iterStringLast = self.iterString();
            
            % do subclass cleanup stuff
            self.onIterFinish();
            
            % Display output
            if self.display
                fprintf("%s\n", self.iterStringLast);
            end
            
            if self.plotInterval~=0 && rem(self.iter, self.plotInterval) == 0
                self.drawPlots(true, true);
            end
            
        end
        
        function onIterInit(self)
        end
        
        function onIterFinish(self)
        end

        function yh = predictLogY(self, icol, yh)
            if nargin < 2, icol = 1:self.nunits; end
            if nargin < 3
                yh = 0;
            else
                % yh must either be scalar 0 or a matrix nt*numel(icol)
                yhsz = size(yh);
                ncol = numel(icol);
                assert(isequal(yh, 0) || isequal(yhsz, [self.nt, ncol]));
            end
            yh = yh + self.onPredictLogY(icol);
        end
        
        function hasPlots = drawPlots(self, refresh, catchErrors)
            if nargin < 3 || isempty(catchErrors), catchErrors = self.catchPlottingErrors; end
            if nargin < 2 || isempty(refresh), refresh = true; end
            nplts = numel(self.plotFcns);
            hasPlots = nplts > 0;
            for p = 1:nplts
                self.drawOnePlot(p, catchErrors);
            end
            if hasPlots && refresh
                drawnow();
            end
        end
        
        function drawOnePlot(self, n, catchErrors)
            if nargin < 3 || isempty(catchErrors), catchErrors = self.catchPlottingErrors; end
            plotFcn = self.plotFcns{n};
            fig = findFig(self, n);

            if isempty(fig)
                % Create new figure
                fig = figure("windowStyle", "docked");
                fig.NumberTitle = 'off';
                figTag = self.figureTag(n);
                fig.Name = figTag;
                fig.Tag = figTag;
            end
            if catchErrors
                try
                    plotFcn(fig, self);
                catch e
                    warning("Failed to draw plot: '%s'", e.message);
                end
            else
                plotFcn(fig, self);
            end
        end
        
        function str = nameString(self)
            if isempty(self.name)
                str = class(self);
            else
                str = self.name;
            end
        end
        
        function str = iterString(self)
            if self.iter == 0, str = ""; return, end
            
            % template method for any specialized subclass-defined output
            strCustom = self.onIterString();
            
            Lstr = sprintf("L=%.7g, ", self.L);
            Lstr = sprintf("%12s", Lstr);
            str = sprintf("%15s iter=%03u, %s%s", self.nameString(), ...
                self.iter, Lstr, strCustom);
        end

        function str = onIterString(self)
            str = "";
        end
        
        function str = figureTag(self, figNum)
            str = sprintf("%s_#%u", self.nameString(), figNum);
        end
        
        function fig = findFig(self, figNum)
            allFigs = findobj(0, "type", "figure");
            
            if isempty(allFigs)
                tags = string([]);
            else
                tags = string({allFigs.Tag});
            end
            
            qtag = self.figureTag(figNum);
            match = find(tags == qtag);
            if isempty(match)
                fig = [];
            else
                fig = allFigs(match(1));
            end
        end
        
        function filenames = saveFigures(self, baseFilename, figNums)
            if nargin < 3 || strcmpi(figNums, "all")
                figNums = 1:numel(self.plotFcns);
            end
            
            filenames = string([]);
            for f = 1:numel(figNums)
                fnum = figNums(f);
                fig = self.findFig(fnum);
                if isempty(baseFilename)
                    prefix = self.name;
                else
                    prefix = baseFilename;
                end
                fn = sprintf("%s_%u.fig", prefix, fnum);
                fprintf("Saving figure '%s' ... ", fn);
                try
                    saveas(fig, fn);
                    fprintf("done\n");
                    filenames(end+1) = fn;
                catch e
                    warning("failed to save figure '%s': '%s'", fn, e.message);
                end
            end
        end
        
        function A = castArr(self, A)
            A = cast(A, self.floatClass);
        end
        
        function set.unitIds(self, val)
            % TODO clean this up when Composite model doesn't depend on
            % this superclass
            if ~isempty(val) && isnumeric(val) && ~isempty(self.nunits)
                assert(numel(val)==self.nunits);
            end
           self.unitIds = val; 
        end
        
        function val = get.nunits(self)
            val = size(self.Y, 2);
        end
        
        function val = get.nt(self)
            val = size(self.Y, 1);
        end
        
        function val = get.nspikes(self)
            val = sum(self.Y);
        end
        
    end
    
end