classdef CompositeModel < BaseModel
    % TODO make this class independent of PoissonModel
    
    properties
        models       (:,1) cell
        modelUnitIds (:,1) cell
        logYh
        logYPredCache
    end
    
    properties (Dependent)
        nmodels
    end
    
    methods
        function self = CompositeModel(Y, models, allUnitIds)
            if nargin
                args = {Y};
            else
                args = {};
            end
            self@BaseModel(args{:});
            if nargin
                assert(numel(allUnitIds)==size(Y, 2), ...
                    "ID must be specified for each column of Y");
                self.unitIds = allUnitIds;
                nmodels = numel(models);
                for m = 1:nmodels
                    self.addModel(models{m});
                end
            end
        end
        
        function addModel(self, model)
            idx = self.nmodels+1;
            unitIds = model.unitIds; % can be logical, numeric or string
            if islogical(unitIds)
                unitIds = find(unitIds);
            end
            assert(model.nunits == numel(unitIds));
            if isempty(self.models)
                names = string([]);
            else
                names = cellfun(@(mdl) mdl.name, self.models);
            end
            assert(~ismember(model.name, names), "All models must have unique names");
            
            self.models{idx, 1} = model;
            self.modelUnitIds{idx, 1} = unitIds;
            model.unitIds = unitIds;
        end
        
        function [mdl, mdlIndex] = getModel(self, modelName)
            [mdlExists, mdlIndex] = self.hasModel(modelName);
            if mdlExists
                mdl = self.models{mdlIndex};
            else
                mdl = [];
            end
        end

        function [tf, idx] = hasModel(self, modelName)
            names = cellfun(@(mdl) mdl.name, self.models);
            modelName = string(modelName);
            idx = find(modelName == names);
            tf = ~isempty(idx);
        end
        
        function L = onStep(self, ~)
            % override
            
            mdls = self.models;
            n = self.nmodels;

            if self.display
                fprintf("iter=%03u\n", self.iter);
            end
            
            % step through models in series
            for m = 1:n

                logYPred = self.getLogYHat(m);
                
                % If any of the models reaches is stall iteration limit,
                % stop optimization
                mdl = mdls{m};
                L = mdl.step(logYPred);
                
                if self.display
                    str = mdl.iterStringLast;
                    fprintf("\t\t%s\n", str);
                end
            end
            
            if self.plotInterval~=0 && rem(self.iter, self.plotInterval) == 0
                self.drawPlots(true, false);
            end
            
        end
        
        function onInit(self)
            nmod = self.nmodels;
            for m = 1:nmod
                mdl = self.models{m};

                % This class will trigger the command-line display/plotting
                % operations during fitting. Therefore we disable this
                % functionality in the invididual model components, to
                % prevent double-calling these functions.
                mdl.display = false;
                mdl.plotInterval = 0; % stop models from triggering their own plotting
                mdl.initialize();
                
                if self.display
                    fprintf("Initialized model '%s', nunits=%u\n", mdl.name, mdl.nunits);
                end
            end
        end

        function str = onIterString(self)
            if self.display
                % for detailed output, this is separately handled during
                % overridden step()
                str = [];
                return;
            end
            str = sprintf("\t\t");
            nmod = self.nmodels;
            for m = 1:nmod
                mdl = self.models{m};
                str = str + mdl.iterString();
                if m < nmod
                    str = str + sprintf(",\t\t");
                end
            end
        end
        
        function logY = onPredictLogY(self, ~, ~)
            error("Not implemented for this class");
        end
        
        function logY = getLogYHat(self, modelId)
            % Calculate estimated Y input for a specified model

            inds = 1:self.nmodels;

            if isnumeric(modelId)
                % modelId is the numeric index
                modelNum = modelId;
            elseif ischar(modelId) || isstring(modelId)
                % modelId is the model's name string
                [~, modelNum] = self.getModel(modelId);
            else
                error("Invalid value for model ID");
            end
            inpModelInds = inds(inds~=modelNum);
            
            % Destination model and unit IDs
            mdlDest = self.models{modelNum};
            uidDest = self.modelUnitIds{modelNum};
            
            sz = [self.nt, mdlDest.nunits];

            % We accumulate the logYHat contributions, one model at a time.
            logY = zeros(sz, "like", self.Y);
            
            for i = 1:numel(inpModelInds)
                iinp = inpModelInds(i);
                uidSrc = self.modelUnitIds{iinp};
                logY = self.getOneLogYHat(iinp, uidSrc, uidDest, logY);  
            end
            
        end
        
        function Yh = getOneLogYHat(self, imdl, iusrc, iudest, Yh)
            ncol = numel(iudest);
            icol = zeros(1, ncol);
            mdl = self.models{imdl};
            [v, loc] = ismember(iudest, iusrc);
            icol(v) = loc(v);
            if nargin < 5
                Yh = 0;
            end
            Yh = Yh + mdl.predictLogY(icol);
        end
        
        function drawPlots(self, refresh, throwErr)
            if nargin < 2 || isempty(refresh), throwErr = true; end
            if nargin < 3 || isempty(throwErr), refresh = true; end
            for m = 1:self.nmodels
                hasPlots(m) = self.models{m}.drawPlots(false, throwErr);
            end
            if any(hasPlots) && refresh
                drawnow();
            end
        end
        
        function val = get.nmodels(self)
            val = numel(self.models);
        end
        
    end
    
end