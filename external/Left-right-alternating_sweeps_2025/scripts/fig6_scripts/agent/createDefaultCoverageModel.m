function mdl = createDefaultCoverageModel(P)

arguments
    P.modules = "off";
    P.footprintType {mustBeMember(P.footprintType, ["gaussian", "vm_invsqd"])} = "vm_invsqd";
    P.params = struct();
    P.posGridRange = [-2, 2];
    P.posGridStep = 0.01;
end

ft = P.footprintType;

mdl = ModularSweepAgentSimulation();
mdl.floatClass = "single";
mdl.sweepProfileType = ft;

if ft=="gaussian"
    sweepLenToSpacingRatio = 0.22; % empirically measured value
    mdl.sigma = 0.15 / sweepLenToSpacingRatio;
    % mdl.sigma = 0.3;               % width of sweep footprint/grid field
    widthArg = mdl.sigma;
    firstModuleSweepLength = 0.1;
elseif ft=="vm_invsqd"
    % sweep width as fraction of length
    mdl.kappa = 5;
    if P.modules=="off"
        firstModuleSweepLength = 0;
        widthArg = Inf;
    else
        widthArg = 0.5;
        firstModuleSweepLength = 0.1;
    end
end

scaleRatio = sqrt(2);
 
if strcmpi(P.modules, "off")
    mdl.configureModules(1, ft, firstModuleSweepLength, widthArg, scaleRatio);
else
    % Create a standard configuration of 3 modules
    mdl.configureModules(3, ft, firstModuleSweepLength, widthArg, scaleRatio);
    if strcmpi(P.modules, "independent")
        mdl.moduleIndependentSweepDirections = true;
    elseif strcmpi(P.modules, "coordinated")
        mdl.moduleIndependentSweepDirections = false;
    else
        error("Value of parameter 'modules' must be 'off', 'independent', or 'coordinated'.");
    end
end

mdl.createSquarePosGrid(P.posGridStep, P.posGridRange);

% finally, assign custom model params (these lower-level settings will 
% override any settings specified by other arguments)
for fd = fieldnamesstr(P.params)
    mdl.(fd) = P.params.(fd);
end

end