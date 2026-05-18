function glmPlots(fig, mdl, plotType, plotData)

handles = getPlotHandles(mdl, plotType);

if nargin == 4 && ~isempty(plotData)
    vinds = plotData.variableInds;
    vnames = string(fieldnames(vinds));
    nvars = numel(vnames);
else
    nvars = 0;
end

% scheme = figColScheme2;
% fgcol = scheme.fgCol;
% bgcol = scheme.bgCol;

switch lower(plotType)
        
    case "glm"
        
        clf(fig);
        nunits = mdl.nunits;
        beta = gather(mdl.F);
        nbeta = size(beta, 1);
        
        
        % GLM HDTC
        for n = 1:nvars
            vname = vnames(n);
            fspace = plotData.fspaces.(vname);
            X = fspace.defaultPlotX{1};
            z = fspace.evaluate(X, 1);
            
            ax = subplot(nvars+1, 1, n, "parent", fig);
            inds = plotData.variableInds.(vname);
            betaVar = beta(inds, :);
            zw = z*betaVar;
%             zw = exp(zw);
            title(ax, "GLM beta");
            imagesc(ax, 1:nunits, X, zw);
            xlabel(ax, "unit #");
            ylabel(ax, "HD / rad");
%             ax.CLim = [0, prctile(zw(:), 99)];
            ax.CLim = prctile(zw(:), [1, 99]);
            colorbar(ax);
            title(ax, vnames(n));
        end
        
        % GLM BETA VALS
        ax = subplot(3, 1, 3, "parent", fig);
        z = reshape(beta, nbeta, nunits);
        z = z(2:end, :);
        if ~isempty(z)
            imagesc(ax, z);
            xlabel(ax, "unit #");
            ylabel(ax, "param #");
            scl = prctile(abs(z(:)), 99);
            ax.CLim = scl*[-1, 1];
            title(ax, "BETA VALS");
            colorbar(ax);
        end
        
end

getPlotHandles(mdl, plotType, handles);


end