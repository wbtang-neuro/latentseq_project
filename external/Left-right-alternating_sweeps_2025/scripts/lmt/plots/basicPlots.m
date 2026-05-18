function basicPlots(fig, mdl)

persistent mdl0
newmdl = isempty(mdl0) || mdl0~=mdl;

H = getPlotHandles(mdl, "basic");
if isfield(H, "axes")
    axs = H.axes;
else
    axs = [];
end

% idat = mdl.iterData;
idat = {struct()};

F = mdl.F;
Fdat = {F};
names = "F";
% FLast = mdl.stateLastIter.F;
% Fdat = {F, FLast-F};
% names = ["F", "dF"];

for n = 1:numel(Fdat)
    newax = ~checkHandle(axs, "image", n) || newmdl;
    z = gather(Fdat{n});
    x = 1:mdl.nunits;
    y = 1:size(z, 1);
    
    if isprop(mdl, "isCircular") && mdl.isCircular
        % special case for a circular 1-D tuning curves
        z = [z; z];
        y = linspace(0, 4*pi, size(z, 1));
        yt = [0, 2, 4] * pi;
        ytlab = ["0", "2\pi", "4\pi"];
    else
        yt = [];
        ytlab = [];
    end
    
    if newax
        ax = subplot(5, 2, n, "parent", fig);
        axs.image(n) = ax;
        H.images(n) = imagesc(ax, x, y, z);
        colorbar(ax);
        axis(ax, 'xy');
        ax.YTick = yt;
        ax.YTickLabels = ytlab;
        title(ax, names(n));
    else
        h = H.images(n);
        h.CData = z;
        ax = h.Parent;
%         v = ~isoutlier(z(:));
        v = ~isnan(z);
        clim = prctile(z(v), [1, 99]);
        clim = max(abs(clim)) * [-1, 1]; % centered around zero
        if clim(2) > clim(1)
            ax.CLim = clim;
        end
    end

    
end

names=["Fopt", "Fstd", "FstdUnit_m", "FL", "FL_poiss_m", "FL_gp", "FL_L1", ...
    "FdL_poiss", "FdL_gp", "FdL_L1", ...
    "Fstepsize", "FstepsizeZ", ...
    "Xopt", "Xstd", "XL", "XL_poiss_m", "XL_poiss_ratio", "XL_prior", ...
    "Xstepsize", "XstepsizeZ"];
nfds = numel(names);

for n = 1:nfds
    
    name = names(n);
    titlestr = name;
    
    if endsWith(name, "_z")
        normfcn = @zscore;
        name = char(name);
        name = string(name(1:end-2));
        zl = [-3, 3];
    elseif endsWith(name, "_m")
        normfcn = @(x) x./mean(x);
        name = char(name);
        name = string(name(1:end-2));
        zl = [0, 2];
    else
        normfcn = [];
    end
    
    if isfield(idat{end}, name) && ~all(isnan(idat{end}.(name)(:)))
        
        hasfd = cellfun(@(s) isfield(s, name), idat);
        if all(hasfd)
            z1 = idat{1}.(name);
            % If the data is a scalar of vector, concatenate values across all steps
            if min(size(z1)) == 1
                z = cellfun(@(s) {double(s.(name)(:)')}, idat);
                z = cat(1, z{:});
                xname = "step";
                yname = "";
            else
                % If the data is a matrix, plot only for the current step
                z = idat{end}.(name)';
                xname = "unit #";
                yname = "F bin";
            end
            z = gather(z);
        else
            continue;
        end
        
        if isempty(normfcn)
            rng = prctile(z(:), [5, 95]);
            zl = mean(rng) + 0.75*diff(rng)*[-1, 1];
%             zl(1) = max(0, zl(1)); % assume no quantities can be negative
        else
            z = normfcn(z);
        end
        
        
        [np, nz] = size(z);
        if nz > 1
            ndz = 2;
        else
            ndz = 1;
        end
        
        if ndz == 1
            pltType = "line";
        elseif ndz == 2
            pltType = "image";
        end
        newplt = ~checkHandle(H, pltType, n);
        
        if newplt
            ax = subplot(5, 5, 5 + n, "parent", fig);
            cla(ax);
            ax.Title.String = titlestr;
            ax.Title.Interpreter = "none";
            if pltType == "line"
                h = line(ax, nan, nan, "color", 'k');
            elseif pltType == "image"
                h = imagesc(ax);
                axis(ax, "xy");
            end
            xlabel(ax, xname);
            ylabel(ax, yname);
            H.(pltType)(n) = h;
        else
            h = H.(pltType)(n);
            ax = h.Parent;
        end
        
        if pltType == "line"
            h.XData = (1:np)';
            h.YData = z; 
            if zl(2) > zl(1)
                ax.YLim = zl;
            end
        elseif pltType == "image"
            h.XData = 1:np;
            h.CData = z';
            ax = h.Parent;
            if zl(2) > zl(1)
                ax.CLim = zl;
            end
            colorbar(ax);
        end
        
    end
end

H.axes = axs;
getPlotHandles(mdl, "basic", H);

mdl0 = mdl;

end