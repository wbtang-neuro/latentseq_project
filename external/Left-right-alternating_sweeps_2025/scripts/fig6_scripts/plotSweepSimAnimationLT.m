function plotSweepSimAnimationLT(posxy, mdl, P)
% Animation showing agent decision-making algorithm

arguments
    posxy
    mdl
    P.inputSweepDirs = []
    P.videoFile = []
    P.frameRate = 30
    P.animSpeed = 2
end

xp = posxy(:,1);
yp = posxy(:,2);

saveVideo = ~isempty(P.videoFile);

frameRate = 30;
animSpeed = 2;

nplt = 15;
cl = [0, 1];
zmax = 1; % sweep trace clipping threshold
scale = 1;

selfDriving = isempty(P.inputSweepDirs);

mdl = mdl.copy();
mdl.initialize();

Fsz = mdl.posGridSize;

dirGrid = gather(mdl.dirGrid);
nDirGrid = numel(dirGrid);
[~, dirGridOrder] = sort(wrapTo2Pi(dirGrid)); % sort such that the first scanned angle is 0 (East)
dirGrid = dirGrid([1:end, 1]); % wrap for plotting

figure("colormap", bone());
ax = gca;
himg = imagesc(mdl.posGridX, mdl.posGridY, nan(Fsz));
F = zeros(Fsz);
hpath = plot(xp(1:nplt), yp(1:nplt), "r--", "lineWidth", 2);
htitle = text(xp(1), 0.3*scale, "", ...
    "fontSize", 14, ...
    "color", "w", ...
    "horizontalAlignment", "center", ...
    "verticalAlignment", "top");
axis("image", "off", scale*[-1.6000    1   -0.5000    0.5000]);

zblank = zeros(Fsz);
ztrace = repmat(zblank, [1, 1, 3]);

clear vwriter
if saveVideo
    vwriter = VideoWriter("animation.mp4", "MPEG-4");
    vwriter.FrameRate = P.frameRate;
    nPauseFrames = ceil(P.frameRate/P.animSpeed);
    nScanFrames = ceil(1/P.animSpeed);
    vwriter.open();
end

colMdl = [0.3, 0.3, 1];
colId = [0, 0.8, 0];

polrad = 0.1*scale;
[cx, cy] = pol2cart(dirGrid, polrad);
col = [0, 0, 0]+0.6;
hpolc = plot(cx, cy, "color", col, "lineWidth", 1.5);
hpolx = plot([-1, 1]*polrad, [0, 0], "color", col, "lineWidth", 1.5);
hpoly = plot([0, 0], [-1, 1]*polrad, "color", col, "lineWidth", 1.5);
hscores = plot(nan, nan, 'w', 'lineWidth', 2);
hscoremx = plot(nan, nan, 'o', 'color', colMdl, 'lineWidth', 2, 'markerSize', 10);
hinputDir = plot(nan, nan, 'o', 'color', colId, 'lineWidth', 2, 'markerSize', 10);

for i = 1:nplt

    hpolc.XData = cx + xp(i);
    hpolc.YData = cy + yp(i);
    hpolx.XData = xp(i) + [-1, 1]*polrad;
    hpoly.XData = xp(i)*[1, 1];
    htitle.Position(1) = xp(i);
    htitle.String = "NEW TIME STEP";
    if saveVideo
        frame = getframe(ax);
        for n = 1:nPauseFrames
            vwriter.writeVideo(frame);
        end
    else
        pause(1);
    end

    if i>1 || selfDriving

        % iterate through all possible sweep angles, calculating overlaps
        overlaps = zeros(nDirGrid, 1);
        wwIter = zeros([Fsz, nDirGrid]);
        for a = 1:nDirGrid
            ww = mdl.calcTotalSweepProfile(xp(i), yp(i), mdl.dirGrid(a), false);
            ww = gather(ww);
            ww = ww ./ mdl.profileMaxWeight;
            ww = reshape(ww, Fsz);
            % ww = reshape(ww, Fsz) / mdl.profileMaxWeight;
            wwIter(:, :, a) = ww;
            overlap = sum(ww(:) .* F(:));
            overlaps(a) = overlap;
        end

        r = overlaps;
        [~, imn] = min(overlaps);
        optimAngle = gather(mdl.dirGrid(imn));
        r = r./max(r); % max 1
        r = r * polrad;
        r = r([1:end, 1]);
        [u, v] = pol2cart(dirGrid, r);
        htitle.String = "COMPUTE OVERLAP";

        for a = 1:nDirGrid
            aidx = dirGridOrder(a);
            Fpen = wwIter(:, :, aidx);
            zsweep = sqrt(Fpen);
            zsweep = zsweep / zmax;
            zsweep(zsweep>1) = 1;
            zsweep = zsweep .* reshape(colMdl, [1, 1, 3]);
            % zsweep = zsweep./cl(2); % map intensity with same upper limit as ztrace

            % z = gather(ztrace + zsweep);
            z = ztrace + zsweep;
            himg.CData = z;

            hscores.XData(a) = xp(i) + u(aidx);
            hscores.YData(a) = yp(i) + v(aidx);
            if saveVideo
                frame = getframe(ax);
                for f = 1:nScanFrames
                    vwriter.writeVideo(frame);
                end
            else
                pause(1/frameRate/animSpeed);
            end
        end

        himg.CData = ztrace;

        % plot optimal angle
        x = xp(i) + cos(optimAngle) * polrad;
        y = yp(i) + sin(optimAngle) * polrad;
        hscoremx.XData = x;
        hscoremx.YData = y;
        plot([xp(i), x], [yp(i), y], 'color', colMdl, 'lineWidth', 2);
        set(hscoremx, "visible", "on");
        htitle.String = "FIND MINIMUM";
        if saveVideo
            frame = getframe(ax);
            for n = 1:nPauseFrames
                vwriter.writeVideo(frame);
            end
        else
            pause(1);
        end

        % update penalty trace with new sweep
        hscores.XData(:) = nan;
        hscores.YData(:) = nan;
        set(hscoremx, "visible", "off");
    end

    % Select sweep dir and add to penalty
    if selfDriving
        % using self-driven mode
        sweepDir = optimAngle;
        str = "ADD CHOSEN SWEEP";
    else
        % externally driven
        sweepDir = P.inputSweepDirs(i);
        x = xp(i) + cos(sweepDir) * polrad;
        y = yp(i) + sin(sweepDir) * polrad;
        hinputDir.XData = x;
        hinputDir.YData = y;
        plot([xp(i), x], [yp(i), y], 'color', colId, 'lineWidth', 2);
        str = "ADD IMPOSED SWEEP";
    end
    Fpen = mdl.calcTotalSweepProfile(xp(i), yp(i), sweepDir, false);
    Fpen = Fpen / mdl.profileMaxWeight;
    Fpen = reshape(Fpen, Fsz);
    F = F + Fpen;
    ztrace = sqrt(F);
    ztrace = ztrace / zmax;
    ztrace(ztrace>1) = 1;
    ztrace = mapcolors(ztrace, cl, bone());
    himg.CData = ztrace;
    htitle.String = str;
    if saveVideo
        frame = getframe(ax);
        for n = 1:nPauseFrames
            vwriter.writeVideo(frame);
        end
    else
        pause(1);
    end


end

if saveVideo
    vwriter.close();
end

end