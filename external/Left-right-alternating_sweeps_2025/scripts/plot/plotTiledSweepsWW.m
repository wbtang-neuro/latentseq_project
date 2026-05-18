function [outputArg1,outputArg2] = plotTiledSweepsWW(tl,decAll, sweeps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
decs = string(fields(decAll))';
ndecs = numel(decs);

%% For each dec, iterate through sweeps and move to next dec
for decname = decs
    for s = 1:numel(sweeps)
        dec = decAll.(decname);
        dec.ax = nexttile(tl);
%         dec.ax = nexttile(tl, 1);
%         if s == 5
%             sweeps(s).iStop=sweeps(s).iStop-1;
%         end
        nsteps = numel(sweeps(s).iStart:sweeps(s).iStop);

        dec = setUpHandles(dec, nsteps=nsteps, ax = dec.ax, istart = sweeps(s).iStart);
        axis image
        dec.ax.XLim = [-1.4,1.1]*.8; dec.ax.YLim = -fliplr(dec.ax.XLim)-0;
        dec.ax.XLim = [-1,1]*.8; dec.ax.YLim = (dec.ax.XLim)-.5-0;
        decAll.(decname) = dec;
        set(gca, 'YDir', 'normal', 'XDir', 'reverse')
        if s ==1
        ylabel(upper(decname), 'Color', 'k');%, 'HorizontalAlignment','center')
        end
    end
end
end

function dec = setUpHandles(dec, kwargs)
arguments
    dec
    kwargs.nsteps = 8;
    kwargs.istart = 0;
    kwargs.ax = gca;
    kwargs.normalizer = .04;
    kwargs.vidPath = [];
end
p = kwargs;
if isempty(p.istart)
    p.istart = p.nsteps;
end


prob = dec.prob;
t = dec.t;
% x = interp1(D.t, D.x, t, 'nearest');
% y = interp1(D.t, D.y, t, 'nearest');
% id = interp1(D.t, D.id, t, 'nearest');
x = dec.x*1; y = dec.y*1; id=dec.id; length = .4;%*0+.5; 
x = gsmooth(x, .4);y = gsmooth(y, .4);
id_x = [x, x+cos(id).*length];
id_y = [y, y+sin(id).*length];


% Set up video, if it exists
if isfield(dec, "vid")
   vid = dec.vid.vid;
   vidt = linspace(0, vid.Duration, vid.numFrames);
   vid.CurrentTime = vidt(p.istart+p.nsteps-0);
   frame = gather(readFrame(vid));
   hImg = imshow(frame+50, dec.vid.referenceFrame);
end
% Set up colors
frame = prob(:, :, 1);
nsteps = p.nsteps;
for s = 1:nsteps
    color =[];
    color(1,1,:) = mapcolors(s, [1,nsteps], 'cool');
    color = repmat(color, [size(frame), 1]);
    hdec(s) = imshow(color, 'XData',dec.grid.x, 'YData',dec.grid.y);
    hdec(s).AlphaData = gather(frame./max(frame(:)));
end
dec.handles.hImg = hImg;
dec.handles.hdec = hdec;

frameidx = p.istart+nsteps;
for stepidx = 1:nsteps
    frame = dec.prob(:, :, frameidx+stepidx-nsteps);
    frame = frame';
    dec.handles.hdec(stepidx).AlphaData = gather(frame./max(frame(:)));
end
hdir = plot(id_x(1, :), id_y(1, :), 'g','LineWidth',1);    
hdir.XData = id_x(frameidx, :); hdir.YData = id_y(frameidx, :);
plot(dec.x(p.istart:end), dec.y(p.istart:end), ':w')
end
