function [res] = plotAvgSweeps(D, dec)
%TRIGSWEEPS2D Summary of this function goes here
%   Detailed explanation goes here
cap = true;
if nargin==2
    cap=1;
end

sweeps = dec.sweeps;
chk = dec.chk;
% chk = D.thetaChunks;
% chk.speed = interp1(D.t, D.speed, D.thetaChunks.tStart, 'linear');
% chk.tCen = movmean(chk.tStart, [0,1]);

if isfield(dec, "poshpf")
    phpf = dec.poshpf;
else
    ptrue = [D.x, D.y];
    decpos = dec.decpos;
    phpf = decpos-ptrue;
    phpf = gsmooth(phpf, 1);
end
%% for right and left
% Trigger
vswp = chk.speed>.2 & [sweeps.nvalid]'>3 & [sweeps.nvalid]'<10 & [sweeps.straight]'>.5;
if isfield(dec, "vchk")
vswp = chk.speed>.2 & [sweeps.nvalid]'>3 & [sweeps.nvalid]'<20 & [sweeps.straight]'>.5 & dec.vchk;
end
% vswp = chk.speed>.15 & [sweeps.nvalid]'>2 & [sweeps.nvalid]'<10 & [sweeps.straight]'>.4;
% vswp(:)=true;
vswp(1:5) =false; vswp(end-5:end)=false;
sd = [sweeps.hpfPosDirection]';
% hd = sd; 
% hd(~vswp) = nan;
% hd(isnan(hd)) = interp1circ(D.t, D.hd, D.thetaChunks.tStart(isnan(hd)), 'linear');
% hd = gsmoothcirc(hd, 1);
% hd = interp1circ(D.t, D.hd, [sweeps.tStart]', 'linear');
hd = chk.hd;

egosd = circshift([circ_diff(sd);nan], 1);
egosd(~vswp)=nan;
egosd(find(abs(egosd) > pi/2)) = nan;

prevego = circshift(egosd, 1);
vswps.right = prevego>0 & vswp;
vswps.left = prevego<0 & vswp;

fds = fieldnamesstr(vswps);
npoints = 50;
for fd = fds
    trigIdx = find(vswps.(fd));
    sweepx = nan(numel(trigIdx),npoints);
    sweepy = nan(numel(trigIdx),npoints);
    for i = 1:numel(trigIdx)
        % pos = sweeps(trigIdx(i)).posHpf;
        % pos = sweeps(trigIdx(i)).pos;
        % pos = phpf(sweeps(trigIdx(i)).iSweep,:);
        pos = phpf(sweeps(trigIdx(i)).iStart:sweeps(trigIdx(i)).iStop,:);
        
        pos = rotate2d(pos, -hd(trigIdx(i))+pi/2);
        
%         pos(pos(:, 2)>1, :)=nan;
        npos = size(pos, 1);
        pos = interp1(1:npos, pos, linspace(1, npos, npoints), "cubic");
        sweepx(i, :) = pos(:, 1);
        sweepy(i, :) = pos(:, 2);
%         disp(i)
    end
    % Get the hpf positions of the next sweep
    if cap
        dist = sqrt(sweepx.^2 + sweepy.^2);
%         invalid = zscore(dist)>1;
        invalid = dist>1000;
        sweepx(invalid) = nan;
        sweepy(invalid) = nan;
    end
    avgy = mean(sweepy, 'omitnan');
    avgx = mean(sweepx, 'omitnan');

    % avgy = median(sweepy, 'omitnan');
    % avgx = median(sweepx, 'omitnan');
    
%     plot(nsamp)
%     yline(prctile(nsamp, 50))
%     avgy = mean(trigposy.*valid, 'omitnan');
%     avgx = mean(trigposx.*valid, 'omitnan');
%     tmp = avgy; tmp(nsamp<prctile(nsamp, 1))=nan;
%     [~, imx] = max(tmp);
%     [~, imn] = min(avgy(1:imx));
% %     imn = 2;
%     avgy = avgy(imn:imx);
%     avgx = avgx(imn:imx);
    res.(fd) = gsmooth([avgx', avgy'], 0);
%     scatter(trigposx(:, imx), trigposy(:, imx), 2, cols.(fd), 'markerFaceAlpha', .5)
%     inds = find(vswps.(fd));

end

end



