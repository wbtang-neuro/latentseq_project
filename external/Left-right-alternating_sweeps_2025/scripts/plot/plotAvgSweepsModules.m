function [res] = plotAvgSweepsModules(D, dec)
%TRIGSWEEPS2D Summary of this function goes here
%   Detailed explanation goes here
cap = true;
if nargin==2
    cap=1;
end

sweeps = dec.sweeps;
chk = D.chk;
chk.speed = interp1(D.t, D.speed, D.thetaChunks.tStart, 'linear');
chk.tCen = movmean(chk.tStart, [0,1]);

%% for right and left
% Trigger
vswp = chk.speed>.2 & [sweeps.nvalid]'>3 & [sweeps.nvalid]'<10 & [sweeps.straight]'>.5;
if isfield(dec, "vchk")
    vswp = chk.speed>.2 & [sweeps.nvalid]'>3 & [sweeps.nvalid]'<20 & [sweeps.straight]'>.5 & dec.vchk;
end
vswp(1:5) =false; vswp(end-5:end)=false;
sd = [sweeps.hpfPosDirection]';
hd = chk.hd; 
hd(~vswp) = nan;

egosd = circ_dist(chk.id, hd);
egosd = circ_dist(egosd, circ_mean(egosd(vswp)));

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
        pos = sweeps(trigIdx(i)).posHpf;
        
        pos = rotate2d(pos, -hd(trigIdx(i))+pi/2);
        
        npos = size(pos, 1);
        pos = interp1(1:npos, pos, linspace(1, npos, npoints), "cubic");
        sweepx(i, :) = pos(:, 1);
        sweepy(i, :) = pos(:, 2);
    end

    avgy = median(sweepy, 'omitnan');
    avgx = median(sweepx, 'omitnan');
    

    res.(fd) = gsmooth([avgx', avgy'], 0);

end

end

