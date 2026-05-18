function sweeps = chunkThetaPosSweepsModules(D, dec, Y)
gridAxes = dec.gridAxes;
maxAngleDiff = 2;
maxDist = .15;%./spacingMeters;

%%
pTrue = [D.x, D.y];
p = gsmooth(dec.decpos, 0);
pSlow = dec.possm;
pHpf = wrapPosToGridTile(p-pSlow, gridAxes);
pHpf = gsmooth(pHpf, .8);
p = pHpf;
t = D.t;

% Compute travel dir and distance between points
pos1 = p(1:end-1, :); pos2 = p(2:end, :);
[dvec, dpos] = hexDistance(pos2,pos1,gridAxes);

travelDir = atan2(dvec(:, 2), dvec(:, 1));
travelDirdiff = circ_diff(travelDir);
travelDirdiff = [0;abs(travelDirdiff);0];
baddirdiff = travelDirdiff>maxAngleDiff;

nextdist = [dpos;nan];
prevdist = [nan;dpos];

% Scrap bins with too few active cells
p(gsmooth(dec.nactive, 1)<eps, :) = nan;
disp(size(p))
%%
Ysum = sum(Y, 2);
Yv = Y>0;

chk = D.thetaChunks;
ncyc = numel(chk.iStart);

%% 
for c = ncyc:-1:1 % avoid re-allocating by running backwards
    %

    i0 = chk.iStart(c);
    if c==ncyc
        i1 = numel(t);
    else
        i1 = chk.iStart(c+1);
    end
    ic = i0 : i1;
    pos = p(ic, :);
    
    Yc = Y(ic, :);
    Ysumc = Ysum(ic);
    Yvc = Yv(ic, :);
        
    [~,imx] = max(Ysumc);
    
    % Find adjacent samples less than 10cm away
    badnxt = find(nextdist(ic)>maxDist|baddirdiff(ic)|isnan(pos(:, 1)));
    badprv = find(prevdist(ic)>maxDist|baddirdiff(ic)|isnan(pos(:, 1)));

    istart = 1;
    iend = numel(ic);
    if ~isempty(badprv)
        istart = max([badprv(badprv<imx); 1]);
    end
    if ~isempty(badnxt)
        iend = min([badnxt(badnxt>=imx); iend]);
    end
    valid = istart:iend;
    invalid = true(size(ic));
    invalid(valid)= false;
    pos(invalid, :)=nan;
    

    %% maximizes distance from theta0
    hpfpos = pHpf(ic, :);
    hpfpos = pos;
    hpfpos(invalid, :) = nan;
    hpfdist = hypot(hpfpos(:, 1), hpfpos(:, 2));
    [~, imx0] = max(hpfdist);


    dpos = pos(1:imx0,:)-hpfpos(imx0,:);
    [~, istart] = max(hypot(dpos(:, 1), dpos(:, 1)));
%     [~, istart] = min(hpfdist);
    if istart>imx0
        istart = 1;
    end
%     dn = d;
%     dn(~hiAct) = nan;
%     [dmx, imx0] = max(dn); %update iend
    imx = imx0 + i0 - 1;
    iswp0 = i0 + istart - 1;
    isweep = iswp0 : imx;

    dpos = pos(imx0,:)-pos(istart,:);
    dmx = hypot(dpos(1), dpos(2));
    
    invalid = true(size(ic));
    invalid(istart:imx0)= false;
    pos(invalid, :)=nan;
    hpfpos(invalid, :)=nan;
    
    % define sweep axis
    ortend = rotate2d(hpfpos(imx0, :), pi/2); % Orthogonal axis
    ortline = [0,0,ortend];
    distsort = distancePointLine(hpfpos, ortline);
    totalvar = sum(var(hpfpos, 'omitnan'), 'omitnan');
    ortvar = var(distsort, 'omitnan');
    varratio = ortvar./totalvar;

    % Turuosity
    diffs = [0; vecnorm(diff(pos), 2,2)];
    pathlength = sum(diffs(istart:imx0), 'omitnan');

    s = struct();
    s.straight = varratio;
    s.pathlength = pathlength;
    s.netlength = dmx;
    s.turtuosity = pathlength./dmx;
    s.iStart = iswp0;
    s.iStop = imx;
    s.iSweep = isweep;
    s.iStop2 = iswp0+iend-1;
    s.tStart = t(s.iStart);
    s.tStop = t(s.iStop);
    s.posStart = pos(istart, :);
    s.posStop = pos(imx0, :);
    s.nvalid = numel(iswp0:imx);
    
    s.pos = p(isweep, :);
    s.posHpf = hpfpos(~invalid, :);
    s.posAll = p(ic, :);
    s.posHpfAll = pHpf(ic, :);

    s.posTrue = pTrue(s.iStart, :);
    s.length = dmx;

    s.activeCellsTotal = sum(Yvc, 1) > 0;       % logical vector indicating cells that fired >=1 spike in the cycle
    s.nActiveCells = uint16(sum(Yvc, 2));        % number of cells firing >=1 spike in current time bin
    s.nActiveCellsTotal = uint16(sum(s.activeCellsTotal)); % number of cells firing >=1 spike in the cycle
    s.nSpikes = uint16(sum(Yc, 2));             % summed spike count across cells, for each sample
    s.nSpikesTotal = uint16(sum(s.nSpikes));    % total summed spike count
    s.thetaCycleInds = [i0, i1];
    
    % Displacement vector of sweep path
    distVec = pos(imx0, :)-pos(istart, :);
    s.travelVector = distVec;
    s.travelDirection = atan2(distVec(2), distVec(1));
       
    % Displacement vector from true to final sweep position
    dpos = s.posHpf(end, :);
    s.hpfPosVector = dpos;
    s.hpfLength = hypot(dpos(1), dpos(2));
    s.hpfPosDirection = atan2(dpos(2), dpos(1));

    sweeps(c, 1) = s;

    
end



end

