function [recs, rez] = filterConnRes(allrecs, p)

%% Filter out invalid pairs 
% for each rec, discard all pairs that do not fulfill crieria
% add opt
recs = allrecs;
nrecs = numel(recs);
[recs.npairs] = dealArr(0);
[recs.nconn] = dealArr(0);
[recs.isconn] = dealArr(0);
[recs.isvalid] = dealArr(0);
[recs.globalConnProb] = dealArr(0);


totalpairs = 0;
totalconns = 0;
totalpairs_inter = 0;
totalconns_inter = 0;

for r = 1:nrecs
    res = recs(r);
    units = res.units;
    isvalid = [units.nspkFull]'>p.min_spikes & [units.nspkFull]>p.min_spikes;
    isvalid = isvalid;% & res.totalCnt>p.min_count;
    if p.only_good
        bothgood = string({units.ks2Label})'=="good" & string({units.ks2Label})=="good";
        isvalid = isvalid & bothgood;
    end
    if p.only_functional % kick all unclassified and non-grid cells
            ctypes = [units.cellType];
            isfunctional = ~ismember(ctypes, ["unclassified", "nongrid"]);
            notsame = ctypes'~=ctypes;
            isvalid = isfunctional & isfunctional' & isvalid;% & notsame;
    end
    isvalid = isvalid(:);
    isvalidnotsame = isvalid & notsame(:);

    % While not correct for n pairs, it works for reducing the dataset
   isconn = res.nspk_i>=p.min_spikes & res.nspk_j>=p.min_spikes;
   isconn = isconn & res.totalCnt>=p.min_count; %Should also be applied to all pairs

   [~, i] = (ismember(res.unitIds(:, 1), [units.id]));
   [~, j] = (ismember(res.unitIds(:, 2), [units.id]));
   
   nu = numel(units);
   inds = sub2ind([nu, nu], i, j);

    if p.only_good
        isconn = isconn & string({units(i).ks2Label})'=="good" & string({units(j).ks2Label})'=="good";
    end
    if p.only_functional
        ctypes_i = [units(i).cellType]';
        isfunctional_i = ~ismember(ctypes_i, ["unclassified", "nongrid"]);
        ctypes_j = [units(j).cellType]';
        isfunctional_j = ~ismember(ctypes_j, ["unclassified", "nongrid"]);
        isconn = isconn & isfunctional_i & isfunctional_j;% & ctypes_i~=ctypes_j;
        notsame = ctypes_j~=ctypes_i;
    end
    
    isconn = isconn & res.pkvalZ>p.min_peakZ & res.conn & res.pkWidth<3 & res.pkExcess>p.min_peakExcess & res.tpeak>.25e-3 & res.tpeak<5.1e-3;
%     isconn = isconn & res.pkvalZ>p.min_peakZ & res.conn & res.pkWidth<4 & res.pkExcess>p.min_peakExcess & res.tpeak>.25e-3 & res.tpeak<5.5e-3;
    % delete all nonconnected pairs
%     isconn = isconn & res.maxacausal<.2*(res.pkExcess) & res.acausal==0;
    isconn = isconn & res.maxacausal<2.5*(res.std) & res.acausal==0 & res.pkPval<.001;
    fds = string(fields(res))';
    for fd = fds
        if size(res.(fd), 1)==numel(isconn)
            res.(fd) = res.(fd)(isconn, :);
        end
    end
    isconnnotsame = isconn(:) & notsame(:);

    % Replace cellype labels for grid cells
    [units.cellTypeAll] = dealArr([units.cellType]);
    ctypesAll = [units.cellTypeAll];
    isbursty = ismember(ctypesAll, ["bursty", "prospective"]);
    isnonbursty = ismember(ctypesAll, ["nonbursty"]);
    [units(isbursty).cellType] = dealArr("gridBursty");
    [units(isnonbursty).cellType] = dealArr("gridNonbursty");
    res.units = units;
    res.nconn = sum(isconn);
    res.npairs = sum(isvalid);
    res.globalConnProb = res.nconn./res.npairs;
    res.isvalid = isvalid;
    res.isconn = isconn;
    totalpairs = totalpairs+sum(isvalid);
    totalconns = totalconns+sum(isconn);

    totalpairs_inter = totalpairs_inter+sum(isvalidnotsame);
    totalconns_inter = totalconns_inter+sum(isconnnotsame);
    recs(r)=res;
%     disp(r)
end
% recs(1).nconn;
%% Report total pairs and connections
str = struct();
str.pairs = sprintf("Total pairs: %u", totalpairs);
str.conns = sprintf("Total conns: %u", totalconns);
str.rate = sprintf("Connectivity: %.2f%", 100*totalconns./totalpairs);
str.nrats = sprintf("N rats: %u", numel(recs));
str.mean = sprintf("Conn mean: %.2f%", 100*mean([recs.globalConnProb]'));
str.std = sprintf("Conn std: %.2f%", 100*std([recs.globalConnProb]'));

str.interclass = sprintf("\nInterclass conns: %u/%u, %.2f%s\n", totalconns_inter, totalpairs_inter, 100*totalconns_inter/totalpairs_inter, "%");

fullstr = "";
for fd = fieldnamesstr(str)
    fullstr = fullstr+newline+str.(fd);
end
disp(fullstr)

rez.totalpairs = totalpairs;
rez.totalconns = totalconns;
rez.connrate = 100*totalconns./totalpairs;
rez.interconnrate = 100*totalconns_inter./totalpairs_inter;
end

