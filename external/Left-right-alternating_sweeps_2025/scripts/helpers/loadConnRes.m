function [recs, allrecs] = loadConnRes(p)
% Loads conn results
S = SweepsSettings;
%% Load recordings
nrecs = numel(p.recNames);
clear recs
if p.loadHeavy
    for r = 1:nrecs
        fname = fullfile(p.basedir, "res_heavy_"+p.recNames(r)+".mat");
        if r==1
            tmp = load(fname);
            animal = strsplit(tmp.recName, "_");
            tmp.animal = animal(1);
            recs(r) = tmp;
        else
             tmp = load(fname);
             fds = fieldnamesstr(tmp);
             fdskeep = fieldnamesstr(recs);
             badfds = ~ismember(fds, fdskeep);
             tmp = rmfield(tmp, fds(badfds));
             animal = strsplit(tmp.recName, "_");
             tmp.animal = animal(1);
             recs(r) = tmp;
        end     
        disp(recs(r).recName)
    end
    allrecs = recs;
else
    for r = 1:nrecs
        fname = fullfile(p.basedir, "res_"+p.recNames(r)+".mat");
        if r==1
            tmp = load(fname);
            animal = strsplit(tmp.recName, "_");
            tmp.animal = animal(1);
            recs(r) = tmp;
        else
             tmp = load(fname);
             fds = fieldnamesstr(tmp);
             fdskeep = fieldnamesstr(recs);
             badfds = ~ismember(fds, fdskeep);
             tmp = rmfield(tmp, fds(badfds));
             animal = strsplit(tmp.recName, "_");
             tmp.animal = animal(1);
             recs(r) = tmp;
        end
        disp(recs(r).recName)
    end
    allrecs = recs;
end
%%
if p.loadTuning    
    [recs.tuning] = dealArr(0);
    for r = 1:nrecs
        fname = fullfile(p.basedir, "tuning_"+p.recNames(r)+".mat");
        recs(r).tuning = load(fname);
        disp(recs(r).recName)

    end
end

end