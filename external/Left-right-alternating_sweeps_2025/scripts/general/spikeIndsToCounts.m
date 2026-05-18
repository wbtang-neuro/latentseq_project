function [spikeCounts, spikeInds1, spikeUnits1] = spikeIndsToCounts(spikeInds, nt)
% Convert vectors of spike indices to vectors of spike counts

% INPUTS:
% spikeInds - a vector of spike time-bin indices, or a cell array
%             containing multiple such vectors.
% nt        - 

if ~iscell(spikeInds), spikeInds = {spikeInds}; end
nu = numel(spikeInds);

spikeUnits = cellfun(@(a, b) {b*ones(size(a))}, spikeInds, num2cell(1:nu));
spikeInds1 = cat(1, spikeInds{:});
spikeUnits1 = cat(1, spikeUnits{:});

ind = sub2ind([nt nu], spikeInds1, spikeUnits1);
spikeCounts = accumarray(ind, single(1), [nt*nu, 1]);
spikeCounts = reshape(spikeCounts, nt, nu);

end