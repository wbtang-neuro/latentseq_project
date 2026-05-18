function [spikeCounts, spikeInds1, spikeUnits1] = reconstructSpikeCounts(spikeInds, nt)
% Generate spike-count matrix from a spike time-bin indices of one or more
% units

if ~iscell(spikeInds), spikeInds = {spikeInds}; end
nu = numel(spikeInds);

spikeUnits = cellfun(@(a, b) {b*ones(size(a))}, spikeInds, num2cell(1:nu));
spikeInds1 = cat(1, spikeInds{:});
spikeUnits1 = cat(1, spikeUnits{:});

ind = sub2ind([nt nu], spikeInds1, spikeUnits1);
spikeCounts = accumarray(ind, single(1), [nt*nu, 1]);
% spikeCounts = accumarray(ind, 1, [nt*nu, 1], [], [], true);
spikeCounts = reshape(spikeCounts, nt, nu);

end