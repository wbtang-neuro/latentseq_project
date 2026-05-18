function prob = decodePv(p)
%DECODEPV General population-vector correlation decoder
arguments
    p.spikeCounts = []; % cols=units
    p.tuning = []; % cols=units
    p.chunked = false;
    p.chunkSize = 5e4;
    p.useGpu = false;
    p.smoothing = 0;
    p.smoothSpikes = [];
    p.vt = "all";
end

if strcmpi(p.vt, "all")
   Y = p.spikeCounts;
else
    Y = p.spikeCounts(p.vt, :);
end

if ~isempty(p.smoothSpikes)
    Y = gsmooth(Y, p.smoothSpikes);
end

tuning = p.tuning./mean(p.tuning+eps, 'omitnan');

if p.useGpu
    tuning = gpuArray(tuning);
    Y = gpuArray(Y);
else
    tuning = gather(tuning);
    Y = gather(Y);
end

nd = size(p.tuning, 1);
nt = size(Y, 1);

if~p.chunked
    prob = corr(Y', tuning');
    if p.smoothing>0
        prob = gsmooth(prob, p.smoothing);
    end
else
    % Run batch decoding
    c = 1;
    nchunks = 0;
    prob = single(zeros(nt, nd));
    CHUNK_SIZE = p.chunkSize;

    while c<nt
        nchunks = nchunks + 1;
        inds = c + (1:CHUNK_SIZE) - 1;
        if inds(end) > nt
            inds = c:nt;
        end
        fprintf("Chunk %u, bins [%u - %u]\n", nchunks, inds([1, end]));
        YChunk = Y(inds, :);
        probChunk = corr(tuning', YChunk');
        probChunk = probChunk';
    
        if p.smoothing>0
            probChunk = gsmooth(probChunk, p.smoothing);
        end
        
        prob(inds, :) = gather(probChunk);
        prob(inds, :) = (probChunk);
        c = c + CHUNK_SIZE;
    end
end

end



