function h = annealHparam(harg, iter, learnRate)
% Anneal 

r = harg;
n = numel(harg);

if n == 2 
    % exponential rise or fall
    f = learnRate^(iter-1);
    if r(2) > r(1)
        % rise
        h = r(1) / f;
        h = min(h, r(2));
    else
        % fall
        h = r(1) * f;
        h = max(h, r(2));
    end
else
    
    if n == 3
        hvec = logspace(log10(r(1)), log10(r(2)), r(3));
    else
        % If hrange is a vector with more than 2 elements, we assume it
        % specifies a curve, with values indexed by the iteraction number.
        hvec = r;
    end
    
    if numel(hvec) >= iter
        % Index the vector with the iteration count
        h = hvec(iter);
    else
        % If the curve has too few elements, return the last element
        h = hvec(end);
    end
    
end

end