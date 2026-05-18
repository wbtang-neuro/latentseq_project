function n = calcNParams(prmInds)
fields = fieldnames(prmInds);
n = 0;
for f = 1:numel(fields)
    val = prmInds.(fields{f});
    if isstruct(val)
        nC = calcNParams(val);
    else
        nC = max(val);
    end
    n = max(n, nC);
end
end