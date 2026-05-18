function Bfun2 = permuteBfun(Bfun)

Bfun2 = @(x,f) permute(reshape(Bfun(reshape(permute(x,[2,1,3]),size(x,2),[]),f),[],size(x,1),size(x,3)),[2,1,3]);

end