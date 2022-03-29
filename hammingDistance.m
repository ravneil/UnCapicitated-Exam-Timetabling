% Hamming distance
function hd = hammingDistance(x, best)
hd=0;
d=length(x);
for i=1:d
 if x(i)~=best(i)
 hd=hd+1;
 end
end
