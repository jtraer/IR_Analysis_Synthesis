function ChckSm(x1,x2,str);
%* == ChckSm.m i.e. Check if two values are the same ==
%* outputs an error if they do not
if x1~=x2
    error(sprintf('%s are not the same.  No good will come of this.',str))
end
