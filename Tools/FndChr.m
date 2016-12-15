function ndx=FndChr(str1,x)
%* == FndChr.m i.e. Find a character within a string ==
%* Searches through a string for appearances of a character (e.g. slashes in a path or "."s in a structure )
ndx=[];
cnt=0;
for jj=1:length(str1);
    if strcmp(str1(jj),x); cnt=cnt+1;
        ndx(cnt)=jj;
    end
end
