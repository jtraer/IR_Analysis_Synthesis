function PthStm=GtPthStm(Pth)
%* == GtPthStm.m i.e. Get Path Stem ==
%* Extracts all the parts of a path up until the last slash.  This is useful because the names returned by the dir function only nclude the parts after this
PthStm=[];
if ~isempty(strcmp(Pth,'/'));
    ndx=find(Pth=='/');
    ndx=max(ndx)-1;
    PthStm=Pth(1:ndx);
end
