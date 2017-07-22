function [V]=IntDrctTrnsFn(C);
% Integrate Directional Transfer Function
% -- input: a structure containing multiple recordings made at different broadcast directions (C) and the polar (Drct_th) and azimuthal (Drct_ph) angles of th direct path between teh speaker and the microphone
% -- output: a structure with the properties of the direct signal (D) and volume integrated signal (V)

%* Form a new IR time series as a weighted sum of the individually measured IRs
%** scroll through the Calibration IRs
for jc=1:length(C);
    %*** => Extract broadcast angles
    th(jc)=str2num(C(jc).Meta.App.PolarAngle_fromTop);
    ph(jc)=str2num(C(jc).Meta.App.AzimuthalAngle_fromFront);
end
%** Compute a mesh of azimuths and polar angles
[Th,Ph]=meshgrid(linspace(0,pi,1e3),linspace(0,2*pi*(1-2e3),2e3));
%** Scroll through each gridpoint and assign it to it's closest measurement
Ndx=zeros(size(Th));
for jg=1:prod(size(Th));
    tdf=mod(Th(jg)-th,pi);
    pdf=mod(Ph(jg)-ph,2*pi);
    df=sqrt(tdf.^2+pdf.^2); % only strictly true for small angles, but should find teh closest
    [~,ndx]=min(abs(df));
    Ndx(jg)=ndx;
end
%** sum the weighted IRs
%*** find the length of the shortest IR
Npts=1e24;
sNpts=1e24;
Nsnps=1e24;
for jc=1:length(C);
    Npts=min([Npts length(C(jc).h)]);
    sNpts=min([sNpts size(C(jc).h_snps,1)]);
    Nsnps=min([Nsnps size(C(jc).h_snps,2)]);
    
end
h=zeros(Npts,1);
h_snps=zeros(sNpts,Nsnps);
for jc=1:length(C);
    %*** find the weight
    w=length(find(Ndx==jc))/length(Ndx);
    th=C(jc).h;
    h=h+w*th(1:Npts);
    th=C(jc).h_snps;
    h_snps=h_snps+w*th(sNpts,Nsnps);
end
V=C(1);
V.h=h;
V.h_snps=h_snps;
%* Rename
ndx=regexp(V.Path,V.Name);
Pth1=V.Path(1:ndx-1);
Pth2=V.Path(ndx+length(V.Name):end);
V.Path=[Pth1 'Omni' Pth2];
V.Name='CAL-Omni';
unix(sprintf('mkdir -p %s',V.Path));
V.Meta.App.PolarAngle_fromTop='Omni';
V.Meta.App.AzimuthalAngle_fromFront='Omni';
%* Update Meta data
mflds=fields(V.Meta);
envFLG=0; objFLG=0;
for jf=1:length(mflds);
    if strcmp(mflds{jf},'Env');
        envFLG=1;
    else strcmp(mflds{jf},'Obj'); 
        objFLG=1;
    end
end
if envFLG==1;
    flds=fields(V.Meta.Env);
    for jf=1:length(flds);
        eval(sprintf('V.Meta.Env.%s=''CAL-Omni'';',flds{jf}));
    end
end
if objFLG==1;
    flds=fields(V.Meta.Obj);
    for jf=1:length(flds);
        eval(sprintf('V.Meta.Obj.%s=''CAL-Omni'';',flds{jf}));
    end
end
