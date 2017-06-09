function [D,V]=IntDrctTrnsFn(C,Drct_th,Drct_ph);
% Integrate Directional Transfer Function
% -- input: a structure containing multiple recordings made at different broadcast directions (C) and the polar (Drct_th) and azimuthal (Drct_ph) angles of th direct path between teh speaker and the microphone
% -- output: a structure with the properties of the direct signal (D) and volume integrated signal (V)

% Normalize to make all recorded spectra have the same length
for jc=1:length(C);
    Nft(jc)=length(C(jc).Spcff);
end
[Nft,ndx]=min(Nft);
for jc=1:length(C);
    C(jc).spc=interp1(abs(C(jc).spc),C(jc).Spcff,C(ndx).Spcff);
    C(jc).Spcff=C(ndx).Spcff;
end

%** interpolate spatial power spectrum over broadcast angles
for jc=1:length(C);
    th(jc)=str2num(C(jc).Meta.App.PolarAngle_fromTop);
    ph(jc)=str2num(C(jc).Meta.App.AzimuthalAngle_fromFront);
%===> DRR
%===> DRR_std
%===> RT60
%===> RT60_std
%===> spc
%===> NsFlr
    DRR(jc,:)=C(jc).DRR;
    DRR_std(jc,:)=C(jc).DRR_std;
    RT60(jc,:)=C(jc).RT60;
    RT60_std(jc,:)=C(jc).RT60_std;
    spc(jc,:)=C(jc).spc;
    NsFlr(jc,:)=C(jc).NsFlr;
end
%*** find direct signals
ndx1=find(th==Drct_th);
ndx2=find(ph==Drct_ph);
Drct_ndx=intersect(ndx1,ndx2);
if length(C)==1; Drct_ndx=1; end
D=C(Drct_ndx);
for jj=1:length(D)
    D(jj).Name='Direct';
end

%** sort polar and azimuthal angles into montonically increasing vectors and make a mesh 
[th2,ndx_th]=sort(unique(th));
[ph2,ndx_ph]=sort(unique(ph));
th2=[0:10:180];
ph2=[0:10:350];
[thh,phh]=meshgrid(th2,ph2);

%*** scroll through the polar angles to compute solid angle of each grid point
dth=180/(length(th2)+1);
dph=360/(length(ph2)+1);
for jth=1:size(thh,2);
    %**** => for each polar angle compute the band over which this value applies - check the edges are 0 and 180
    th_min=thh(1,jth)-dth/2;
    th_max=thh(1,jth)+dth/2;
    if th_min<0; 
        th_min=th_min+dth/2; 
    end
    if th_max>pi; 
        th_max=th_max-dth/2; 
    end
    %**** => for each band avergae over interpolated gridpoints 
    tmp=sum(sin(pi/180*linspace(th_min,th_max,1e3)))/1e3;
    %**** => store a matrix of solid angles
    dS(:,jth)=tmp*ones(size(thh,1),1);
end
SldAngl=dth*dph*dS*(pi/180)^2;

%** scroll through grid and assign values for each grid point.  Do this for 
%===> DRR
%===> DRR_std
%===> RT60
%===> RT60_std
%===> spc
%===> NsFlr
for jgrd=1:length(thh(:));
    th_j=thh(jgrd);
    ph_j=phh(jgrd);
    %*** => compute distance between this gridpoint and all the data
    th_df=abs(th-th_j);
    ph_df=abs(ph-ph_j);
    ph_df(find(ph_df)>180)=180-ph_df(find(ph_df)>180);
    %*** if this aligns with a measurement use that
    sm_ndx=find(th_df+ph_df==0);
    sm_ndx2=find(th_df==0);
    sm_ndx=unique([sm_ndx sm_ndx2]);
    if ~isempty(sm_ndx)
        df_ndx=sm_ndx;
        df=zeros(size(df_ndx));
    end
    %*** if not average the closest 3 measurements
    if isempty(sm_ndx);
        %**** find the closest three measurements
        df=sqrt((th_df).^2+(ph_df*sin(pi/180*th_j)).^2);
        [df,df_ndx]=sort(df);
        Nct=min([5 length(df)]);
        ndx=find(df==df(Nct));
        Nct=max(ndx);
        df=df(1:Nct); 
        df_ndx=df_ndx(1:Nct);
    end
%===> DRR
%===> DRR_std
%===> RT60
%===> RT60_std
%===> spc
%===> NsFlr
    ndx=(180-df(:))*ones(1,length(D.ff))/sum(180-df);
    vDRR(:,jgrd)=sum(DRR(df_ndx,:).*ndx,1);
    vDRR_std(:,jgrd)=sum(DRR_std(df_ndx,:).*ndx,1);
    vRT60(:,jgrd)=sum(RT60(df_ndx,:).*ndx,1);
    vRT60_std(:,jgrd)=sum(RT60_std(df_ndx,:).*ndx,1);
    vNsFlr(:,jgrd)=sum(NsFlr(df_ndx,:).*ndx,1);
    
    % for spectrum we need a longer vector for more frequencies 
    ndx=(180-df(:))*ones(1,length(D.Spcff))/sum(180-df);
    vspc(:,jgrd)=sum(spc(df_ndx,:).*ndx,1);
end

%** scroll through frequencies and interpolate calibration measures for each frequency
V=C(Drct_ndx); 
Nf=length(D.ff);
V.DRR=sum(vDRR.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));
V.DRR_std=sum(vDRR_std.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));
V.RT60=sum(vRT60.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));
V.RT60_std=sum(vRT60_std.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));
V.NsFlr=sum(vNsFlr.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));

Nf=length(D.Spcff);
V.spc=sum(vspc.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));

if length(C)==1;
    V=D;
end
V.Name='Omnidirectional';
