function H=hPrp(H,C,Nbnds,flm,sb_fs,ftp);
%* == hPrp.m i.e. IR properties ==
%* Takes a structure of IR properties (i.e. H, output by hExtrct.m), filters into cochlear subbands (No. of bands specified by Nbnds and frequency limits by flm), and measures the decay properties of each subband. A structure of calibration IRs (i.e. C) can be used to remove speaker/microphone effects form the IR. 

set(0,'DefaultFigureVisible','off');
fntsz=15;

% first we check that downsampling won't leave us with too little data for a sensible fitter
nL=0;
if nL<1e3;
    nL=length(H.h)*ceil(H.fs/sb_fs);
    sb_fs=2*sb_fs;
end
if sb_fs>H.fs; sb_fs=H.fs; end

% find data for Direct (D) and Omnidirectional (V) Speaker IRs
if ~isempty(C);
    for jc=1:length(C);
        dndx_t(jc)=length(regexp(C(jc).Name,'Front'));
        vndx_t(jc)=length(regexp(C(jc).Name,'Omni'));
    end
    dndx=find(dndx_t~=0);
    D=C(dndx);
    D=D(find([D.Channel]==H.Channel));
    vndx=find(vndx_t~=0);
    V=C(vndx);
    V=V(find([V.Channel]==H.Channel));
    H.CalibrationFiles=V(1).Path;
else
    D=[]; V=[];
    H.CalibrationFiles=[];
end


%* === Compute kurtosis in windows ===
BnL=5; % [Length of window in ms]
%** Find the numbr of points in the bin and make sure it is an even number
Nbn=ceil(BnL/1e3*H.fs); 
if Nbn>length(H.h)*2;
    Nbn=length(H.h)/4;
end
Nbn=Nbn+rem(Nbn,2); 
%** Repeat the first and last sections
tmp=[H.h(1:Nbn/2); H.h; H.h(end-Nbn/2+1:end)]; 
%** Scroll through data points
krt=zeros(length(H.h),1); 
stndx=0; Nstp=1; 
while stndx<(length(tmp)/Nstp-Nbn); stndx=stndx+1; 
    sc=tmp((stndx-1)*Nstp+[1:Nbn]); 
    krt(stndx,:)=kurtosis(sc); 
end; 
H.krt=krt;
%** Compute the expected variance of kurtosis for samples of Gaussian noise  
VrKrt=24*Nbn*(Nbn-1)^2/((Nbn-3)*(Nbn-2)*(Nbn+3)*(Nbn+5)); 
%** Classify data points as "Sparse" or "Noise-like" 
Sndx=find(krt>3+2*VrKrt);  %Sparse
Nndx=find(krt<=3+2*VrKrt); %Noise-like
H.Tail_ndx=Nndx;
%** Compute the crossover to Gaussian statistics as the point at which there has been as many Gaussian points as sparse (this is a stable measure but it is also arbitrary and crude)
%*** Find the maximum (preumably this is near the first arrival)
[~,mxndx]=max(abs(H.h));
Sndx(find(Sndx<=mxndx))=[];
Nndx(find(Nndx<=mxndx))=[];
NGs=0; NER=1; cnt=mxndx; 
while (NGs<=NER&&cnt<length(krt)); cnt=cnt+1; 
    NGs=length(find(Nndx<=cnt)); 
    NER=length(find(Sndx<=cnt)); 
end; 
H.Tgs=cnt/H.fs;

%* === remove the direct and omnidirectional speaker transfer function from IR time series ===
H.h_before_removing_speaker_TF=H.h;
if ~isempty(V);
    h=H.h;
    vSpc=V.Attck(3).RwSpc;
    vff=V.Attck(3).ff;
    h_Tail_Calibrated=RmvSpkTrnsFn(h,H.fs,vSpc,vff);
    dSpc=D.Attck(3).RwSpc;
    dff=D.Attck(3).ff;
    h_Direct_Calibrated=RmvSpkTrnsFn(h,H.fs,dSpc,dff);
    % splice the two calibrated IR together with a crossfade
    CrssL=10; %crossfade length in ms
    if CrssL>2*H.Tgs;
        CrssL=H.Tgs;
    end
    Ncrss=ceil(CrssL/1e3*H.fs);
    Ncrss=Ncrss+rem(Ncrss,2);
    N1=cnt-Ncrss/2;
    N2=length(H.h)-N1-Ncrss;
    wn1=[ones(N1,1); linspace(1,0,Ncrss).'; zeros(N2,1)];
    wn2=[zeros(N1,1); linspace(0,1,Ncrss).'; ones(N2,1)];
    if length(wn1)>length(h); % as might be the case if Tgs is the length of the IR
        wn1=wn1(1:length(h));
        wn2=wn2(1:length(h));
    end
    H.h=wn1.*h_Direct_Calibrated+wn2.*h_Tail_Calibrated;
    % repeat this for all the snapshots
    Nsnps=size(H.h_snps,2);
    for js=1:Nsnps;
        h=H.h_snps(:,js);
        th_Tail_Calibrated=RmvSpkTrnsFn(h,H.fs,vSpc,vff);
        th_Direct_Calibrated=RmvSpkTrnsFn(h,H.fs,dSpc,dff);
        H.h_snps(:,js)=wn1.*th_Direct_Calibrated+wn2.*th_Tail_Calibrated;
    end
end

%* Remove the noise floor in subbands
%** == Compute the cochleagram ==
%*** zeropad to avoid edge effects
Npts=length(H.h);
[fltbnk,ff,erbff]=make_erb_cos_filters(3*Npts,H.fs,Nbnds,flm(1),flm(2));
Cgrm=generate_subbands([zeros(Npts,1); H.h; zeros(Npts,1)].',fltbnk);
Cgrm=Cgrm(Npts+[1:Npts],:).'; 
%** Remove the extreme bands
Cgrm=Cgrm([2:(end-1)],:);
H.off=ff;
H.ff=ff([2:(end-1)]);
%** Repeat this for the snapshots
Nsnps=size(H.h_snps,2);
SnpCgrm=zeros(size(Cgrm,1),size(Cgrm,2),Nsnps);
for jsnp=1:Nsnps;
    sCgrm=generate_subbands([zeros(Npts,1); H.h_snps(:,jsnp); zeros(Npts,1)].',fltbnk);
    sCgrm=sCgrm(Npts+[1:Npts],:).'; 
    sCgrm=sCgrm([2:(end-1)],:);
    SnpCgrm(:,:,jsnp)=sCgrm;
end

%* == Scroll through cochlear channels ==
unix(sprintf('! mkdir -p %s/Subbands_%d',H.Path,Nbnds));
BdBndsFlg=zeros(1,Nbnds);
for jbn=1:Nbnds; 
    fprintf('%s: Band %d/%d\n',H.Path,jbn,Nbnds);
    % Extract the subband
    tmp=Cgrm(jbn,:);
    %** record the subband peak amplitude
    H.SbAmp(jbn)=20*log10(max(abs(tmp)));
    % rescale the ERs relative to the diffuse tail according to the face and volume speaker transfer functions
    tmp2=tmp; 
    % take envelope
    tmp2=abs(hilbert([zeros(1,Npts) tmp2 zeros(1,Npts)]));  
    tmp2=tmp2(Npts+[1:Npts]);
    % resample
    tmp3=resample(tmp2,sb_fs,H.fs);  
    N2ndx=ceil(Nndx*sb_fs/H.fs);
    % Fit an exponential decay model
    tt=[1:length(tmp3)]/sb_fs; 
    [Pft,NsFlr,Test,FVE]=FtPlyDcy(tmp3(N2ndx),tt(N2ndx),1,1);
    alph=Pft(2); bt=-Pft(1);
    % Do this for all the snapshots
    for jsnp=1:Nsnps
        snp=SnpCgrm(jbn,:,jsnp); 
        snp2=abs(hilbert([zeros(1,Npts) snp zeros(1,Npts)]));  
        snp2=snp2(Npts+[1:Npts]);
        snp3=resample(snp2,sb_fs,H.fs);  
        [sPft,snp_NsFlr,snp_Test,snp_FVE]=FtPlyDcy(snp3(N2ndx),tt(N2ndx),1,1);     
        snpB(jsnp)=-sPft(1);
        snpRT60(jsnp)=60/-sPft(1);
        snpDRR(jsnp)=sPft(2);
    end
    %** Get variances
    sdB=std(snpB);
    sdRT60=std(snpRT60);
    sdDRR=std(snpDRR);
    % Plot
    figure; set(gcf,'visible','off')
    BdFLG=PltIRSbbnd(H,V,jbn,tmp,tmp2,tmp3,sb_fs,Pft,Test,sdDRR,sdB);
    saveas(gcf,sprintf('%s/Subbands_%d/%03d',H.Path,Nbnds,jbn),'jpg'); 
    close all
    % compute a new subband with the noise floor removed
    infndx=ceil(Test*H.fs);
    if infndx>length(tmp); infndx=length(tmp)-1; end
    ntmp=tmp.*([ones(1,infndx) 10.^((-bt*[1:(length(tmp)-infndx)]/H.fs)/20)]); 
    nCgrm(jbn,:)=ntmp;
    %** compute subband properties
    %*** spectrum of the early reflections and diffuse section 
    %*** (for only the sections where the IR is above the noise floor)
    tmpERndx=Sndx; tmpERndx(find(tmpERndx)>Test*H.fs)=[];
    tmpGsndx=Nndx; tmpGsndx(find(tmpGsndx)>Test*H.fs)=[];
    H.spcER(jbn)=rms(ntmp(tmpERndx));
    H.spcGs(jbn)=rms(ntmp(tmpGsndx));
    H.spcAllGs(jbn)=rms(ntmp(Nndx));
    %*** decay statistics
    H.DRR(jbn)=Pft(2); 
    H.DRR_std(jbn)=sdDRR; 
    H.RT60(jbn)=60/(-Pft(1));
    H.RT60_std(jbn)=sdRT60; 
    %*** diagnostic information
    H.NsFlr(jbn)=alph-bt*Test;
    H.TTest(jbn)=Test;
    H.BdBndsFlg(jbn)=BdBndsFlg(jbn);
end
H.BdBndsFlg=find(H.BdBndsFlg==1);
% resynthesize the new denoised IR estimate
nCgrm=[zeros(1,size(nCgrm,2)); nCgrm; zeros(1,size(nCgrm,2))];
nh=collapse_subbands([zeros(size(nCgrm)) nCgrm zeros(size(nCgrm))].',fltbnk);
nh=nh(Npts+[1:Npts]);
nCgrm=nCgrm(2:(end-1),:);
H.h_before_removing_noisefloor=H.h;
H.h=nh;

%* measure broadband properties
%** Spectrum
tmp=[zeros(size(H.h)); H.h; zeros(size(H.h))];
nft=2^ceil(log2(length(tmp)));
spc=fft(tmp,nft);
H.spc=spc(1:end/2);
H.Spcff=[1:nft/2]*H.fs/nft;
%* Compute the spectrum in time windows
ndx=min(find(abs(H.h)>prctile(abs(nh),90)));
cnt=0;
for jj=1:2:11; cnt=cnt+1;
    Nft=2^(jj+3);
    tmp=[nh; zeros(2*Nft,1)];
    Bgspc=zeros(Nft/2,1);
    for jstrt=1:(Nft/4);
        spc_t=tmp(ndx+jstrt-1+[0:(Nft-1)]);
        spc_t=spc_t.*hann(length(spc_t));
        spc=fft(spc_t,Nft); 
        Bgspc=Bgspc+abs(spc(1:Nft/2))/(Nft/4); 
    end
    Attck(cnt).RwSpc=Bgspc(1:Nft/2);
    Attck(cnt).SpcIntrp=interp1([1:Nft/2]*H.fs/Nft,Bgspc,ff,'spline');
    Attck(cnt).ff=[1:Nft/2]*H.fs/Nft;
    Attck(cnt).T=Nft/H.fs;
end
H.Attck=Attck;

%Search for modes
fprintf('searching for Modes...\n')
%H.Modes=hExtrctMds(H,2048);
H.Modes=hExtrctMds(H,1024); % for now this is faster
fprintf('%d modes found.\n',length(H.Modes))
unix(sprintf('! mkdir -p %s/Modes_%d',H.Path,Nbnds));
%** fit decay properties of modes
Npts=length(H.h);
for jm=1:length(H.Modes);
    Md=H.Modes(jm);
    [fltbnk,ff,erbff]=make_erb_cos_filters(3*Npts,H.fs,1,Md.cf-Md.bw/2,Md.cf+Md.bw/2);
    md=generate_subbands([zeros(Npts,1); H.h; zeros(Npts,1)].',fltbnk).';
    md=md(2,:);
    md=abs(hilbert(md));
    md=md(Npts+[1:Npts]);
    md=interp1([1:Npts]/H.fs,md,tt);
    % measure mode decay properties
    [Pft,NsFlr,Test,FVE]=FtPlyDcy(md,tt,1,1);
    H.Modes(jm).OnPwr=Pft(2);
    H.Modes(jm).RT60=60/abs(Pft(1));
    H.Modes(jm).MnPwr=mean(20*log10(md));

    figure;
    plot(tt,20*log10(abs(md))); axis tight;
    xlabel('Time (s)')
    ylabel('Power (dB)')
    title(sprintf('%2.2kHz, FV=%2.2f',Md.cf/1e3,Md.FV))
    saveas(gcf,sprintf('%s/Modes_%d/%03d',H.Path,Nbnds,jm),'jpg'); 
end

% and compute spectrograms to find modes
%[NsSgrm,Nsff,Nstt]=spectrogram(nh,32,16,32,H.fs);
%if ~isempty(D)
%    tD=D(find([D.Channel]==H.Channel));
%    ClSgrm=mean(abs(tD.Ns.Sgrm),2);
%    ClSgrm=medfilt2(ClSgrm,[2 1],'Symmetric');
%    %NsSgrm=NsSgrm./(mean(abs(tD.Ns.Sgrm),2)*ones(1,length(Nstt)));
%end
%Nft=512; Nbn=Nft;
%if Nbn>length(nh)/8;
%    Nbn=ceil(length(nh)/10);
%    Nbn=Nbn+rem(Nbn,2);
%end
%[MdSgrm,Mdff,Mdtt]=spectrogram(nh,Nbn,Nbn/2,Nft,H.fs);
%if ~isempty(D)
%    tD=D(find([D.Channel]==H.Channel));
%    ClSgrm=mean(abs(tD.Md.Sgrm),2);
%    ClSgrm=medfilt2(ClSgrm,[15 1],'Symmetric');
%    %MdSgrm=MdSgrm./(ClSgrm*ones(1,length(Mdtt)));
%end
%% compute modes
%for jm=1:length(H.Modes);
%    [~,ndx]=min(abs(Mdff-H.Modes(jm).cf));
%    md=abs(MdSgrm(ndx,:));
%    p=sum(md);
%    %** find mode width
%    cnt=0;
%    FLG=0;
%    while FLG==0; cnt=cnt+1;
%        p2=sum(abs(MdSgrm(ndx+cnt,:)));
%        if p2<p/2;
%            FLG=1;
%        end
%        if cnt>1;
%            if p2>=plst;
%                FLG=1;
%            end
%        end
%        if ndx+cnt==size(MdSgrm,1);
%            FLG=1;
%        end
%        plst=p2;
%    end
%    undx=ndx+cnt-1;
%    % and the lower frequencies
%    if ndx>1;
%        cnt=0;
%        FLG=0;
%        while FLG==0; cnt=cnt+1;
%            p2=sum(abs(MdSgrm(ndx-cnt,:)));
%            if p2<p/2;
%                FLG=1;
%            end
%            if cnt>1;
%                if p2>=plst;
%                    FLG=1;
%                end
%            end
%            if ndx-cnt==1;
%                FLG=1;
%            end
%            plst=p2;
%        end
%        lndx=ndx-cnt+1;
%    else
%        lndx=ndx;
%    end
%    H.Modes(jm).Wd=max([(Mdff(undx)-Mdff(lndx)) (Mdff(2)-Mdff(1))]);
%    % compute decay properties
%    [Pft,NsFlr,Test,FVE]=FtPlyDcy(md,Mdtt,1,1);

%* == Plot ==
close all
fcnt=0;

%** => plot IR
fcnt=fcnt+1; figure(fcnt)
PltIRts(H);
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Raw_IR',Pth),'fig');
saveas(gcf,sprintf('%s/IR',H.Path),ftp);

%** => plot kurtosis
fcnt=fcnt+1; figure(fcnt); 
fprintf('%s: Plotting Kurtosis\n');
PltIRKrt(H,C,VrKrt);
saveas(gcf,sprintf('%s/Kurtosis',H.Path),ftp);

%** => plot spectrum
fcnt=fcnt+1; figure(fcnt); 
fprintf('%s: Plotting Spectrum\n');
PltIRSpc(H,C);
saveas(gcf,sprintf('%s/Spc',H.Path),ftp);

%** Plot Cochleagram
fcnt=fcnt+1; figure(fcnt)
PltIRCgrm(H.h,H);
colormap(othercolor('Blues9',64));
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Cgram',H.Path),ftp);
saveas(gcf,sprintf('%s/Cgram',H.Path),'jpg');


%** plot subband properties
%*** Rtt
fcnt=fcnt+1; figure(fcnt)
PltIRRT60(H,C);
%set(gca,'fontsize',fntsz);
saveas(gcf,sprintf('%s/RT60',H.Path),ftp);
%*** DRR
fcnt=fcnt+1; figure(fcnt)
PltIRDRR(H,C);
saveas(gcf,sprintf('%s/DRR',H.Path),ftp);

%** plot mode properties
%*** Rtt
fcnt=fcnt+1; figure(fcnt)
PltIRMds(H,C);
%set(gca,'fontsize',fntsz);
saveas(gcf,sprintf('%s/Modes',H.Path),ftp);

%** => plot IR phase
fcnt=fcnt+1;figure(fcnt)
PltIRDyn(H);
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Raw_IR',Pth),'fig');
saveas(gcf,sprintf('%s/Dyn',H.Path),'jpg');
%saveas(gcf,sprintf('%s/Dyn',H.Path),ftp);

%======= Below here is a mess of old plots we don't use anymore ======
%** => Figure for Reverb paper
%fcnt=fcnt+1; figure(fcnt)
%subplot(1,6,1);
%PltIRts(H);
%subplot(1,6,2);
%PltIRPkAmp(H);
%subplot(1,6,3);
%PltIRKrt(H);
%subplot(1,6,4);
%PltIRSpc(H);
%subplot(1,6,5);
%PltIRDRR(H);
%subplot(1,6,6);
%PltIRRT60(H);
%saveas(gcf,sprintf('%s/Fg4Ppr',H.Path),'jpg');   
%saveas(gcf,sprintf('%s/Fg4Ppr',H.Path),ftp);   

%** Plot Spectrogram
%fcnt=fcnt+1; figure(fcnt)
%PltIRSpc(H,C,NsSgrm,Nsff,MdSgrm,Mdff);
    %for jm=1:length(H.Modes);
    %    [~,mndx]=min(abs(H.Spcff-H.Modes(jm).cf+H.Modes(jm).bw));
    %    [~,mxdx]=min(abs(H.Spcff-H.Modes(jm).cf-H.Modes(jm).bw));
    %    [~,f0ndx]=min(abs(H.Spcff-H.Modes(jm).cf));
    %    plot(20*log10(abs(H.spc(f0ndx)/min(abs(H.spc))))*sin(linspace(0,pi,(mxdx-mndx+1)))+20*log10(min(abs(H.spc))),H.Spcff(mndx:mxdx)/1e3,'r');
    %    plot(20*log10(abs(H.spc(f0ndx))),H.Spcff(f0ndx)/1e3,'ro')
    %end
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Spc',H.Path),'jpg');
%saveas(gcf,sprintf('%s/Spc',H.Path),ftp);

%subplot(2,1,1);
%pcolor(Nstt,Nsff/1e3,plt);
%axis xy; shading flat
%xlabel('Time (s)');
%ylabel('Frequency (kHz)');
%title([H.Path ': Noise Spectrogram']);
%set(gca,'clim',max(max(plt))+[-80 0]);
%set(gca,'yscale','log');
%colorbar
%colormap(othercolor('Blues9',64));
%saveas(gcf,sprintf('%s/NsSgram',H.Path),ftp);
%%** Plot Spectrogram
%fcnt=fcnt+1; figure(fcnt)
%subplot(2,1,1)
%plt=20*log10(abs(MdSgrm)); 
%if size(plt,2)==1;
%    plt=[plt plt];
%    Mdtt=[Mdtt 2*Mdtt];
%end
%pcolor(Mdtt,Mdff/1e3,plt);
%axis xy; shading flat
%xlabel('Time (s)');
%ylabel('Frequency (kHz)');
%title([H.Path ': Mode Spectrogram']);
%set(gca,'clim',max(max(plt))+[-80 0]);
%set(gca,'yscale','log');
%colorbar
%colormap(othercolor('Blues9',64));
%for jplt=1:5;
%    if jplt==1; 
%        for jf=1:size(plt,1);
%            [mx(jf),mxndx(jf)]=max(plt(jf,1:ceil(size(plt,2)/10)));
%        end
%        ndx=ceil(mean(mx.*mxndx)/mean(mx));
%    elseif jplt==5; ndx=ceil(size(plt,2)/2);
%    else ndx=ceil((jplt-1)*size(plt,2)/8);
%    end
%    ndx=max([ndx 1]);
%    if ndx<size(plt,2)-1;
%        subplot(2,1,1); hold on;
%        plot(Mdtt(ndx)*ones(1,2),Mdff([2 end])/1e3,'k--'); 
%        subplot(2,5,5+jplt);
%        plot(mean(plt(:,ndx+[0:1]),2),Mdff);
%        axis tight
%        set(gca,'xlim',max(max(plt))+[-80 0]);
%        set(gca,'yscale','log')
%    end
%end
%saveas(gcf,sprintf('%s/MdSgram',H.Path),'jpg');
%saveas(gcf,sprintf('%s/MdSgram',H.Path),ftp);

%%** Plot Spectra of attack
%%** => plot spectrum
%fcnt=fcnt+1; figure(fcnt)
%for jplt=1:length(H.Attck);
%    Hspc=H.Attck(jplt).Spc;
%    if ~isempty(C)
%        Cspc=C(1).Attck(jplt).Spc;
%        Hspc=Hspc./Cspc;
%    end
%    Hspc=20*log10(abs(Hspc));
%    plot(Hspc,H.Attck(jplt).ff);
%    hold on;
%    lgnd{jplt}=sprintf('%2.1fms',H.Attck(jplt).T*1e3);
%end
%hold off
%xlabel('Power (db)');
%ylabel('Frequency (kHz)')
%set(gca,'yscale','log')
%legend(lgnd); 
%title([H.Path ': IR Attack spectra'])
%saveas(gcf,sprintf('%s/IR_AttckSpc',H.Path),'jpg');
%saveas(gcf,sprintf('%s/IR_AttckSpc',H.Path),ftp);
%%** => plot spectrum interpolated to Cgrm resolution
%fcnt=fcnt+1; figure(fcnt)
%for jplt=1:length(H.Attck);
%    Hspc=H.Attck(jplt).SpcIntrp;
%    if ~isempty(C)
%        Cspc=C(1).Attck(jplt).SpcIntrp;
%        Hspc=Hspc./Cspc;
%    end
%    Hspc=20*log10(abs(Hspc));
%    plot(Hspc,ff);
%    hold on;
%    lgnd{jplt}=sprintf('%2.1fms',H.Attck(jplt).T*1e3);
%end
%hold off
%xlabel('Power (db)');
%ylabel('Frequency (kHz)')
%set(gca,'yscale','log')
%legend(lgnd); 
%title([H.Path ': IR Attack spectra'])
%saveas(gcf,sprintf('%s/IR_AttckSpcIntrp',H.Path),'jpg');
%saveas(gcf,sprintf('%s/IR_AttckSpcIntrp',H.Path),ftp);

