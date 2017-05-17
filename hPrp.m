function H=hPrp(H,C,Nbnds,flm,sb_fs,ftp);
%* == hPrp.m i.e. IR properties ==
%* Takes a structure of IR properties (i.e. H, output by hExtrct.m), filters into cochlear subbands (No. of bands specified by Nbnds and frequency limits by flm), and measures the decay properties of each subband. A structure of calibration IRs (i.e. C) can be used to remove speaker/microphone effects form the IR. 

set(0,'DefaultFigureVisible','off');
fntsz=15;

%* == Compute kurtosis in 10ms windows ==
%** Find the numbr of points in a 2ms window and make sure it is an even number
Nbn=ceil(0.002*H.fs); 
if Nbn>length(H.h)*2;
    Nbn=length(H.h)/4;
end
Nbn=Nbn+rem(Nbn,2); 
%** Repeat the first and last 5ms sections
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

%* == Compute the cochleagram ==
%** zeropad to avoid edge effects
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
eval(sprintf('! mkdir -p %s/Subbands_%d',H.Path,Nbnds));
eval(sprintf('! mkdir -p %s/Modes_%d',H.Path,length(H.Modes)));
BdBndsFlg=zeros(1,Nbnds);
for jbn=1:Nbnds; 
    fprintf('%s: Band %d/%d\n',H.Path,jbn,Nbnds);
    % Extract the subband
    tmp=Cgrm(jbn,:); 
    % rescale the ERs relative to the diffuse tail according to the face and volume speaker transfer functions
    tmp2=tmp; 
    % Compute spectral amplitude of the Sparse and Gaussian patches of the entire time series (for calibration only) -- shouldn't this be after we remove the noise floor?
    Nspc(jbn)=rms(gather(tmp2(Nndx)));
    Sspc(jbn)=rms(gather(tmp2(Sndx)));
    % take envelope
    tmp2=abs(hilbert([zeros(1,Npts) tmp2 zeros(1,Npts)]));  
    tmp2=tmp2(Npts+[1:Npts]);
    % resample
    tmp3=resample(tmp2,sb_fs,H.fs);  
    % Fit an exponential decay model
    [Pft,NsFlr,Test,FVE]=FtPlyDcy(tmp3,[1:length(tmp3)]/sb_fs,1,1);     
    % Do this for all the snapshots
    for jsnp=1:Nsnps
        snp=SnpCgrm(jbn,:,jsnp); 
        snp2=abs(hilbert([zeros(1,Npts) snp zeros(1,Npts)]));  
        snp2=snp2(Npts+[1:Npts]);
        snp3=resample(snp2,sb_fs,H.fs);  
        [sPft,snp_NsFlr,snp_Test,snp_FVE]=FtPlyDcy(snp3,[1:length(tmp3)]/sb_fs,1,1);     
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
    plot([1:length(tmp)]/H.fs,20*log10(abs(tmp))); 
    hold on
    plot([1:length(tmp2)]/H.fs,20*log10(abs(tmp2)),'c');
    plot([0:(length(tmp3)-1)]/sb_fs,20*log10(abs(tmp3)),'g:');
    plot([1:length(tmp)]/H.fs,Pft(2)+Pft(1)*[1:length(tmp)]/H.fs,'r--');
    plot([1:length(tmp)]/H.fs,(Pft(2)+sdDRR/2)+(Pft(1)-sdB/2)*[1:length(tmp)]/H.fs,'r:');
    plot([1:length(tmp)]/H.fs,(Pft(2)-sdDRR/2)+(Pft(1)+sdB/2)*[1:length(tmp)]/H.fs,'r:');
    plot([1 length(tmp)]/H.fs,Pft(2)+Pft(1)*Test*ones(1,2),'k--');
    if ~isempty(C)
        plot([1:length(tmp)]/H.fs,tmp3(1)-(60/C(2).RT60(jbn))*[1:length(tmp)]/H.fs,'m--');
        %** Check if recorded decay is less than or equal to the speaker-microphone IR
        if C(2).RT60(jbn)>(-60/Pft(1))*0.75;
            BdBndsFlg(jbn)=1;
            text(0.5*Test,Pft(2),1.001,sprintf('Danger: Speaker-Microphone IR RT60 is %d%% of recorded',round(100*C(2).RT60(jbn)/(-60/Pft(1)))));
        end
    end
    hold off
    title(sprintf('%s: Band %d',H.Path,jbn));
    set(gca,'xlim',[0 3*Test]);
    set(gca,'ylim',[Pft(2)+Pft(1)*Test-20 Pft(2)+20]);
    xlabel('Time (s)')
    ylabel('Power (dB)')
    saveas(gcf,sprintf('%s/Subbands_%d/%03d',H.Path,Nbnds,jbn),'jpg'); 
    close all
    % Compute the spectrum of the early reflections and diffuse section for only the sections where the IR is above the noise floor
    tmpERndx=Sndx; tmpERndx(find(tmpERndx)>Test*H.fs)=[];
    tmpGsndx=Nndx; tmpGsndx(find(tmpGsndx)>Test*H.fs)=[];
    spcER(jbn)=rms(gather(tmp2(tmpERndx)));
    spcGs(jbn)=rms(gather(tmp2(tmpGsndx)));
    % record subband values
    bbt(jbn)=-Pft(1); 	bt=-Pft(1);
    aa(jbn)=Pft(2);
    sdR(jbn)=sdRT60;
    sda(jbn)=sdDRR;
    raw_aa(jbn)=aa(jbn);
    if ~isempty(C)
        aa(jbn)=Pft(2)-(C(1).DRR(jbn)-mean(C(2).DRR)); 
    end
    alph=aa(jbn);
    Rtt(jbn)=60/bt;
    dttm(jbn)=(60+alph)/bt;
    NNs(jbn)=alph-bt*Test;
    TTest(jbn)=Test;
    %FVVE(jbn,:)=FVE;
    % compute a new subband with the noise floor removed
    infndx=ceil(Test*H.fs);
    if infndx>length(tmp); infndx=length(tmp)-1; end
    nCgrm(jbn,:)=tmp.*([ones(1,infndx) 10.^((-bt*[1:(length(tmp)-infndx)]/H.fs)/20)]);
    if ~isempty(C)
        nCgrm(jbn,:)=10^((-C(1).DRR(jbn)+mean(C(2).DRR))/20)*nCgrm(jbn,:);
    end
    % compute higher order models
    %for jp=1:Nprms; %fprintf('poly %d\n',jp)
    %    [Pft,NsFlr,Test,FVE]=FtPlyDcy(tmp3,tt3,jp,Qswtch);     
    %    FVVE(jbn,jp)=gather(FVE);
    %    PFFt(jp).Prms(jbn,:)=[Test Pft];
    %end
end

% resynthesize the new denoised IR estimate
nCgrm=[zeros(1,size(nCgrm,2)); nCgrm; zeros(1,size(nCgrm,2))];
nh=collapse_subbands([zeros(size(nCgrm)) nCgrm zeros(size(nCgrm))].',fltbnk);
nh=nh(Npts+[1:Npts]);
nCgrm=nCgrm(2:(end-1),:);
% rescale peak
Rscl=1/max(abs(nh));
nh=nh*Rscl;

%* Compute the spectrum of the attack
ndx=min(find(abs(nh)>prctile(abs(nh),90)));
cnt=0;
for jj=1:2:11; cnt=cnt+1;
    Nft=2^(jj+3);
    tmp=[nh; zeros(2*Nft,1)];
    Bgspc=zeros(Nft/2,1);
    for jstrt=1:(Nft/4);
        spc=fft(tmp(ndx+jstrt-1+[0:(Nft-1)]),Nft); 
        Bgspc=Bgspc+abs(spc(1:Nft/2))/(Nft/4); 
    end
    Attck(cnt).Spc=Bgspc(1:Nft/2);
    Attck(cnt).SpcIntrp=interp1([1:Nft/2]*H.fs/Nft,Bgspc,ff,'spline');
    Attck(cnt).ff=[1:Nft/2]*H.fs/Nft;
    Attck(cnt).T=Nft/H.fs;
end

% and compute spectrograms to find modes
[NsSgrm,Nsff,Nstt]=spectrogram(nh,32,16,32,H.fs);
if ~isempty(C)
    ClSgrm=mean(abs(C(1).Ns.Sgrm),2);
    ClSgrm=medfilt2(ClSgrm,[2 1],'Symmetric');
    NsSgrm=NsSgrm./(mean(abs(C(1).Ns.Sgrm),2)*ones(1,length(Nstt)));
end
Nft=512; Nbn=Nft;
if Nbn>length(nh)/8;
    Nbn=ceil(length(nh)/10);
    Nbn=Nbn+rem(Nbn,2);
end
[MdSgrm,Mdff,Mdtt]=spectrogram(nh,Nbn,Nbn/2,Nft,H.fs);
if ~isempty(C)
    ClSgrm=mean(abs(C(1).Md.Sgrm),2);
    ClSgrm=medfilt2(ClSgrm,[15 1],'Symmetric');
    MdSgrm=MdSgrm./(ClSgrm*ones(1,length(Mdtt)));
end
% compute modes
for jm=1:length(H.Modes);
    [~,ndx]=min(abs(Mdff-H.Modes(jm).cf));
    md=abs(MdSgrm(ndx,:));
    p=sum(md);
    %** find mode width
    cnt=0;
    FLG=0;
    while FLG==0; cnt=cnt+1;
        p2=sum(abs(MdSgrm(ndx+cnt,:)));
        if p2<p/2;
            FLG=1;
        end
        if cnt>1;
            if p2>=plst;
                FLG=1;
            end
        end
        if ndx+cnt==size(MdSgrm,1);
            FLG=1;
        end
        plst=p2;
    end
    undx=ndx+cnt-1;
    % and the lower frequencies
    if ndx>1;
        cnt=0;
        FLG=0;
        while FLG==0; cnt=cnt+1;
            p2=sum(abs(MdSgrm(ndx-cnt,:)));
            if p2<p/2;
                FLG=1;
            end
            if cnt>1;
                if p2>=plst;
                    FLG=1;
                end
            end
            if ndx-cnt==1;
                FLG=1;
            end
            plst=p2;
        end
        lndx=ndx-cnt+1;
    else
        lndx=ndx;
    end
    H.Modes(jm).Wd=max([(Mdff(undx)-Mdff(lndx)) (Mdff(2)-Mdff(1))]);
    % compute decay properties
    [Pft,NsFlr,Test,FVE]=FtPlyDcy(md,Mdtt,1,1);

    figure;
    plot(20*log10(abs(md))); axis tight;
    saveas(gcf,sprintf('%s/Modes_%d/%03d',H.Path,length(H.Modes),jm),'jpg'); 
    H.Modes(jm).OnPwr=Pft(2); 
    H.Modes(jm).RT60=60/abs(Pft(1));
    H.Modes(jm).MnPwr=mean(20*log10(abs(MdSgrm(ndx,:))));
end

%% save basic data to structure
H.nh=gather(nh);
H.krt=krt;
% and the channel values
H.spcER=spcER/mean([spcGs]);
H.spcGs=spcGs/mean([spcGs]);
H.spcAllGs=Nspc/mean(Nspc);
H.DRR=aa; %+20*log10(Rscl);
H.DRR_std=sda; 
H.RT60=Rtt;
H.RT60_std=sdR; 
H.NsFlr=NNs;
H.TTest=TTest;
H.BdBndsFlg=find(BdBndsFlg);
%H.FVE=FVVE;
%H.PFFt=PFFt;

% Spectrograms
%H.Ns=Ns;
H.Ns.Sgrm=NsSgrm;
H.Ns.ff=Nsff;
H.Ns.tt=Nstt;
%H.Md=Md;
H.Md.Sgrm=MdSgrm;
H.Md.ff=Mdff;
H.Md.tt=Mdtt;
H.Attck=Attck;

%* == Plot ==
save(sprintf('%s/H',H.Path),'H');
close all
fcnt=0;

%** => plot kurtosis
fcnt=fcnt+1; figure(fcnt); 
fprintf('%s: Plotting Kurtosis\n');
PltIRKrt(H,C,VrKrt,ftp);
saveas(gcf,sprintf('%s/Kurtosis',H.Path),'jpg');
saveas(gcf,sprintf('%s/Kurtosis',H.Path),ftp);

%** Plot Cochleagram
fcnt=fcnt+1; figure(fcnt)
subplot(2,1,1);
PltIRCgrm(H.h,H);
subplot(2,1,2);
PltIRCgrm(H.nh,H);
colormap(othercolor('Blues9',64));
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Cgram',H.Path),ftp);
saveas(gcf,sprintf('%s/Cgram',H.Path),'jpg');

%** => plot IR
fcnt=fcnt+1; figure(fcnt)
PltIRts(H);
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Raw_IR',Pth),'fig');
saveas(gcf,sprintf('%s/IR',H.Path),'jpg');
saveas(gcf,sprintf('%s/IR',H.Path),ftp);

%** => plot IR phase
fcnt=fcnt+1; figure(fcnt)
PltIRDyn(H);
%set(gca,'fontsize',fntsz);
%saveas(gcf,sprintf('%s/Raw_IR',Pth),'fig');
saveas(gcf,sprintf('%s/Dyn',H.Path),'jpg');
saveas(gcf,sprintf('%s/Dyn',H.Path),ftp);

%** plot subband properties
%*** Rtt
fcnt=fcnt+1; figure(fcnt)
PltIRRT60(H,C);
%set(gca,'fontsize',fntsz);
saveas(gcf,sprintf('%s/RT60',H.Path),'jpg');
saveas(gcf,sprintf('%s/RT60',H.Path),ftp);
%*** DRR
fcnt=fcnt+1; figure(fcnt)
PltIRDRR(H,C,raw_aa);
saveas(gcf,sprintf('%s/DRR',H.Path),'jpg');
saveas(gcf,sprintf('%s/DRR',H.Path),ftp);

%** plot mode properties
%*** Rtt
fcnt=fcnt+1; figure(fcnt)
PltIRMds(H,C);
%set(gca,'fontsize',fntsz);
saveas(gcf,sprintf('%s/Modes',H.Path),'jpg');
saveas(gcf,sprintf('%s/Modes',H.Path),ftp);

%** => Figure for Reverb paper
fcnt=fcnt+1; figure(fcnt)
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
fcnt=fcnt+1; figure(fcnt)
%PltIRSpc(H,C,NsSgrm,Nsff,MdSgrm,Mdff);
%set(gca,'fontsize',fntsz);
saveas(gcf,sprintf('%s/Spc',H.Path),'jpg');
saveas(gcf,sprintf('%s/Spc',H.Path),ftp);

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

