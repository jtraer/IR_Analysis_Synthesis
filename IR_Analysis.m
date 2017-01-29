%* == IR_Analysis.m == 
%** Analyzes the IRs extracted by IR_Extract.m
%** This makes use of the following functions:
%*** - hExtrct.m	: extracts the IR time series
%*** - hPrp.m 	: measures the properties of the IR time series
%*** - GtMtDt.m	: Reads metadata about recording from a textfile

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')

%* == Specify Inputs == 

%** = File with IRs extracted by IR_Extract.m =
Hpth='H_raw_RoomReverb_2IRs_19-Dec-2016';
Hpth='H_raw_Foam_Board_2IRs_24-Jan-2017';
Hpth='H_raw_Board_2IRs_25-Jan-2017';
Hpth='H_raw_Slt_Crt_1IRs_29-Jan-2017';
%Hpth='H_raw_FirstTest_Marble_1IRs_06-Dec-2016'
%** = File with Calibration IRs (Optional, but recommended) =
% These are re-recorded broadcasts (as in the measurements) but in as close to anechoic conditions as possible. These are required to ensure we are recording the properties of the space, not of the speaker/microphone/soundcard.
%Cpth='H_raw_Cal_1IRs_12-Dec-2016';
%Cpth='H_raw_Cal_Zipp-DR40_14IRs_19-Dec-2016';
Cpth='H_raw_Cal_Contact_2IRs_24-Jan-2017';
Cpth='H_raw_CAL_Crt_1IRs_29-Jan-2017';
%Cpth=[];
%** = Number of cochlear subbands for analysis =
Nbnds=[40];
%** = Frequency limits in Hz =
flm=[100 20e3];
%** = Frequency of subband envelopes in Hz =
Sb_fs=1e3;

%** = Specify MetaData we want to record (these can be added or removed arbitrarily) =
%mcnt=0;
%mcnt=mcnt+1;Mt{mcnt}='App.Mic';

%* == Calibrate Apparatus ==
%** search for calibration files 
if ~isempty(Cpth)
    load(Cpth); %H
    %*** => Check to see if the IR is already analyzed
    Pth=GtPthStm(GtPthStm(H(1).Path)); % A path to save the calibration files
    Pth=sprintf('%s/Cal_%03dBnds_%05d-%05dHz_ft%05dHz',Pth,Nbnds,flm(1),flm(2),Sb_fs);
    if exist(sprintf('%s/H.mat',Pth));
        load(sprintf('%s/H.mat',Pth));
    %*** => If not scroll through IRs
    else
        ccnt=0;
        for jc=1:length(H);
            tC=H(jc);
            %*** TODO add this so we don't repeat analyses unnecessarily
            %*** => Analyze
            tC=hPrp(tC,[],Nbnds,flm,Sb_fs);
            %*** => compile into a structure
            if jc==1;
                C=tC;
                %*** save some useful info
                Nf=length(tC.ff); % number of frequency bins
                OtPth=GtPthStm(GtPthStm(tC.Path)); % A path to save the calibration files
                OtPth=sprintf('%s/Cal_%03dBnds_%05d-%05dHz_ft%05dHz',OtPth,Nbnds,flm(1),flm(2),Sb_fs);
                unix(sprintf('mkdir -p %s',OtPth));
            else
                C=[C tC];
            end
        end % for jc=1:length(Dc)
        C_Left=C(find([C.Channel]==1));
        C=C_Left; % This is a HACK!!!!
    %	%** interpolate spatial power spectrum
        for jc=1:length(C);
            th(jc)=str2num(C(jc).Meta.App.PolarAngle_fromTop);
            ph(jc)=str2num(C(jc).Meta.App.AzimuthalAngle_fromFront);
            DRR(jc,:)=C(jc).DRR;
        end
        %*** find direct 
        ndx1=find(th==0);
        ndx2=find(ph==0);
        Drct_ndx=intersect(ndx1,ndx2);
        if length(C)==1; Drct_ndx=1; end
        D=C(Drct_ndx);
        %** sort these into montonically increasing vectors and make a mesh of polar and azimuthal angles
        [th2,ndx_th]=sort(unique(th));
        [ph2,ndx_ph]=sort(unique(ph));
        th2=[0:10:180];
        ph2=[0:10:350];
        [thh,phh]=meshgrid(th2,ph2);
        %** compute the solid angle associated with each grid-point
        dth=180/(length(th2)+1);
        dph=360/(length(ph2)+1);
        %*** scroll through the polar angles
        for jph=1:size(thh,2);
            %**** => for each polar angle compute the band over which this value applies - check the edges are 0 and 180
            th_min=thh(1,jph)-dth/2;
            th_max=thh(1,jph)+dth/2;
            if th_min<0; 
                th_min=th_min+dth/2; 
            end
            if th_max>pi; 
                th_max=th_max-dth/2; 
            end
            %**** => for each band avergae over interpolated gridpoints 
            tmp=sum(sin(pi/180*linspace(th_min,th_max,1e3)))/1e3;
            %**** => store a matrix of solid angles
            dS(:,jph)=tmp*ones(size(thh,1),1);
        end
        SldAngl=dth*dph*dS*(pi/180)^2;
        %** scroll through grid and assign values for each grid point
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
            vDRR(:,jgrd)=sum(DRR(df_ndx,:).*((180-df(:))*ones(1,Nf)),1)/sum(180-df);
        end

        %** scroll through frequencies and interpolate calibration measures for each frequency
        V=C(Drct_ndx); D=C(Drct_ndx); 
        NaNcols=sum(isnan(vDRR));
        vDRR=vDRR(:,find(NaNcols==0));
        SldAngl=SldAngl(find(NaNcols==0));
        V.DRR=sum(vDRR.*(ones(Nf,1)*SldAngl(:).'),2)/sum(SldAngl(:));
        if length(C)==1;
            V.DRR=mean(V.DRR)*ones(size(V.DRR));
        end
    %	%** save plots of calibration information
        clear C
        C(1)=D; 
        C(1).Name='Direct';
        C(2)=V;
        C(2).Name='Omnidirectional';
        for jj=1:2;
            C(jj).Path=OtPth;
        end
        save(sprintf('%s/H.mat',OtPth),'C');
        %*** Plot direct vs volume DRRs
        figure(1);
        plot(D.DRR,D.ff);
        hold on
        plot(V.DRR,D.ff,'r--');
        hold off
        legend({'Direct';'Omnidirectional'})
        xlabel('DRR (s)');
        ylabel('Frequency (kHz)');
        title([OtPth ': Cochleagram']);
        saveas(gcf,sprintf('%s/Cgram',OtPth),'jpg');
    end % if ~exist(Precomputed Data File)
else
    C=[];
end % if ~isempty(Dc)
clear H 

%* == Extract IRs == 
%** load IRs
load(Hpth); %H
%** Scroll through them 
for jh=1:length(H);
    tH=H(jh);
    %*** => Analyze
    tH=hPrp(tH,C,Nbnds,flm,Sb_fs);
    %*** => save plots of IR information
    %*** => save audio
    h=tH.nh;
    MaxAmp=max(abs(h));
    h=h/MaxAmp*(1-1e-6);
    h=[zeros(ceil(tH.fs/5),1); h];
    audiowrite(sprintf('%s/h_denoised_%03d.wav',tH.Path,Nbnds),h,tH.fs,'BitsPerSample',24);
    save(sprintf('%s/H_%03d.mat',tH.Path,Nbnds),'tH');
    if jh==1;
        nH=tH;
        %*** save some useful info
        OtPth=sprintf('%s/Cal_%03dBnds_%05d-%05dHz_ft%05dHz',tH.Path,Nbnds,flm(1),flm(2),Sb_fs); % A path to save the calibration files
    else
        nH=[nH tH];
    end
end
% rename
H=nH;

%* == Save details about CPU run time

%=================================
%
%% load raw golay codes
%[tmp,fh]=audioread([Gpth{1}]);
%g=tmp(:,1);
%Npts=length(g);
%
%% load recorded audio
%DAB=dir([Rpth '*.wav'])
%Nrc=length(DAB);
%if Nrc>=1;
%    % iterate through files
%    h=zeros(NGly,Nrc*8,2);         
%    Tshft=zeros(1,2*Nrc*8,2);
%    hcnt=0; snpcnt=[0 0];
%    for jrc=1:Nrc; hcnt=hcnt+1;
%        % load measured data
%        fnm=['RecordedAudio/' DAB(jrc).name];
%        [tmp,fs]=audioread(fnm); if size(tmp,2)==1; tmp=[tmp tmp]; end
%        if length(tmp)>Npts
%            ab=tmp(1:Npts,:);
%        else
%            ab=[tmp; mean(mean(tmp))*10^(-18)*randn(Npts-length(tmp),2)];
%        end            
%        % scroll through left and right channels
%        for jch=1:2
%             fprintf('Computing file %s: %d/%d, channel %d\n',Rpth,jrc,Nrc,jch)
%             % estimate IR from audio
%             [th,shft]=FreqDomainProcessing7(g,ab(:,jch),fs,NGly,SpkrSpc,Nbnds); 
%             % catenate
%             tsnp=size(th,2); 
%             h(:,snpcnt(jch)+[1:tsnp],jch)=th;  Tshft(1,snpcnt(jch)+[1:2*tsnp],jch)=shft;
%             snpcnt(jch)=snpcnt(jch)+tsnp;
%             %plot
%             close force all; pause(1); figure(152); 
%             plt=h(:,1:snpcnt(jch),jch);  plt=20*log10(abs(plt)); imagesc(plt); set(gca,'ylim',[0 1e3]); caxis([max(max(plt))-60 max(max(plt))]); hold on; plot([0 snpcnt(jch)],[(1/2)*fs/340]*ones(1,2),'k--'); title([Rpth ': Golay ' num2str(jrc) '/' num2str(Nrc)]);
%             drawnow
%             close force all; pause(1); figure(153)
%             subplot(2,1,1);plot(Tshft(1,1:2*snpcnt(jch),1)); axis tight; subplot(2,1,2); plot(Tshft(1,1:2*snpcnt(jch),2)); axis tight; xlabel('snapshot'); title('Estimated offset between the Golay A and B sequences')
%             close all
%        end
%    end 
%    % remove everything after the maximum time limit
%    Nmax=ceil(Tlim*fs);
%    h=h(1:Nmax,:,:);   
%    
%    % roughly assess the quality of each snapshot with the ratio of the
%    % peak to standard deviation
%    pk2sg=max(abs(h),[],1)./std(h,0,1);
%    for jch=1:2
%        % sort in order of decreasing goodness
%        [~,Isrt]=sort(pk2sg(1,:,jch),'descend');
%        hsrt=h(:,Isrt,jch);
%        % see how this 'goodness' changes if we remove the worst offenders
%        % before computing the mean
%% % %       ctff=snpcnt(jch); tp2so=0; tp2sn=1;
%% % %       while (tp2sn>tp2so&&ctff>1);
%% % %            ctff=ctff-1;
%% % %            thn=mean(hsrt(:,1:ctff),2);         tho=mean(hsrt(:,1:ctff+1),2);     
%% % %            tp2sn=max(abs(thn))./std(thn)*sqrt(ctff);      tp2so=max(abs(tho))./std(tho)*sqrt(ctff+1);
%% % %       end
%        ctff=size(hsrt,2)-1;
%        hopt=hsrt(:,1:ctff+1);   
%        % average across snapshots for the best possible estimate and
%        % compute properties
%        hbst=mean(hopt,2);
%        %eval(['!mkdir -p ' Rpth ';']);
%        %eval(['!mkdir -p ' Rpth '/ch' num2str(jch) ';']);
%        mkdir(Rpth);
%        mkdir([Rpth '/ch' num2str(jch)]);
%        Hbst=IRprp(hbst,fs,Nbnds,[20 2e4],SpkrSpc,100,1,1,[Rpth '/ch' num2str(jch)]);
%        Hbst.hop=hopt;
%        H(jch)=Hbst;
%    end
%    
%    if SvFLG==1;
%         fprintf('saving to %s\n',Rpth)
%         save([Rpth '/H'],'H','h','Rpth')
%         % print time series
%         close force all; pause(1); figure(199);
%         for jch=1:2; subplot(2,2,jch); plot(H(jch).tt,H(jch).h); axis([0 2 -1 1]);  hold on; subplot(2,2,2+jch); plot(H(jch).tt,20*log10(abs(H(jch).h))); axis([0 2 -80 0]); hold on;  end
%         print(gcf,'-dpng',[Rpth '/ts'])
%         % write audio
%         It0=ceil(2*max((H(1).Rtt(:,end)))*H(1).fs); if It0==0; It0=length(H(1).h); end; if (It0>=length(H(1).h)); It0=length(H(1).h)-10; end;
%         audiowrite([Rpth '/h0L.wav'],H(1).h(1:It0)/max(abs(H(1).h(1:It0)))* ...
%                    0.999,fh,'BitsPerSample',24)
%         if length(H)>1;
%          It0=ceil(2*max((H(2).Rtt(:,end)))*H(2).fs); if isempty(It0); It0=length(H(2).h); end; if (It0>=length(H(2).h)); It0=length(H(2).h)-10; end;
%          audiowrite([Rpth '/h0R.wav'],H(2).h(1:It0)/max(abs(H(2).h(1:It0)))* ...
%                     0.999,fh,'BitsPerSample',24)
%        end
%         % and the cleaned versions
%         for jch=1:2; subplot(2,2,jch); plot(H(jch).tt,H(jch).nh,'r'); axis([0 2 -1 1]); subplot(2,2,2+jch); plot(H(jch).tt,20*log10(abs(H(jch).nh)),'r'); axis([0 2 -80 0]);  end
%         print(gcf,'-dpng',[Rpth '/ts'])
%         It0=ceil(2*max((H(1).Rtt(:,end)))*H(1).fs); if It0==0; It0=length(H(1).nh); end; if (It0>=length(H(1).nh)); It0=length(H(1).nh)-10; end;
%         audiowrite([Rpth '/nh0L.wav'],H(1).nh(1:It0)/max(abs(H(1).nh(1:It0)))*0.999,fh,'BitsPerSample',24)
%         if length(H)>1;
%             It0=ceil(2*max((H(2).Rtt(:,end)))*H(2).fs); if isempty(It0); It0=length(H(2).nh); end; if (It0>=length(H(2).nh)); It0=length(H(2).nh)-10; end;
%              audiowrite([Rpth '/nh0R.wav'],H(2).nh(1:It0)/max(abs(H(2).nh(1:It0)))*0.999,fh,'BitsPerSample',24)
%         end
%    end 
%end
