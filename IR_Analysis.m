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
%Hpth='RecordedAudio/RoomReverb/H_raw_RoomRvrb_RmDst_24IRs_08-May-2017'; Cpth='CalibrationRecordings/H_raw_CAL_ZP-TSCM_12IRs_08-May-2017';
%Hpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDst_14IRs_08-May-2017'; Cpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDstCAL_12IRs_08-May-2017';
%Hpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDst_30IRs_08-May-2017'; Cpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDstCAL_12IRs_08-May-2017';
%Hpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDst_30IRs_09-May-2017'; Cpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDstCAL_12IRs_08-May-2017'; % test with th other Golay code

Hpth='RecordedAudio/RoomReverb/4078/H_raw_Office_Test_14IRs_16-May-2017'; Cpth='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDstCAL_12IRs_08-May-2017'; % test with th other Golay code

%Hpth='H_raw_LED4-Center_21IRs_05-Feb-2017'; Cpth='RecordedAudio/Boards/LED4-Center/CAL/H_raw_CAL_LED4_1IRs_05-Feb-2017';

%Hpth='RecordedAudio/Tweeter/H_raw_Twt_2IRs_22-Mar-2017'; Cpth='';
%Hpth='RecordedAudio/Boards/LED4-Center/SbSt/H_raw_LED4-Center_Sb_9IRs_25-Apr-2017'; Cpth='RecordedAudio/Boards/LED4-Center/CAL/H_raw_CAL_LED4_1IRs_22-Apr-2017';
%Hpth='RecordedAudio/Balls/Smll/H_raw_Balls_22IRs_22-Apr-2017.mat'; Cpth='RecordedAudio/Balls/Smll/H_raw_Balls_CAL_3IRs_22-Apr-2017';
%Hpth='RecordedAudio/Boards/LED4-Center/CAL/H_raw_CAL_LED4_1IRs_05-Feb-2017'; Cpth='';

%** = Number of cochlear subbands for analysis =
Nbnds=[30];
%** = Frequency limits in Hz =
flm=[50 20e3];
%** = Frequency of subband envelopes in Hz =
Sb_fs=1e4;
%** filetype
%ftp='epsc';
ftp='jpg';

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
            tC=hPrp(tC,[],Nbnds,flm,Sb_fs,ftp);
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
    %keyboard
    tH=hPrp(tH,C,Nbnds,flm,Sb_fs,ftp);
    %*** => Make synthetic IRs
    %tH=hSynth(tH,ftp);
    %*** => save plots of IR information
    %*** => save audio
    h=tH.nh;
    MaxAmp=max(abs(h));
    h=h/MaxAmp*(1-1e-6);
    h=[zeros(ceil(tH.fs/5),1); h; zeros(ceil(tH.fs/2),1)];
    audiowrite(sprintf('%s/h_denoised_%03d.wav',tH.Path,Nbnds),h,tH.fs,'BitsPerSample',24);
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
% save
save(sprintf('%s_%dBnds',Hpth,Nbnds),'H','C')   
fprintf('Data saved to %s_%dBnds\n',Hpth,Nbnds);   

figure(200)
for jh=1:length(H);
    CoM=mean(H(jh).RT60.*[1:Nbnds])/mean(H(jh).RT60);
    hold on;
    plot(1,CoM,'o');
    text(1,CoM,1.001,sprintf(H(jh).Name));
end


%* == Save details about CPU run time

