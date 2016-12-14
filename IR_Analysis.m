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
Hpth='H_raw_TryBox_1IRs_12-Dec-2016'
%Hpth='H_raw_FirstTest_Marble_1IRs_06-Dec-2016'
%** = File with Calibration IRs (Optional, but recommended) =
% These are re-recorded broadcasts (as in the measurements) but in as close to anechoic conditions as possible. These are required to ensure we are recording the properties of the space, not of the speaker/microphone/soundcard.
Cpth='H_raw_Cal_1IRs_12-Dec-2016';
%** = Number of cochlear subbands for analysis =
Nbnds=[30];
%** = Frequency limits in Hz =
flm=[100 20e3];
%** = Frequency of subband envelopes in Hz =
Sb_fs=100;

%** = Specify MetaData we want to record (these can be added or removed arbitrarily) =
%mcnt=0;
%mcnt=mcnt+1;Mt{mcnt}='App.Mic';
%mcnt=mcnt+1;Mt{mcnt}='App.Recorder';
%mcnt=mcnt+1;Mt{mcnt}='App.Gain';
%mcnt=mcnt+1;Mt{mcnt}='App.Speaker';
%mcnt=mcnt+1;Mt{mcnt}='App.Volume';
%mcnt=mcnt+1;Mt{mcnt}='Env.Class';
%mcnt=mcnt+1;Mt{mcnt}='Env.Size';
%mcnt=mcnt+1;Mt{mcnt}='Env.Material';

%* == Calibrate Apparatus ==
%** search for calibration files 
Dc=dir(sprintf('%s*.mat',Cpth));
%** scroll through the available files
if ~isempty(Dc)
	ccnt=0;
	for jc=1:length(Dc);
        load(Dc(jc).name);
        tC=H;
        %*** => Check to see if the IR is already analyzes...
        %*** TODO add this so we don't repeat analyses unnecessarily
        %*** => Analyze
        tC=hPrp(tC,[],Nbnds,flm,Sb_fs);
        
        
        if jc==1;
            C=tC;
        else
            C=[C tC];
        end
%		%*** => if metadata exists open it, otherwise query user
%		H=GtMtDt(H,sprintf('%s/%s/Meta.txt',Cpth,cnm),'Mic','Gain','Speaker','Vol','Recorder','PolarAngle','Azimuth','Distance');
%		%*** => add to a structure of all calibration measurements
%		C(ccnt)=H;
	end % for jc=1:length(Dc)
%	%** interpolate spatial power spectrum
%	th=[C.PolarAngle];
%	ph=[C.Azimuth];
%	%*** scroll through all calibration measurements
%	for jc=1:length(C)
%		%*** => record subband power for each measurement
%		Pwr(jc,:)=C(jc).P_k.';	
%	end
%	%*** scroll through subbands
%	nth=linspace(0,180,100);
%	nph=linspace(-180,180,200);
%	for jbnd=1:(Nbds+2);
%		%*** => interpolate a map of directional power
%		SpkrPwr(:,:,jbnd)=interp2(th,ph,Pwr(:,jbnd),nth,nph);
%	end
%	%** save plots of calibration information
%	%*** Scroll through calibration measurements
%	%*** plot Equipment IRs
%	%*** 
end % if ~isempty(Dc)

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
    audiowrite(sprintf('%s/h_denoised_%03d.wav',tH.Path,Nbnds),h,tH.fs,'BitsPerSample',24);
    save(sprintf('%s/H_%03d.mat',tH.Path,Nbnds),'tH');
    if jh==1;
        H=tH;
    else
        H=[H tH];
    end
end


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
