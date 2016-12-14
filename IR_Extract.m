%* == IR_Extract.m == 
% Extracts the IRs from recordings of Golay sequence broadcast into a space and re-recorded. This code only outputs audio files of IRs.  See IR_Analysis.m to analyze them or to calibrate them to remove the speaker-microphoen transfer function.
% This code makes use of the following functions:
% - hExtrct.m	: extracts the IR time series
% - GtMtDt.m	: Reads metadata about recording from a textfile

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')

%* == Specify Inputs == 

%** = Name =
Nm='TryBox';
%Nm='Cal';
%** = Path to recording of recorded broadcast =
% (wildcards accepted to process multiple files in a single run)
Rpth='RecordedAudio/*FR6*.wav';
%Rpth='CalibrationRecordings/*Woofit*Rode*.wav'
%** = Path to Golay Code used in the broadcast =
Gpth='RawGolay';
%** = Name of golay code =
Gnm='golay_44kHz_N16_2min_24bits';
%** = Specify MetaData we want to record (these can be added or removed arbitrarily) =
mcnt=0;
mcnt=mcnt+1;Mt{mcnt}='App.Mic';
mcnt=mcnt+1;Mt{mcnt}='App.Recorder';
mcnt=mcnt+1;Mt{mcnt}='App.Gain';
mcnt=mcnt+1;Mt{mcnt}='App.Speaker';
mcnt=mcnt+1;Mt{mcnt}='App.Volume';
mcnt=mcnt+1;Mt{mcnt}='Env.Class';
mcnt=mcnt+1;Mt{mcnt}='Env.Size';
mcnt=mcnt+1;Mt{mcnt}='Env.Material';

%* == Load golay sequence ==
load(sprintf('%s/%s.mat',Gpth,Gnm)); %G

%* == Extract IRs == 
%** Find all audio files
Dh=dir(Rpth);
%** Save stem of path 
PthStm=GtPthStm(Rpth);
%** Scroll through them
hcnt=0;
for jh=1:length(Dh);
	%*** => get filename 
	fnm=Dh(jh).name(1:end-4);
    Fllnm=sprintf('%s/%s',PthStm,fnm);
	eval(sprintf('!mkdir -p %s',Fllnm));
    fprintf('Extracting %s\n',Fllnm);
	%*** => read the audio file
	[rc,fr]=audioread(sprintf('%s/%s.wav',PthStm,fnm));
	%*** => Check golay and recorded audio have the same sampling frequency
	ChckSm(fr,G.fs,'fs of Golay and recording');
	%*** => Extract
	for jch=1:size(rc,2); hcnt=hcnt+1;
        Fllnm_ch=sprintf('%s/ch%d',Fllnm,jch);
        eval(sprintf('!mkdir -p %s',Fllnm_ch));
		tH=hExtrct(rc(:,jch),G,Fllnm_ch);
		tH.Name=fnm;
		tH.Path=sprintf('%s/ch%d',Fllnm,jch);
		tH.Channel=jch;
		%*** => Load (or query) MetaData
		M=GtMtDt([Fllnm '/Meta.txt'],Mt);
        %*** => re-order fields for consistency
        M=orderfields(M);
        tH.Meta=M;
		%*** => normalize and save audio
        h=tH.h;
        tH.MaxAmp=max(abs(h));
        h=h/tH.MaxAmp*(1-1e-6);
		audiowrite(sprintf('%s/h.wav',Fllnm_ch),h,tH.fs,'BitsPerSample',24);
		%*** => save structure
        save(sprintf('%s/H.mat',Fllnm_ch),'tH');
		%*** => save plot of time series
		%figure(1);
		%h=tH.h;
		%plot([1:length(h)]/tH.fs,sign(h).*abs(h).^0.3);
		%xlabel('Time (s)');
		%ylabel('Compressed amplitude');
		%title(Fllnm)
		%set(gca,'xscale','log');
		%saveas(gcf,sprintf('%s/ts.jpg',Fllnm_ch));
		%set(gca,'xscale','lin');
		%saveas(gcf,sprintf('%s/ts_LinTime.epsc',tpth));
		%figure(2);
		%plot([1:length(h)]/tH.fs,20*log10(abs(h)));
		%xlabel('Time (s)');
		%ylabel('Amplitude (dB)');
		%title(fnm)
		%set(gca,'xscale','log');
		%saveas(gcf,sprintf('%s/ts_dB.epsc',tpth));
		%set(gca,'xscale','lin');
		%saveas(gcf,sprintf('%s/ts_dB_LinTime.epsc',tpth));
        
        %** => add to big structure
        H(hcnt)=tH;
        clear tH
	end % jch
end % jh=1:length(Dh);

%* == Save all data ==
save(sprintf('H_raw_%s_%dIRs_%s',Nm,length(H),date),'H')
fprintf('Data saved to H_raw_%s_%dIRs_%s\n',Nm,length(H),date);

%* == Save details about code and CPU run time ==

SummarizeCode(mfilename('fullpath'))
