%* == IR_Extract.m == 
% Extracts the IRs from recordings of Golay sequence broadcast into a space and re-recorded. This code only outputs audio files of IRs.  See IR_Analysis.m to analyze them or to calibrate them to remove the speaker-microphoen transfer function.
% This code makes use of the following functions:
% - hExtrct.m	: extracts the IR time series
% - GtMtDt.m	: Reads metadata about recording from a textfile

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')

%* == Specify Inputs == 

%** = Path to recording of recorded broadcast =
% (wildcards accepted to process multiple files in a single run)
Rpth='RecordedAudio/*Hallway*.wav'
%Rpth='CalibrationRecordings'
%** = Path to Golay Code used in the broadcast =
Gpth='RawGolay';
%** = Name of golay code =
Gnm='golay_N19_44kHz_24bits';
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
for jh=1:length(Dh);
	%*** => get filename 
	fnm=Dh(jh).name(1:end-4);
	eval(spriintf('!mkdir -p %s/%s',PthStm,fnm))
	%*** => read the audio file
	[rc,fr]=audioread(sprintf('%s/%s.wav',PthStm,fnm));
	%*** => Check golay and recorded audio have the same sampling frequency
	ChckSm(fr,G.fs,'fs of Golay and recording');
	%*** => Extract
	for jch=1:2;
		tH(jch)=hExtrct(rc(:,jch),G);
		tH.Name=fnm;
		tH.Channel=jch;
		%*** => Load (or query) MetaData
		tH=GtMtDt(tH,sprintf('%s/%s/Meta.txt',Cpth,cnm),Mt);
		%*** => save audio
		tpth=sprintf('%s/%s/ch%d',jch);
		eval(sprintf('!mkdir -p %s',tpth))
		audiowrite(sprintf('%s/h.wav',tpth),tH.h,tH.fs,'BitsPerSample',24));
		%*** => save plot of time series
		figure(1);
		h=tH.h;
		plot([1:length(h)]/tH.fs,sign(h).*abs(sign).^0.3);
		xlabel('Time (s)');
		ylabel('Compressed amplitude');
		title(fnm)
		set(gca,'xscale','log');
		saveas(gcf,sprintf('%s/ts.epsc',tpth,fnm));
		set(gca,'xscale','lin');
		saveas(gcf,sprintf('%s/ts_LinTime.epsc',tpth));
		figure(2);
		plot([1:length(h)]/tH.fs,20*log10(abs(h)));
		xlabel('Time (s)');
		ylabel('Amplitude (dB)');
		title(fnm)
		set(gca,'xscale','log');
		saveas(gcf,sprintf('%s/ts_dB.epsc',tpth));
		set(gca,'xscale','lin');
		saveas(gcf,sprintf('%s/ts_dB_LinTime.epsc',tpth));
	end
end % jh=1:length(Dh);

%* == Save details about CPU run time
