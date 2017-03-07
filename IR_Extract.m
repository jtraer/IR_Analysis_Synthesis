%* == IR_Extract.m == 
% Extracts the IRs from recordings of Golay sequence broadcast into a space and re-recorded. This code only outputs audio files of IRs.  See IR_Analysis.m to analyze them or to calibrate them to remove the speaker-microphoen transfer function.
% This code makes use of the following functions:
% - hExtrct.m	: extracts the IR time series
% - GtMtDt.m	: Reads metadata about recording from a textfile

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')

%* == Specify Inputs == 

%** = Name and Path to recording of recorded broadcast =
% (wildcards accepted to process multiple files in a single run)

Rpth='RecordedAudio/Boards/Crt-Edge/CAL/*.wav'; Nm='CAL_Crt';
Rpth='RecordedAudio/Boards/Crt-Edge/*.wav'; Nm='Crt-Edge';
Rpth='RecordedAudio/Boards/LED4-Center/CAL/*.wav'; Nm='CAL_LED4';
Rpth='RecordedAudio/Boards/LED4-Center/CAL/*.wav'; Nm='CAL_LED4';
Rpth='RecordedAudio/Boards/LED4-Center/SbSt/*.wav'; Nm='LED4-Center_Sb';
%Rpth='RecordedAudio/Boards/LED8-Corner/CAL/*.wav'; Nm='CAL_LED8';
%Rpth='RecordedAudio/Boards/LED8-Corner/*.wav'; Nm='LED8-Corner';
%Rpth='RecordedAudio/Boards/Rd-Ext/CAL/CAL*.wav'; Nm='CAL_Rd';
%Rpth='RecordedAudio/Boards/Rd-Ext/*.wav'; Nm='Rd-Ext';

%Rpth='RecordedAudio/RoomReverb/Tst/*.wav'; Nm='RoomRvrb_tst';
%Rpth='RecordedAudio/RoomReverb/Distance_RoomSize/*.wav'; Nm='RoomRvrb_RmDst';
%Rpth='RecordedAudio/RoomReverb/Location/*.wav'; Nm='RoomRvrb_Loc';
%Rpth='RecordedAudio/RoomReverb/Empty_vs_Full/*.wav'; Nm='RoomRvrb_EmptyOrFull';
%Rpth='CalibrationRecordings/ZIPP-TASCAM/*.wav'; Nm='CAL_ZP-TSCM';

%** = Path to Golay Code used in the broadcast =
Gpth='RawGolay';
%** = Name of golay code =
Gnm='golay_44kHz_N16_2min_24bits';
%Gnm='golay_44kHz_N19_3min_24bits';
%Gnm='golay_44kHz_N16_1min_24bits';
%Gnm='golay_44kHz_N16_3min_24bits';
%** = Specify MetaData we want to record (these can be added or removed arbitrarily) =
mcnt=0;
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Mic';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Recorder';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Gain';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Speaker';
%mcnt=mcnt+1;Mt{mcnt}='Meta.App.Volume';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Class';


%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Size';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Location';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Distance';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.NoPeople';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Door';

%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Size';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Material';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Location';
%mcnt=mcnt+1;Mt{mcnt}='Meta.Env.Damping';

mcnt=mcnt+1;Mt{mcnt}='Meta.App.PolarAngle_fromTop';
mcnt=mcnt+1;Mt{mcnt}='Meta.App.AzimuthalAngle_fromFront';


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
        %tH=rmfield(tH,'h_snps');
		tH.Name=fnm;
		tH.Path=sprintf('%s/ch%d',Fllnm,jch);
		tH.Channel=jch;
		%*** => Load (or query) MetaData
		M=GtMtDt([Fllnm '/Meta.txt'],Mt);
        %*** => re-order fields for consistency
        M=orderfields(M);
        tH.Meta=M.Meta;
        tH.Meta.Path=M.Path;
		%*** => normalize and save audio
        h=tH.h;
        tH.MaxAmp=max(abs(h));
        h=h/tH.MaxAmp*(1-1e-6);
        h=[zeros(ceil(tH.fs/5),1); h];
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
