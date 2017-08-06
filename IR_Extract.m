%* == IR_Extract.m == 
% Extracts the IRs from recordings of Golay sequence broadcast into a space and re-recorded. This code only outputs audio files of IRs.  See IR_Analysis.m to analyze them or to calibrate them to remove the speaker-microphoen transfer function.
% This code makes use of the following functions:
% - hExtrct.m	: extracts the IR time series
% - GtMtDt.m	: Reads metadata about recording from a textfile

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')

%* == Specify Inputs == 
Input_file='Input_IRSurvey_NatStats';
Input_file='Input_Survey_2';
%Input_file='Input_IR_Survey_2';
Input_file='Input_IR_ControlION'; Nm='CntrlION'
Input_file='Input_IR_Control'; Nm='CntrlZpp'
%Input_file='Input_ACvsBth';
%Input_file='Input_ShrtvsLng';
%Input_file='Input_UtahReverb';
%Input_file='Input_ObjectIRs';
%Input_file='Input_ObjectIRs_Ext';
eval(sprintf('[R,C,Mt]=%s;',Input_file));
%==> R is a structure of recordings
%==> C is a structure of calibration recordings
%==> Mt is a structure of Metadata about the recordings

%* == Extract IRs == 
%** Find all audio files and save a structure of directories
Dh=[];
for jr=1:length(R);
    tDh=dir([R(jr).Pth '/*.wav']);
    for jj=1:length(tDh);
        %==> save path to golay codes 
        tDh(jj).Gpth=R(jr).Gpth;
        %==> save path stems
        %tDh(jj).PthStm=GtPthStm(R(jr).Pth);
        tDh(jj).PthStm=R(jr).Pth;
    end
    Dh=[Dh; tDh];
end
for jc=1:length(C);
    tDh=dir([C(jc).Pth '/*.wav']);
    for jj=1:length(tDh);
        %==> save path to golay codes 
        tDh(jj).Gpth=C(jc).Gpth;
        %==> save path stems
        %tDh(jj).PthStm=GtPthStm(C(jc).Pth);  
        tDh(jj).PthStm=C(jc).Pth;  
    end
    Dh=[tDh; Dh];
end
%** Scroll through recordings
hcnt=0;
for jh=1:length(Dh);
	%*** => get filename 
	fnm=Dh(jh).name(1:end-4);
    Fllnm=sprintf('%s/%s',Dh(jh).PthStm,fnm);
	eval(sprintf('!mkdir -p %s',Fllnm));
    fprintf('Extracting %s\n',Fllnm);
	%*** => read the audio file
	[rc,fr]=audioread(sprintf('%s/%s.wav',Dh(jh).PthStm,fnm));
    %* == Load golay sequence ==
    load(sprintf('%s.mat',Dh(jh).Gpth)); %G
	%*** => Check golay and recorded audio have the same sampling frequency
	ChckSm(fr,G.fs,'fs of Golay and recording');
	%*** => Scroll through channels of recording
	for jch=1:size(rc,2); hcnt=hcnt+1;
        Fllnm_ch=sprintf('%s/ch%d',Fllnm,jch);
        eval(sprintf('!mkdir -p %s',Fllnm_ch));
		H=hExtrct(rc(:,jch),G,Fllnm_ch);
		H.Name=fnm;
		H.Path=sprintf('%s/ch%d',Fllnm,jch);
		H.Channel=jch;
		%*** => Load (or query) MetaData
		M=GtMtDt([Fllnm '/Meta.txt'],Mt);
        if (jh==1&&jch==1);
            jnk=input('Now is a good time to copy the Meta Data file and add the photos\n');
        end
        %*** => re-order fields for consistency
        M=orderfields(M);
        H.Meta=M.Meta;
        H.Meta.Path=M.Path;
		%*** => normalize and save audio
        h=H.h;
        H.MaxAmp=max(abs(h));
        h=h/H.MaxAmp*(1-1e-6);
        h=[zeros(ceil(H.fs/5),1); h];
		audiowrite(sprintf('%s/h.wav',Fllnm_ch),h,H.fs,'BitsPerSample',24);
		%*** => save structure
        save(sprintf('%s/H.mat',Fllnm_ch),'H');
        % ===> Plots
        
	end % jch
end % jh=1:length(Dh);
fprintf('%d IRs extracted\n',hcnt); 

%* == Save details about code and CPU run time ==
