%* == IR_Statistics.m == 
% Load mat-files of analyzed IRs and plots their statistics
% This makes use of the following functions:
% - 

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')
set(0,'DefaultFigureVisible','On');

%* == Specify Inputs == 

%** = Path to IRs =
% (wildcards accepted to process multiple files in a single run)
pcnt=0;
%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/RoomReverb/H_raw_RoomRvrb_tst_12IRs_26-Feb-2017_5Bnds';
%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/RoomReverb/H_raw_RoomRvrb_RmDst_24IRs_27-Mar-2017_30Bnds';
%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/FullSurvey/H_raw_Srvy_390IRs_29-Mar-2017_5Bnds';
%%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDst_30IRs_08-May-2017_10Bnds';
pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/RoomReverb/OfficeDistanceTest/H_raw_OffcDst_30IRs_09-May-2017_10Bnds';
pltcnt=0;
pltcnt=pltcnt+1; PltPrms{pltcnt}='Meta.Env.Distance';
%pltcnt=pltcnt+1; PltPrms{pltcnt}='Meta.Env.Size';

%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/RoomReverb/H_raw_RoomRvrb_Loc_18IRs_27-Feb-2017_35Bnds';
%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/RoomReverb/H_raw_RoomRvrb_Loc_26IRs_27-Feb-2017_35Bnds';

%pcnt=pcnt+1; IRpth{pcnt}='RecordedAudio/Boards/LED4-Center/SbSt/H_raw_LED4-Center_Sb_9IRs_25-Apr-2017_40Bnds';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Material';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.YoungsMod';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Density';

%** Specify Rejection Criteria
rcnt=0;
rcnt=rcnt+1; Rjct(rcnt).Expr='datenum(H(cnt).DateCreated)<datenum(''12-May-2016'')';
%** Specify Paramters to be plotted
pcnt=0;
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Size';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Distance';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Location';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.NoPeople';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Door';
%pcnt=pcnt+1; PltPrms{pcnt}='Meta.Env.Class';

%* ==== Load data  ====

%** clear the output folders
eval('! rm -rf IRMAudio/*');
%** load IR data
fprintf('Loading data... \n'); dstrng=[]; tic;

%** scroll through IR folders
tH=[]; 
for jh=1:length(IRpth);
    load(IRpth{jh});
    %H=orderfields(H);
    % May need to do dramatic cleaning here
    tH=[tH; H];
end
H=tH; clear tH;

%**  remove IRs we don't want to analyze today 
for jr=1:length(Rjct);
    cnt=0; rcnt=0;
    while cnt<length(H); cnt=cnt+1;
        eval(sprintf('if %s; H(cnt)=[]; cnt=cnt-1; rcnt=rcnt+1; end',Rjct(jr).Expr));
    end
    fprintf('%d/%d IRs rejected for %s\n',rcnt,length(H)+rcnt,Rjct(jr).Expr);
end

%* Check if all the IRs are labelled correctly
for jh=1:length(H)
    M=GtMtDt(sprintf('%s',H(jh).Meta.Path),PltPrms);
    M=orderfields(M);
    Mflds=fields(M);
    for jfld=1:length(Mflds)
        if ~strcmp(Mflds{jfld},'Path')
            eval(sprintf('H(jh).%s=M.%s;',Mflds{jfld},Mflds{jfld}));
        end
    end
end

%** Normalize amplitudes
Mx=max([H.MaxAmp]);
for jh=1:length(H);
    H(jh).MaxAmp=H(jh).MaxAmp/Mx;
end

%* == add synthetic IRs
%* == Plot Data
%** == TODO: specify a colormap for the global IRs to be used in plots
%*** == TODO: rank and label
%** Run PltStts to plot the statistics of IRs
%*** save these plots
%set(gcf,'inverthardcopy','off');
%set(gcf,'paperpositionmode','auto'); nw=date; %exportfig(gcf,sprintf('IRMFigs/Lgnd%dIRs_%s',length(BH),nw));
%%print(gcf,'-depsc',sprintf('IRMFigs/Lgnd%dIRs_%s',length(BH),nw));
%saveas(gcf,sprintf('IRMFigs/Lgnd%dIRs_%s',length(BH),nw));

%** Write Data
eval('! rm -rf IRMAudio/*');
eval('! mkdir IRMAudio/Audio')
H=hPltStts(H,PltPrms);
%* == Write an html file to display all the data
%** clear the output folders
WrtDt2HTML(H,'IRMAudio/IRdata',PltPrms);
fprintf('Data written to:\n\n %s/IR_Data_Summary.html\n\n',pwd)
eval('! zip IRMAudio/Audio.zip IRMAudio/Audio')

%* == TODO: Collect details about code runtime
%* == TODO: Save this code to a summary file
%eval(sprintf('! grep "%%\\*" %s.m > tmp.org',cfl))
%eval('! sed ''s/^ *//g'' < tmp.org > tmp2.org');  % remove whitespace
%eval('! sed ''s/^[ \\t]+//g'' < tmp2.org > tmp.org');  % remove tabs too
%eval(sprintf('! sed ''s/^.//'' tmp.org > %s.org',cfl))   % remove '%' so emacs can read the indenting

%* == TODO: Save and Archive everything

