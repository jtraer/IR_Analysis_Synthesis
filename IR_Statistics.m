%* == IR_Statistics.m == 
% Load mat-files of analyzed IRs and plots their statistics
% This makes use of the following functions:
% - 

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')
set(0,'DefaultFigureVisible','On');

%* == Specify Inputs == 
Input_file='Input_IR_Survey_2';
eval(sprintf('[Rpth,Cpth,Mt]=%s;',Input_file));

Nbnds=4;

pth=pwd;
if strcmp(pth(1:3),'/om')
    tmtpth='../timit'
else
    tmtpth='/Users/jtraer/LabBook/Workshop/PerceptualExperiments/Sandbox/timit';
end

%** Specify Rejection Criteria
rcnt=0;
rcnt=rcnt+1; Rjct(rcnt).Expr='datenum(H.DateCreated)<datenum(''12-May-2016'')';

%* ==== Load data  ====

%** clear the output folders
eval('! rm -rf IRMAudio/*');
%** load IR data
fprintf('Loading data... \n'); dstrng=[]; tic;

%** scroll through IR folders and create a structure of Path stems
Dh=[];
for jr=1:length(Rpth);
    tDh=dir([Rpth(jr).Pth '/*.wav']);
    for jj=1:length(tDh);
        ttDh=dir(sprintf('%s/%s/ch*',Rpth(jr).Pth,tDh(jj).name(1:end-4)));
        for jch=1:length(ttDh);
            t3Dh=dir(sprintf('%s/%s/ch%d/H_%03dbnds.mat',Rpth(jr).Pth,tDh(jj).name(1:end-4),jch,Nbnds)); 
            %==> save path stems
            if length(t3Dh)>0;
                t3Dh.PthStm=sprintf('%s/%s/ch%d',Rpth(jr).Pth,tDh(jj).name(1:end-4),jch); 
                Dh=[Dh; t3Dh];
            end
        end
    end
end

%**  remove IRs we don't want to analyze today 
for jr=1:length(Rjct);
    cnt=0; rcnt=0;
    while cnt<length(Dh); cnt=cnt+1;
        load(sprintf('%s/%s',Dh(cnt).PthStm,Dh(cnt).name))
        eval(sprintf('if %s; Dh(cnt)=[]; cnt=cnt-1; rcnt=rcnt+1; end',Rjct(jr).Expr));
    end
    fprintf('%d/%d IRs rejected for %s\n',rcnt,length(Dh)+rcnt,Rjct(jr).Expr);
end

%* Label the Path structures
for jh=1:length(Dh)
    load(sprintf('%s/%s',Dh(jh).PthStm,Dh(jh).name))
    M=GtMtDt(sprintf('%s/Meta.txt',Dh(jh).PthStm(1:end-4)),Mt);
    M=orderfields(M);
    Dh(jh).Meta=M.Meta;
    %Mflds=fields(M);
    %for jfld=1:length(Mflds)
    %    if ~strcmp(Mflds{jfld},'Path')
    %        eval(sprintf('H(jh).%s=M.%s;',Mflds{jfld},Mflds{jfld}));
    %    end
    %end
end

%** Normalize amplitudes
for jh=1:length(Dh)
    load(sprintf('%s/%s',Dh(jh).PthStm,Dh(jh).name))
    Dh(jh).MaxAmp=H.MaxAmp;
end
Mx=max([Dh.MaxAmp]);
for jh=1:length(Dh);
    Dh(jh).MaxAmp=Dh(jh).MaxAmp/Mx;
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
hPltStts(Dh,Mt);
%* == Write an html file to display all the data
%** clear the output folders
WrtDt2HTML(Dh,'IRMAudio/IRdata',Mt,tmtpth);
fprintf('Data written to:\n\n %s/IR_Data_Summary.html\n\n',pwd)
unix('cp IR_Data_Summary.html IRMAudio/')
unix('! zip -r IRMAudio.zip IRMAudio/*')

%* == TODO: Collect details about code runtime
%* == TODO: Save this code to a summary file
%eval(sprintf('! grep "%%\\*" %s.m > tmp.org',cfl))
%eval('! sed ''s/^ *//g'' < tmp.org > tmp2.org');  % remove whitespace
%eval('! sed ''s/^[ \\t]+//g'' < tmp2.org > tmp.org');  % remove tabs too
%eval(sprintf('! sed ''s/^.//'' tmp.org > %s.org',cfl))   % remove '%' so emacs can read the indenting

%* == TODO: Save and Archive everything

