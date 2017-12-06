%* == IR_Statistics.m == 
% Load mat-files of analyzed IRs and plots their statistics
% This makes use of the following functions:
% - 

%* == Preamble ==
clear all; close all; clc
path(path,'Tools')
set(0,'DefaultFigureVisible','Off');

%* == Specify Inputs == 
Input_file='Input_Survey_2'; Nm='Tst';
Input_file='Input_IRSurvey_NatStats'; Nm='NtStts';
%Input_file='Input_IR_Survey_2_OM'; Nm='RvrbStrct'
%Input_file='Input_IR_ControlION'; Nm='CntrlION'
%Input_file='Input_Spkr'; Nm='Spkr';
%Input_file='Input_ShrtvsLng';  Nm='ShrtvsLng'
%Input_file='Input_IR_Control'; Nm='CntrlZpp'
%Input_file='Input_ObjIRs'; Nm='ObjIRs'
%Input_file='Input_ObjIRs_Ext'; Nm='ObjIRs_Ext'
eval(sprintf('[Rpth,Cpth,Mt,Amnd,html_tmp]=%s;',Input_file));

Nbnds=4;
H_FLG=1; % if 1 this copies the H.mat files into the output folder (warning this might make a huge and unwieldy directory)

pth=pwd;
if strcmp(pth(1:3),'/om')
    tmtpth='../timit'
else
    tmtpth='~/multigore/Labbook/Projects/ReverbAndPerception/timit';
end

%** Specify Rejection Criteria
rcnt=0;
rcnt=rcnt+1; Rjct(rcnt).Expr='datenum(H.DateCreated)<datenum(''23-July-2015'')';
%rcnt=rcnt+1; Rjct(rcnt).Expr='length(H.Attck(3).Spc)~=256';

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
for jr=1:length(Cpth);
    tDh=dir([Cpth(jr).Pth '/*.wav']);
    % add the omni file
    ttDh=tDh(1);
    ttDh.name='Omni.wav';
    tDh(length(tDh)+1)=ttDh;
    for jj=1:length(tDh);
        ttDh=dir(sprintf('%s/%s/ch*',Cpth(jr).Pth,tDh(jj).name(1:end-4)));
        for jch=1:length(ttDh);
            t3Dh=dir(sprintf('%s/%s/ch%d/H_%03dbnds.mat',Cpth(jr).Pth,tDh(jj).name(1:end-4),jch,Nbnds)); 
            %==> save path stems
            if length(t3Dh)>0;
                t3Dh.PthStm=sprintf('%s/%s/ch%d',Cpth(jr).Pth,tDh(jj).name(1:end-4),jch); 
                Dh=[t3Dh; Dh];
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
fprintf('%d IRs left\n',length(Dh));

%* Label the Path structures
for jh=1:length(Dh)
    load(sprintf('%s/%s',Dh(jh).PthStm,Dh(jh).name))
    M=GtMtDt(sprintf('%s/Meta.txt',Dh(jh).PthStm(1:end-4)),Mt);
    M=orderfields(M);
    Dh(jh).Meta=M.Meta;
    % add a Name
    %Dh(jh).Meta.FileName=Dh(jh).PthStm;
    %Mflds=fields(M);
    %for jfld=1:length(Mflds)
    %    if ~strcmp(Mflds{jfld},'Path')
    %        eval(sprintf('H(jh).%s=M.%s;',Mflds{jfld},Mflds{jfld}));
    %    end
    %end
end
%Mt{length(Mt)+1}='Meta.FileName';

%** Normalize amplitudes
for jh=1:length(Dh)
    load(sprintf('%s/%s',Dh(jh).PthStm,Dh(jh).name))
    Dh(jh).MaxAmp=H.MaxAmp;
end
Mx=max([Dh.MaxAmp]);
for jh=1:length(Dh);
    Dh(jh).MaxAmp=Dh(jh).MaxAmp/Mx;
end

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
fNm=sprintf('IRstts_%s_%03d',Nm,Nbnds);
unix(sprintf('mkdir -p %s',fNm));
unix(sprintf('rm -rf %s/*',fNm));
hPltStts(Dh,Mt,Amnd,fNm);

%Mt=Mt(1:end-1); % this is because WrtDt2HTML bugs if we have FileName in it

%* == Write an html file to display all the data
%** clear the output folders
%fcnt=0;
%fcnt=fcnt+1; Flds{fcnt}='Meta.Env.Class';
%fcnt=fcnt+1; Flds{fcnt}='Meta.Env.SpaceName';
WrtDt2HTML(Dh,sprintf('%s',fNm),html_tmp,sprintf('%s',fNm),Mt,Mt,tmtpth,H_FLG);
if length(Cpth)>0;
    unix(sprintf('cp %s/*.jpg %s/',Cpth(1).Pth,fNm));
end
%unix(sprintf('cp IR_Data_%s.html IRMAudio/',Nm));
%unix(sprintf('mv IRMAudio IRMAudio_%s_%03dbnds',Nm,Nbnds));
unix(sprintf('rm %s.zip',fNm));
unix(sprintf('zip -r %s.zip %s/*',fNm,fNm));
fprintf('Data written to:\n\n %s/%s.zip\n\n',pwd,fNm)

%* == TODO: Collect details about code runtime
%* == TODO: Save this code to a summary file
%eval(sprintf('! grep "%%\\*" %s.m > tmp.org',cfl))
%eval('! sed ''s/^ *//g'' < tmp.org > tmp2.org');  % remove whitespace
%eval('! sed ''s/^[ \\t]+//g'' < tmp2.org > tmp.org');  % remove tabs too
%eval(sprintf('! sed ''s/^.//'' tmp.org > %s.org',cfl))   % remove '%' so emacs can read the indenting

%* == TODO: Save and Archive everything

