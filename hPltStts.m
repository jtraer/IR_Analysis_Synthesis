function hPltStts(Dh,PltPrms,Amnd,fNm);

ftp='png';
%ftp='epsc';
Txt=0;

% scroll through all the parameters we want to investigate
fcnt=0;
for jPrm=1:length(PltPrms)
    % make a directory
    eval(sprintf('!mkdir -p %s/%s',fNm,PltPrms{jPrm}));

    % define the parameter variables and a colorscheme to plot them
    for jj=1:length(Dh); eval(sprintf('Vr{jj}=Dh(jj).%s;',PltPrms{jPrm})); end
    Vr=unique(Vr);
    Vr(find(strcmp('NA',Vr)))=[]; 
    MrkVr=repmat(['o','s','*','^','v','+'],[1 100]);
    cmpVr=othercolor('Dark28',32);
    clear V
    for jv=1:length(Vr)
        V(jv).name=Vr{jv};
        V(jv).cmp=cmpVr(ceil(length(cmpVr)*jv/length(Vr)),:);
        V(jv).mrk=MrkVr(jv);
    end
    % scroll through the directory and remove files as necessary
    tDh=Dh;
    ndx=[]; for ja=1:length(Amnd); if strcmp(Amnd(ja).Prm,PltPrms{jPrm}); ndx=ja; end; end
    if ~isempty(ndx)
        for jv=1:length(Amnd(ndx).Var);
            hcnt=0; rcnt=0;
            while hcnt<length(tDh); hcnt=hcnt+1; Flg=0;
                eval(sprintf('vr=tDh(hcnt).%s;',Amnd(ndx).Var(jv).Exp));
                if isstr(vr);
                    Flg=strcmp(vr,Amnd(ndx).Val(jv).Exp);
                else
                    if vr==Amnd(ndx).Val(jv).Exp;
                        Flg=1;
                    end
                end
                if Flg~=1;
                    tDh(hcnt)=[];
                    hcnt=hcnt-1;
                    rcnt=rcnt+1;
                end
            end %hcnt
            fprintf('Rejected %d IRs: wrong %s\n',rcnt,Amnd(ndx).Var(jv).Exp)
        end %jv
    else
        tDh=Dh;
    end % isempty(ndx)
    if isempty(tDh);
        tDh=Dh;
    end
    % now remove any classes that are no longer relevant
    vcnt=0;
    while vcnt<length(V); vcnt=vcnt+1;
        cnt=0;
        for jh=1:length(tDh);
            eval(sprintf('if strcmp(tDh(jh).%s,V(vcnt).name); cnt=cnt+1; end',PltPrms{jPrm}));
        end
        if cnt==0;
            V(vcnt)=[];
            vcnt=vcnt-1;
        end
    end
    % plot numbers of IRs per class
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_NoIRs(tDh,PltPrms{jPrm},V); 
    saveas(gcf,sprintf('%s/%s/NoIRs',fNm,PltPrms{jPrm}),ftp);
    % plot peak amplitude
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_PkAmp(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/Amp_m',fNm,PltPrms{jPrm}),ftp);
    % plot subband peak amplitude
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_SbAmp(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/Amp',fNm,PltPrms{jPrm}),ftp);
    % plot spectra
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_Spcs(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/Spc',fNm,PltPrms{jPrm}),ftp);
    % plot RT60
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_RT60(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/RT60',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_mRT60(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/RT60_m',fNm,PltPrms{jPrm}),ftp);
    % plot DRR
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_DRR(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/DRR',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_mDRR(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/DRR_m',fNm,PltPrms{jPrm}),ftp);
    % plot RT60 vs DRR
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_RT60vDRR(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/RT60vDRR',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_mRT60vDRR(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/RT60vDRR_m',fNm,PltPrms{jPrm}),ftp);
    % plot Kurtosis
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_Krt(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/Krt',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_KrtLn(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/Krt_ln',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_Tgs(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/Tgs',fNm,PltPrms{jPrm}),ftp);
    % plot ER spectrum
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_SpcER(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/SpcER',fNm,PltPrms{jPrm}),ftp);
    % plot tail spectrum
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_SpcTl(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('%s/%s/SpcTl',fNm,PltPrms{jPrm}),ftp);
    % make a plot worthy for the paper
    %fcnt=fcnt+1; hfg=figure(fcnt);
    %ps=get(gcf,'position');
    %set(gcf,'position',[ps(1:2) ps(3)*3 ps(4)/2]);
    %hfg=PltIRStts(hfg,tDh,PltPrms{jPrm},V);
    %set(gcf,'PaperPositionMode','auto')
    %saveas(gcf,sprintf('%s/%s/Fg4Ppr',fNm,PltPrms{jPrm}),'epsc');
    %saveas(gcf,sprintf('%s/%s/Fg4Ppr',fNm,PltPrms{jPrm}),'png');
    % Plot Modes
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdRT60(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeRT60',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdmRT60(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/Mode_mRT60',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdRT60Pk(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeRT60Pk',fNm,PltPrms{jPrm}),ftp);

    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdOnPwr(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeOnPwr',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdmOnPwr(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/Mode_mOnPwr',fNm,PltPrms{jPrm}),ftp);

    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdFrqHst(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeHst',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdOPvsRT60(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeOPvsRT60',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdmOvR(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/Mode_mOPvsRT60',fNm,PltPrms{jPrm}),ftp);


    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdSpcs(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeSpc',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdJntDts(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/ModeJntDst',fNm,PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
%    PltIRStts_MdJntDtsSb(tDh,PltPrms{jPrm},V);
%    saveas(gcf,sprintf('%s/%s/ModeJntDstSb',fNm,PltPrms{jPrm}),ftp);

    % Check the recording quality
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_NsFlr(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('%s/%s/NsFlr',fNm,PltPrms{jPrm}),ftp);

%    % plot spectral entropy vs mean(RT60)
%    % Plot spcER
%    % Plot spcGs
%    % plot T-gauss    
%    % Plot spectrum
%    % plot spectral centroid
%    % plot mean Tgs against DRR
    
end
