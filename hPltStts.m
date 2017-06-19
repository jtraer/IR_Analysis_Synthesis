function hPltStts(Dh,PltPrms,Amnd);

ftp='png';
ftp='epsc';
Txt=0;

% scroll through all the parameters we want to investigate
fcnt=0;
for jPrm=1:length(PltPrms)
    % make a directory
    eval(sprintf('!mkdir -p IRMAudio/%s',PltPrms{jPrm}));

    % define the parameter variables and a colorscheme to plot them
    for jj=1:length(Dh); eval(sprintf('Vr{jj}=Dh(jj).%s;',PltPrms{jPrm})); end
    Vr=unique(Vr);
    Vr(find(strcmp('NA',Vr)))=[]; 
    MrkVr=repmat(['o','s','*','^','v','+'],[1 100]);
    cmpVr=othercolor('Dark28',32);
    clear V
    for jv=1:length(Vr)
        V(jv).name=Vr{jv};
        V(jv).cmp=cmpVr(ceil(length(cmpVr)*jv/length(Vr)),:)
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
    for jv=1:length(V);
        cnt=0;
        for jh=1:length(Dh);
            eval(sprintf('if strcmp(Dh(jh).%s,V(jv).name); cnt=cnt+1; end',PltPrms{jPrm}));
        end
        if cnt==0;
            V(jv)=[];
        end
    end
    % plot numbers of IRs per class
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_NoIRs(tDh,PltPrms{jPrm},V); 
    saveas(gcf,sprintf('IRMAudio/%s/NoIRs',PltPrms{jPrm}),ftp);
    % plot peak amplitude
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_PkAmp(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/Amp',PltPrms{jPrm}),ftp);
    % plot RT60
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_RT60(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/RT60',PltPrms{jPrm}),ftp);
    % plot DRR
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_DRR(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/DRR',PltPrms{jPrm}),ftp);
    % plot Kurtosis
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_Krt(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/Krt',PltPrms{jPrm}),ftp);
    % plot ER spectrum
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_SpcER(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/SpcER',PltPrms{jPrm}),ftp);
    % plot tail spectrum
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_SpcTl(tDh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/SpcTl',PltPrms{jPrm}),ftp);
    % make a plot worthy for the paper
    fcnt=fcnt+1; hfg=figure(fcnt);
    ps=get(gcf,'position');
    set(gcf,'position',[ps(1:2) ps(3)*3 ps(4)/2]);
    hfg=PltIRStts(hfg,tDh,PltPrms{jPrm},V);
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,sprintf('IRMAudio/%s/Fg4Ppr',PltPrms{jPrm}),'epsc');
    %saveas(gcf,sprintf('IRMAudio/%s/Fg4Ppr',PltPrms{jPrm}),'png');
    % Plot Modes
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdRT60(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60',PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdRT60Pk(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60Pk',PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdOnPwr(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeOnPwr',PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdFrqHst(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeHst',PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdOPvsRT60(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeOPvsRT60',PltPrms{jPrm}),ftp);
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdSpcs(tDh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeSpc',PltPrms{jPrm}),ftp);

%    % plot spectral entropy vs mean(RT60)
%    % Plot spcER
%    % Plot spcGs
%    % Plot kurtosis
%    % plot T-gauss    
%    % Plot spectrum
%    % plot spectral centroid
%    % plot mean Tgs against DRR
    
end
