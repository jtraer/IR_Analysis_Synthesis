function hPltStts(Dh,PltPrms);

Txt=0;

% scroll through all the parameters we want to investigate
fcnt=0;
for jPrm=1:length(PltPrms)
    % make a directory
    eval(sprintf('!mkdir -p IRMAudio/%s',PltPrms{jPrm}));

    % define the parameter variables and a colorscheme to plot them
    for jj=1:length(Dh); eval(sprintf('Vr{jj}=Dh(jj).%s;',PltPrms{jPrm})); end
    Vr=unique(Vr);
    MrkVr=repmat(['o','s','*','^','v','+'],[1 100]);
    cmpVr=othercolor('Dark28',32);
    clear V
    for jv=1:length(Vr)
        V(jv).name=Vr{jv};
        V(jv).cmp=cmpVr(ceil(length(cmpVr)*jv/length(Vr)),:)
        V(jv).mrk=MrkVr(jv);
    end

    % plot numbers of IRs per class
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_NoIRs(Dh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/NoIRs',PltPrms{jPrm}),'png');
    % plot peak amplitude
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_PkAmp(Dh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/Amp',PltPrms{jPrm}),'png');
    % plot RT60
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_RT60(Dh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/RT60',PltPrms{jPrm}),'png');
    % plot DRR
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_DRR(Dh,PltPrms{jPrm},V)
    saveas(gcf,sprintf('IRMAudio/%s/DRR',PltPrms{jPrm}),'png');
    % make a plot worthy for the paper
    fcnt=fcnt+1; hfg=figure(fcnt);
    ps=get(gcf,'position');
    set(gcf,'position',[ps(1:2) ps(3)*3 ps(4)/2]);
    hfg=PltIRStts(hfg,Dh,PltPrms{jPrm},V);
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,sprintf('IRMAudio/%s/Fg4Ppr',PltPrms{jPrm}),'epsc');
    %saveas(gcf,sprintf('IRMAudio/%s/Fg4Ppr',PltPrms{jPrm}),'png');
    % Plot Modes
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdRT60(Dh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60',PltPrms{jPrm}),'png');
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdRT60Pk(Dh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60Pk',PltPrms{jPrm}),'png');
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdOnPwr(Dh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeOnPwr',PltPrms{jPrm}),'png');
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdFrqHst(Dh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeHst',PltPrms{jPrm}),'png');
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdOPvsRT60(Dh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeOPvsRT60',PltPrms{jPrm}),'png');
    fcnt=fcnt+1; figure(fcnt);
    PltIRStts_MdSpcs(Dh,PltPrms{jPrm},V);
    saveas(gcf,sprintf('IRMAudio/%s/ModeSpc',PltPrms{jPrm}),'png');

%    % plot spectral entropy vs mean(RT60)
%    % Plot spcER
%    % Plot spcGs
%    % Plot kurtosis
%    % plot T-gauss    
%    % Plot spectrum
%    % plot spectral centroid
%    % plot mean Tgs against DRR
    
end
