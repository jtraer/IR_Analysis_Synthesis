function H=hPltStts(H,PltPrms);

Txt=0;

% scroll through all the parameters we want to investigate
for jPrm=1:length(PltPrms)
    % make a directory
    eval(sprintf('!mkdir -p IRMAudio/%s',PltPrms{jPrm}));

    % check the Property labels
    %for jj=1:length(H);
    %    eval(sprintf('if size(H(jj).%s,1)>1;',PltPrms{jPrm}));
    %        eval(sprintf('H(jj).%s=H(jj).%s.'';',PltPrms{jPrm},PltPrms{jPrm}));
    %    eval('end');
    %end
    % define colorscheme for materials
    for jj=1:length(H); eval(sprintf('Mt{jj}=H(jj).%s;',PltPrms{jPrm})); end
    Mt=unique(Mt);
    cmpMt=othercolor('Dark28',32);
    MrkMt=repmat(['o','s','*','^','v','+'],[1 100]);
    MrkMt=MrkMt(1:length(Mt));

    % plot numbers of IRs per class
    figure
    for jj=1:length(Mt); hp=plot(-1,-1,MrkMt(jj)); hold on
        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
        set(hp,'linewidth',3,'markersize',10);
    end; legend(Mt); 
    nIRs=zeros(length(Mt),1);
    for jj=1:length(Mt); 
        % x-coordinate is category index
        tH=[];
        for jh=1:length(H);
            eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrms{jPrm}));
        end
        nIRs(jj)=length(tH);
        % plot
        nm=Mt(jj);
        hp=plot(jj,nIRs(jj),MrkMt(find(strcmp(nm,Mt)))); hold on
        set(hp,'linewidth',3,'markersize',6);
        set(hp,'color',cmpMt(ceil(length(cmpMt)*find(strcmp(nm,Mt))/length(Mt)),:));
    end; hold off 
    axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
    set(gca,'xlim',[0 1.2*xlm(2)]);
    set(gca,'ylim',[0 1.2*ylm(2)]);
    xlabel('IR category')
    ylabel('No. of IRs')
    title(PltPrms{jPrm})
    %saveas(gcf,sprintf('IRMAudio/%s/NoIRs',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/NoIRs',PltPrms{jPrm}),'png');

%    % plot spectral entropy vs mean(RT60)
%    tH=H; tMt=Mt; tMrkMt=MrkMt; clear Err
%    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
%    %for jrm=1:length(RmvLst)
%    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
%    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %end
%    fprintf('we have %d/%d IRs left\n',length(tH),length(H));
%    figure
%    for jj=1:length(tMt); hp=plot(-1,-1,tMrkMt(jj)); hold on
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        set(hp,'linewidth',3,'markersize',10);
%    end; legend(tMt); 
%    for jj=1:length(tH); 
%        % x-coordinate is median decay Rate
%        pltx(jj)=median(tH(jj).RT60);
%        % y-coordinate is "Entropy" of spectrum
%        [Spc,ff]=pwelch(tH(jj).nh,4096,2048,4096,tH(jj).fs);
%        Spc=Spc/sum(Spc); 
%        SpcH=-sum(Spc.*log2(Spc))/length(Spc);
%        plty(jj)=SpcH;
%        % plot
%        eval(sprintf('nm=tH(jj).%s;',PltPrms{jPrm}));
%        hp=plot(pltx(jj),plty(jj),MrkMt(find(strcmp(nm,tMt)))); hold on
%        set(hp,'linewidth',3,'markersize',6);
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*find(strcmp(nm,tMt))/length(tMt)),:));
%        if Txt==1;
%            ht=text(1.001*pltx(jj),1.001*plty(jj),1.1,sprintf('%d',jj));
%            set(ht,'fontsize',10)
%            set(ht,'color',cmpMt(ceil(length(cmpMt)*find(strcmp(nm,tMt))/length(tMt)),:));
%        end
%    end; hold off 
%    axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%    set(gca,'xlim',median(pltx)+2*std(pltx)*[-1 1]);
%    set(gca,'ylim',median(plty)+2*std(plty)*[-1 1]);
%    xlabel('Median RT60 (s)')
%    ylabel('Spectral "Entropy"')
%    title(PltPrms{jPrm})
%    saveas(gcf,sprintf('IRMAudio/%s/H_vs_T',PltPrms{jPrm}),'epsc');
%    saveas(gcf,sprintf('IRMAudio/%s/H_vs_T',PltPrms{jPrm}),'png');
    
    % Plot amplitude
    tH=H; tMt=Mt; tMrkMt=MrkMt; clear Err
    figure
    tmp=zeros(1,length(tH)); 
    for jj=1:length(tH); tmp(jj)=length(tH(jj).ff); end; 
    Nbnds=mode(tmp);
    ndx=find(tmp==Nbnds);
    ff=tH(ndx(1)).ff;
    for jj=1:length(tMt)
        jmt=tMt(jj);
        ttH=[];
        for jh=1:length(tH);
            eval(sprintf('if strcmp(tH(jh).%s,Mt(jj)); ttH=[ttH tH(jh)]; end',PltPrms{jPrm}));
        end
        tmp=zeros(1,length(ttH));
        for jj2=1:length(ttH);
            tmp(jj2)=ttH(1,jj2).MaxAmp;
        end
        Err(jj,:)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
        hp=plot(jj,median(tmp,2),sprintf('%s-',tMrkMt(jj))); hold on
        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
    end; legend(tMt);
    for jj=1:length(tMt)
        hp2=plot(jj*ones(1,2),Err(jj,:)); hold on
        set(hp2,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
    end
    hold off
    %set(gca,'yscale','log')
    %set(gca,'ylim',[50 15e3])
    set(gca,'xlim',[0 (jj+1)]);
    xlabel('Class Index')
    ylabel('Amplitude (re Max Measured)')
    title(PltPrms{jPrm})
    saveas(gcf,sprintf('IRMAudio/%s/Amp',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/Amp',PltPrms{jPrm}),'png');
        
    % Plot RT60
    tH=H; tMt=Mt; tMrkMt=MrkMt; clear Err
    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
    %for jrm=1:length(RmvLst)
    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
    %end
    %H=H(strcmp({H.Shape},'Cave')); Mt={'Cave'};
    %ndx=[]; for jj=1:length(tH); if length(tH(jj).ff)<20; ndx=[ndx; jj]; end; end
    %tH(ndx)=[]; 
    figure
    tmp=zeros(1,length(tH)); 
    for jj=1:length(tH); tmp(jj)=length(tH(jj).ff); end; 
    Nbnds=mode(tmp);
    ndx=find(tmp==Nbnds);
    ff=tH(ndx(1)).ff;
    for jj=1:length(tMt)
        jmt=tMt(jj);
        ttH=[];
        for jh=1:length(tH);
            eval(sprintf('if strcmp(tH(jh).%s,Mt(jj)); ttH=[ttH tH(jh)]; end',PltPrms{jPrm}));
        end
        tmp=zeros(Nbnds,length(ttH));
        for jj2=1:length(ttH);
            tmp(:,jj2)=interp1(ttH(jj2).ff,ttH(jj2).RT60,ff);
        end
        Err(:,:,jj)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
        hp=plot(median(tmp,2),ff,sprintf('%s-',tMrkMt(jj))); hold on
        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
    end; legend(tMt);
    for jj=1:length(tMt)
        for jf=[1:floor(length(ff)/5):length(ff)];
            hp=plot(Err(jf,:,jj),ff(jf)*ones(1,2));
            set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
        end
    end
    hold off
    set(gca,'yscale','log')
    set(gca,'ylim',[50 15e3])
    set(gca,'xlim',[0 2]);
    xlabel('RT60 (s)')
    ylabel('Frequency (Hz)')
    title(PltPrms{jPrm})
    saveas(gcf,sprintf('IRMAudio/%s/RT60',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/RT60',PltPrms{jPrm}),'png');

    % Plot Modes
    figure
    PltIRStts_MdRT60(H,PltPrms{jPrm},Mt,MrkMt,cmpMt);
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60',PltPrms{jPrm}),'png');
    PltIRStts_MdRT60Pk(H,PltPrms{jPrm},Mt,MrkMt,cmpMt);
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60Pk',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/ModeRT60Pk',PltPrms{jPrm}),'png');
    figure
    figure
    PltIRStts_MdOnPwr(H,PltPrms{jPrm},Mt,MrkMt,cmpMt);
    saveas(gcf,sprintf('IRMAudio/%s/ModeOnPwr',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/ModeOnPwr',PltPrms{jPrm}),'png');
    figure
    PltIRStts_MdFrqHst(H,PltPrms{jPrm},Mt,MrkMt,cmpMt);
    saveas(gcf,sprintf('IRMAudio/%s/ModeHst',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/ModeHst',PltPrms{jPrm}),'png');
    figure
    PltIRStts_MdOPvsRT60(H,PltPrms{jPrm},Mt,MrkMt,cmpMt);
    saveas(gcf,sprintf('IRMAudio/%s/ModeOPvsRT60',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/ModeOPvsRT60',PltPrms{jPrm}),'png');
    figure
    PltIRStts_MdSpcs(H,PltPrms{jPrm},Mt,MrkMt,cmpMt);
    saveas(gcf,sprintf('IRMAudio/%s/ModeSpc',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/ModeSpc',PltPrms{jPrm}),'png');
    
    % Plot DRR
    tH=H; tMt=Mt; tMrkMt=MrkMt;
    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
    %for jrm=1:length(RmvLst)
    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
    %end
    %H=H(strcmp({H.Shape},'Cave')); Mt={'Cave'};
    %ndx=[]; for jj=1:length(tH); if length(tH(jj).ff)<20; ndx=[ndx; jj]; end; end
    %tH(ndx)=[]; 
    figure
    tmp=zeros(1,length(tH)); 
    for jj=1:length(tH); tmp(jj)=length(tH(jj).ff); end; 
    Nbnds=mode(tmp);
    ndx=find(tmp==Nbnds);
    ff=tH(ndx(1)).ff;
    for jj=1:length(tMt)
        jmt=tMt(jj);
        ttH=[];
        for jh=1:length(tH);
            eval(sprintf('if strcmp(tH(jh).%s,Mt(jj)); ttH=[ttH tH(jh)]; end',PltPrms{jPrm}));
        end
        tmp=zeros(Nbnds,length(ttH));
        for jj2=1:length(ttH);
            tmp(:,jj2)=interp1(ttH(jj2).ff,-ttH(jj2).DRR,ff);
        end
        Err(:,:,jj)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
        hp=plot(median(tmp,2),ff,sprintf('%s-',tMrkMt(jj))); hold on
        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
    end; legend(tMt);
    for jj=1:length(tMt)
        for jf=[1:floor(length(ff)/5):length(ff)];
            hp=plot(Err(jf,:,jj),ff(jf)*ones(1,2));
            set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
        end
    end
    hold off
    axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
    set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
    set(gca,'yscale','log')
    set(gca,'ylim',[50 15e3])
    xlabel('DRR (dB)')
    ylabel('Frequency (Hz)')
    title(PltPrms{jPrm})
    saveas(gcf,sprintf('IRMAudio/%s/DRR',PltPrms{jPrm}),'epsc');
    saveas(gcf,sprintf('IRMAudio/%s/DRR',PltPrms{jPrm}),'png');

%    
%    % Plot spcER
%    tH=H; tMt=Mt; tMrkMt=MrkMt;
%    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
%    %for jrm=1:length(RmvLst)
%    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
%    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %end
%    %H=H(strcmp({H.Shape},'Cave')); Mt={'Cave'};
%    %ndx=[]; for jj=1:length(tH); if length(tH(jj).ff)<20; ndx=[ndx; jj]; end; end
%    %tH(ndx)=[]; 
%    figure
%    tmp=zeros(1,length(tH)); 
%    for jj=1:length(tH); tmp(jj)=length(tH(jj).ff); end; 
%    Nbnds=mode(tmp);
%    ndx=find(tmp==Nbnds);
%    ff=tH(ndx(1)).ff;
%    for jj=1:length(tMt)
%        jmt=tMt(jj);
%        eval(sprintf('ttH=tH(find(strcmp({tH.%s},jmt)));',PltPrms{jPrm}));
%        tmp=zeros(Nbnds,length(ttH));
%        for jj2=1:length(ttH);
%            tmp(:,jj2)=interp1(ttH(jj2).ff,ttH(jj2).spcER,ff);
%        end
%        Err(:,:,jj)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
%        hp=plot(median(tmp,2),ff,sprintf('%s-',tMrkMt(jj))); hold on
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%    end; legend(tMt);
%    for jj=1:length(tMt)
%        for jf=[1:floor(length(ff)/5):length(ff)];
%            hp=plot(Err(jf,:,jj),ff(jf)*ones(1,2));
%            set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        end
%    end
%    hold off
%    axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%    set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%    set(gca,'yscale','log')
%    set(gca,'ylim',[50 15e3])
%    xlabel('Power (dB)')
%    ylabel('Frequency (Hz)')
%    title(PltPrms{jPrm})
%    saveas(gcf,sprintf('IRMAudio/%s/spcER',PltPrms{jPrm}),'epsc');
%    saveas(gcf,sprintf('IRMAudio/%s/spcER',PltPrms{jPrm}),'png');
%
%    % Plot spcGs
%    tH=H; tMt=Mt; tMrkMt=MrkMt;
%    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
%    %for jrm=1:length(RmvLst)
%    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
%    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %end
%    %H=H(strcmp({H.Shape},'Cave')); Mt={'Cave'};
%    %ndx=[]; for jj=1:length(tH); if length(tH(jj).ff)<20; ndx=[ndx; jj]; end; end
%    %tH(ndx)=[]; 
%    figure
%    tmp=zeros(1,length(tH)); 
%    for jj=1:length(tH); tmp(jj)=length(tH(jj).ff); end; 
%    Nbnds=mode(tmp);
%    ndx=find(tmp==Nbnds);
%    ff=tH(ndx(1)).ff;
%    for jj=1:length(tMt)
%        jmt=tMt(jj);
%        eval(sprintf('ttH=tH(find(strcmp({tH.%s},jmt)));',PltPrms{jPrm}));
%        tmp=zeros(Nbnds,length(ttH));
%        for jj2=1:length(ttH);
%            tmp(:,jj2)=interp1(ttH(jj2).ff,ttH(jj2).spcGs,ff);
%        end
%        Err(:,:,jj)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
%        hp=plot(median(tmp,2),ff,sprintf('%s-',tMrkMt(jj))); hold on
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%    end; legend(tMt);
%    for jj=1:length(tMt)
%        for jf=[1:floor(length(ff)/5):length(ff)];
%            hp=plot(Err(jf,:,jj),ff(jf)*ones(1,2));
%            set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        end
%    end
%    hold off
%    axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%    set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%    set(gca,'yscale','log')
%    set(gca,'ylim',[50 15e3])
%    xlabel('Power (dB)')
%    ylabel('Frequency (Hz)')
%    title(PltPrms{jPrm})
%    saveas(gcf,sprintf('IRMAudio/%s/spcGs',PltPrms{jPrm}),'epsc');
%    saveas(gcf,sprintf('IRMAudio/%s/spcGs',PltPrms{jPrm}),'png');
%
%    % Plot kurtosis
%    tH=H; tMt=Mt; tMrkMt=MrkMt;
%    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
%    %for jrm=1:length(RmvLst)
%    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
%    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %end
%    %H=H(strcmp({H.Shape},'Cave')); Mt={'Cave'};
%    %ndx=[]; for jj=1:length(tH); if length(tH(jj).ff)<20; ndx=[ndx; jj]; end; end
%    %tH(ndx)=[]; 
%    figure
%    tmp=zeros(1,length(tH)); clear Err
%    for jj=1:length(tH); tmp(jj)=length(tH(jj).tt); end; 
%    Nt=max(tmp);
%    ndx=find(tmp==Nt);
%    tt=tH(ndx(1)).tt;
%    for jj=1:length(tMt)
%        jmt=tMt(jj);
%        eval(sprintf('ttH=tH(find(strcmp({tH.%s},jmt)));',PltPrms{jPrm}));
%        tmp=zeros(Nt,length(ttH));
%        for jj2=1:length(ttH);
%            tmp(:,jj2)=interp1(ttH(jj2).tt,ttH(jj2).krt,tt);
%        end
%        Err(:,:,jj)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
%        hp=plot(tt,median(tmp,2),sprintf('%s-',tMrkMt(jj))); hold on
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        % get the average value of Tgauss
%        tmp2(jj)=median([ttH.Tgs]);
%        Err2(jj,:)=median([ttH.Tgs])+std([ttH.Tgs],[],2)*[-1 1]/2;
%    end; legend(tMt);
%    for jj=1:length(tMt)
%        for jf=[1:floor(length(tt)/5):length(tt)];
%            hp=plot(tt(jf)*ones(1,2),Err(jf,:,jj));
%            set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        end
%    end
%    % plot T-gauss    
%    for jj=1:length(tMt)
%        hp=plot(tmp2(jj)*ones(1,2),[2 1e1],sprintf('%s-',tMrkMt(jj)));
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        hp=plot(Err2(jj,1)*ones(1,2),[2 1e1],sprintf('%s--',tMrkMt(jj)));
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        hp=plot(Err2(jj,2)*ones(1,2),[2 1e1],sprintf('%s--',tMrkMt(jj)));
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%    end
%    % plot the limits of Gaussianity
%    plot([0 1],3*ones(1,2),'k:')
%    hold off
%    set(gca,'yscale','log')
%    set(gca,'ylim',[1 100])
%    set(gca,'xlim',[0 0.05]);
%    xlabel('time (s)')
%    ylabel('Kurtosis')
%    title(PltPrms{jPrm})
%    saveas(gcf,sprintf('IRMAudio/%s/Krt',PltPrms{jPrm}),'epsc');
%    saveas(gcf,sprintf('IRMAudio/%s/Krt',PltPrms{jPrm}),'png');
%
%    % Plot spectrum
%    tH=H; tMt=Mt; tMrkMt=MrkMt;
%    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
%    %for jrm=1:length(RmvLst)
%    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
%    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %end
%    %H=H(strcmp({H.Shape},'Cave')); Mt={'Cave'};
%    %ndx=[]; for jj=1:length(tH); if length(tH(jj).ff)<20; ndx=[ndx; jj]; end; end
%    %tH(ndx)=[]; 
%    figure
%    tmp=zeros(1,length(tH)); clear Err
%    for jj=1:length(tH); tmp(jj)=length(tH(jj).tt); end; 
%    Nt=max(tmp);
%    ndx=find(tmp==Nt);
%    ff=logspace(1,4,256);
%    for jj=1:length(tMt)
%        jmt=tMt(jj);
%        eval(sprintf('ttH=tH(find(strcmp({tH.%s},jmt)));',PltPrms{jPrm}));
%        tmp=zeros(length(ff),length(ttH));
%        for jj2=1:length(ttH);
%            tmp(:,jj2)=10*log10(pwelch(ttH(jj2).nh,1024,512,ff,ttH(jj2).fs));
%        end
%        Err(:,:,jj)=median(tmp,2)*ones(1,2)+std(tmp,[],2)*[-1 1]/2;
%        %hp=plot(median(tmp,2),ff,sprintf('%s-',tMrkMt(jj))); hold on
%        hp=plot(median(tmp,2),ff); hold on
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        % get the average spectral centroid
%        spc=median(tmp,2);
%        tmp2(jj)=sum((spc).*[1:length(ff)].')/sum((spc));
%    end; legend(tMt);
%    for jj=1:length(tMt)
%        for jf=[1:floor(length(ff)/5):length(ff)];
%            hp=plot(Err(jf,:,jj),ff(jf)*ones(1,2),sprintf('%s-',tMrkMt(jj)));
%            set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        end
%    end
%    % plot spectral centroid
%    for jj=1:length(tMt)
%        hp=plot([min(spc) max(spc)],ff(ceil(tmp2(jj)))*ones(1,2),sprintf('%s:',tMrkMt(jj)));
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%    end
%    hold off
%    set(gca,'yscale','log')
%    set(gca,'ylim',[10 1e4])
%    set(gca,'xlim',[-120 -60]);
%    xlabel('Power (dB)')
%    ylabel('Frequency (Hz)')
%    title(PltPrms{jPrm})
%    saveas(gcf,sprintf('IRMAudio/%s/Spc',PltPrms{jPrm}),'epsc');
%    saveas(gcf,sprintf('IRMAudio/%s/Spc',PltPrms{jPrm}),'png');
%
%    % plot mean Tgs against DRR
%    tH=H; tMt=Mt; tMrkMt=MrkMt; clear Err
%    %RmvLst={'??';'Crevice';'Open';'OpenSpace';'Plank';'Handle';'Other';'Rocky';'RockyBowl';'Canyon';'Alcove';'Forest'};
%    %for jrm=1:length(RmvLst)
%    %    tH(find(strcmp({tH.Shape},RmvLst{jrm})))=[];
%    %    tMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %    tMrkMt(find(strcmp(tMt,RmvLst{jrm})))=[];
%    %end
%    fprintf('we have %d/%d IRs left\n',length(tH),length(H));
%    figure
%    for jj=1:length(tMt); hp=plot(-1,-1,tMrkMt(jj)); hold on
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(tMt)),:));
%        set(hp,'linewidth',3,'markersize',10);
%    end; legend(tMt); 
%    for jj=1:length(tH); 
%        % x-coordinate is median DRR
%        pltx(jj)=median(tH(jj).DRR);
%        % y-coordinate is Tgs
%        plty(jj)=tH(jj).Tgs;
%        % plot
%        eval(sprintf('nm=tH(jj).%s;',PltPrms{jPrm}));
%        hp=plot(pltx(jj),plty(jj),MrkMt(find(strcmp(nm,tMt)))); hold on
%        set(hp,'linewidth',3,'markersize',6);
%        set(hp,'color',cmpMt(ceil(length(cmpMt)*find(strcmp(nm,tMt))/length(tMt)),:));
%        if Txt==1;
%            ht=text(1.001*pltx(jj),1.001*plty(jj),1.1,sprintf('%d',jj));
%            set(ht,'fontsize',10)
%            set(ht,'color',cmpMt(ceil(length(cmpMt)*find(strcmp(nm,tMt))/length(tMt)),:));
%        end
%    end; hold off 
%    axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%    set(gca,'xlim',[1.2*xlm(1) 1.2*xlm(2)]);
%    set(gca,'ylim',[1.2*ylm(1) 1.2*ylm(2)]);
%    xlabel('Median DRR (s)')
%    ylabel('T_{Gauss}')
%    title(PltPrms{jPrm})
%    saveas(gcf,sprintf('IRMAudio/%s/Tgs_vs_mDRR',PltPrms{jPrm}),'epsc');
%    saveas(gcf,sprintf('IRMAudio/%s/Tgs_vs_mDRR',PltPrms{jPrm}),'png');

    % make a plot worthy for the paper
    tH=H; tMt=Mt; tMrkMt=MrkMt; clear Err
    hfg=figure;
    ps=get(gcf,'position');
    set(gcf,'position',[ps(1:2) ps(3)*3 ps(4)/2]);
    hfg=PltIRStts(hfg,tH,tMt,tMrkMt,PltPrms,jPrm,cmpMt);
    set(gcf,'PaperPositionMode','auto')
    saveas(gcf,sprintf('IRMAudio/%s/Fg4Ppr',PltPrms{jPrm}),'epsc');
    %saveas(gcf,sprintf('IRMAudio/%s/Fg4Ppr',PltPrms{jPrm}),'png');
end
