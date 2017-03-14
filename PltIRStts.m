function hfg=PltIRStts(hfg,H,Mt,MrkMt,PltPrms,jPrm,cmpMt);

figure(hfg);

%* Plot Histogram of amplitudes 
subplot(1,6,1)
for jj=1:length(Mt); hp=plot(-1,-1,MrkMt(jj)); hold on
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    set(hp,'linewidth',3,'markersize',10);
end; legend(Mt); 

for jj=1:length(Mt);
    jmt=Mt(jj);
    ttH=[];
    for jh=1:length(H); 
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); ttH=[ttH H(jh)]; end',PltPrms{jPrm}));
    end
    amp=zeros(1,length(ttH));
    for jh=1:length(ttH)
        % x-coordinate is room property
        amp(jh)=ttH(jh).MaxAmp;
    end
    % plot
    [hst,xx]=hist(amp,linspace(0,1,ceil(length(H)/length(Mt)))); 
    hp=plot(xx,hst); hold on;
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
end
hold off;
xlabel('Peak Amplitude ');
ylabel('No. of IRs');
xlm=get(gca,'xlim'); set(gca,'xlim',[0 xlm(2)]);
ylm=get(gca,'ylim'); set(gca,'ylim',[0 ylm(2)]);

%* Plot Kurtosis
subplot(1,6,2)
for jj=1:length(Mt);
    jmt=Mt(jj);
    ttH=[];
    MxL=0;
    for jh=1:length(H); 
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); ttH=[ttH H(jh)]; end',PltPrms{jPrm}));
        MxL=max([MxL length(H(jh).krt)]);
    end
    krt=zeros(MxL,1);
    for jh=1:length(ttH)
        tmp=ttH(jh).krt;
        krt(1:length(tmp))=krt(1:length(tmp))+tmp/length(ttH);
    end
    % plot
    hp=plot([1:length(krt)]/H(1).fs*1e3,krt); hold on
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
end
hold on;
plot([1 length(krt)]/H(1).fs*1e3,3*ones(1,2),'k--');
%plot(H.Tgs*1e3,3,'k+')
%hold off;
set(gca,'xlim',[0 100]);
set(gca,'yscale','log')
xlabel('Time (ms)');
ylabel('Kurtosis');

%%* Plot ER spectrum
subplot(1,6,3)
for jj=1:length(Mt);
    jmt=Mt(jj);
    ttH=[];
    MxL=0;
    for jh=1:length(H); 
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); ttH=[ttH H(jh)]; end',PltPrms{jPrm}));
    end
    plt=zeros(size(ttH(1).Attck(3).Spc));
    for jh=1:length(ttH)
        tmp=ttH(jh).Attck(3).Spc;
        plt=plt+tmp/length(ttH);
    end
    % plot
    hp=plot(plt,ttH(1).Attck(3).ff/1e3); hold on
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
end
hold on;
hold off;
set(gca,'yscale','log')
xlabel('Power (dB)');
ylabel('Freq (kHz)');
title(sprintf('ER '));

%%* Plot tail spectrum
subplot(1,6,4)
for jj=1:length(Mt);
    jmt=Mt(jj);
    ttH=[];
    MxL=0;
    for jh=1:length(H); 
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); ttH=[ttH H(jh)]; end',PltPrms{jPrm}));
    end
    plt=zeros(size(ttH(1).Attck(6).Spc));
    for jh=1:length(ttH)
        tmp=ttH(jh).Attck(6).Spc;
        plt=plt+tmp/length(ttH);
    end
    % plot
    hp=plot(plt,ttH(1).Attck(6).ff/1e3); hold on
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
end
hold on;
hold off;
set(gca,'yscale','log')
xlabel('Power (dB)');
ylabel('Freq (kHz)');
title(sprintf('Tail '));

%%* Plot DRR
subplot(1,6,5)
for jj=1:length(Mt);
    jmt=Mt(jj);
    ttH=[];
    for jh=1:length(H); 
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); ttH=[ttH H(jh)]; end',PltPrms{jPrm}));
    end
    plt=zeros(size(ttH(1).DRR));
    for jh=1:length(ttH)
        tmp=ttH(jh).DRR;
        plt=plt+tmp/length(ttH);
    end
    % plot
    hp=plot(plt,ttH(1).ff/1e3); hold on
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
end
set(gca,'yscale','log')
xlabel('DRR (dB)');
ylabel('Freq (kHz)');

%%* Plot RT60
subplot(1,6,6)
for jj=1:length(Mt);
    jmt=Mt(jj);
    ttH=[];
    for jh=1:length(H); 
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); ttH=[ttH H(jh)]; end',PltPrms{jPrm}));
    end
    plt=zeros(size(ttH(1).RT60));
    for jh=1:length(ttH)
        tmp=ttH(jh).RT60;
        plt=plt+tmp/length(ttH);
    end
    % plot
    hp=plot(plt,ttH(1).ff/1e3); hold on
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
end
set(gca,'yscale','log');
xlabel('RT60 (s)');
ylabel('Freq (kHz)');
