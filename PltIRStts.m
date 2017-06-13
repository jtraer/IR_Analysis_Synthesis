function hfg=PltIRStts(hfg,Dh,PltPrm,V);

figure(hfg);

%* Plot Histogram of amplitudes 
subplot(1,6,1)
% preallocate one data point for each class for the legend
MkLgnd(V)
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    for jh=1:length(tH)
        % x-coordinate is room property
        amp(jh)=tH(jh).MaxAmp;
    end
    % plot
    [hst,xx]=hist(amp,linspace(0,1.2*max(amp),ceil(length(Dh)/length(V)))); 
    hp=plot(xx,hst); hold on;
    set(hp,'color',V(jj).cmp);
end
hold off;
xlabel('Peak Amplitude ');
ylabel('No. of IRs');
xlm=get(gca,'xlim'); set(gca,'xlim',[0 xlm(2)]);
ylm=get(gca,'ylim'); set(gca,'ylim',[0 ylm(2)]);
set(gca,'xscale','log')

%* Plot Kurtosis
subplot(1,6,2)
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    MxL=0;
    for jh=1:length(tH);
        MxL=max([MxL length(tH(jh).krt)]);
    end
    krt=zeros(MxL,1);
    for jh=1:length(tH)
        tmp=tH(jh).krt;
        krt(1:length(tmp))=krt(1:length(tmp))+tmp/length(tH);
    end
    % plot
    hp=plot([1:length(krt)]/H.fs*1e3,krt); hold on
    set(hp,'color',V(jj).cmp);
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
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    plt=zeros(size(tH(1).Attck(3).Spc));
    for jh=1:length(tH)
        tmp=tH(jh).Attck(3).Spc;
        plt=plt+tmp/length(tH);
    end
    % plot
    hp=plot(20*log10(abs(plt)),tH(1).Attck(3).ff/1e3); hold on
    set(hp,'color',V(jj).cmp);
end
hold on;
hold off;
set(gca,'yscale','log')
set(gca,'ylim',[20 20e3]/1e3)
xlabel('Power (dB)');
ylabel('Freq (kHz)');
title(sprintf('ER '));

%%* Plot tail spectrum
subplot(1,6,4)
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    plt=zeros(size(tH(1).Attck(6).Spc));
    for jh=1:length(tH)
        tmp=tH(jh).Attck(6).Spc;
        plt=plt+tmp/length(tH);
    end
    % plot
    hp=plot(20*log10(abs(plt)),tH(1).Attck(6).ff/1e3); hold on
    set(hp,'color',V(jj).cmp);
end
hold on;
hold off;
set(gca,'yscale','log')
set(gca,'ylim',[20 20e3]/1e3)
xlabel('Power (dB)');
ylabel('Freq (kHz)');
title(sprintf('Tail '));

%%* Plot DRR
subplot(1,6,5)
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    plt=zeros(size(tH(1).DRR));
    for jh=1:length(tH)
        tmp=tH(jh).DRR;
        plt=plt+tmp/length(tH);
    end
    % plot
    hp=plot(plt,tH(1).ff/1e3); hold on
    set(hp,'color',V(jj).cmp);
end
set(gca,'yscale','log')
xlabel('DRR (dB)');
ylabel('Freq (kHz)');

%%* Plot RT60
subplot(1,6,6)
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    plt=zeros(size(tH(1).RT60));
    for jh=1:length(tH)
        tmp=tH(jh).RT60;
        plt=plt+tmp/length(tH);
    end
    % plot
    hp=plot(plt,tH(1).ff/1e3); hold on
    set(hp,'color',V(jj).cmp);
end
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('RT60 (s)');
ylabel('Freq (kHz)');
