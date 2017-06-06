function PltIRStts_MdFrqHst(H,PltPrm,Mt,MrkMt,cmpMt)

% Mode frequencies
subplot(1,3,1);
% scroll through classes
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    Mdf=[];
    MdR=[];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.OnPwr]];
    end
    % compute histogram
    [hst,ff]=histcounts(Mdf);
    hst=[0 hst 0];
    ff=[ff ff(end)+(ff(2)-ff(1))];
    ff(find(ff==0))=20;
    % plot this class
    hp=plot(hst,ff/1e3,sprintf('%s-',MrkMt(jj)));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    hold on
end; legend(Mt);
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
xlabel('No of Modes histogram')
ylabel('Mode Frequency (kHz)')

% Mode spacing frequencies
subplot(1,3,2);
%scroll through classes
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    Mdf=[];
    MdR=[];
    for jh=1:length(tH);
        ttH=tH(jh);
        df=[];
        for jmd=2:length(ttH.Modes);
            df=[df ([ttH.Modes(jmd:end).cf]-ttH.Modes(jmd-1).cf)];
        end
        df=abs(df);
        Mdf=[Mdf df];
        MdR=[MdR [tH(jh).Modes.OnPwr]];
    end
    % compute histogram

    [hst,ff]=histcounts(Mdf);
    hst=[0 hst 0];
    ff=[ff ff(end)+(ff(2)-ff(1))];
    ff(find(ff==0))=20;
    % plot this class
    hp=plot(hst,ff/1e3,sprintf('%s-',MrkMt(jj)));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    hold on
end; legend(Mt);
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
xlabel('No of Modes histogram')
ylabel('Mode Frequency Difference (kHz)')
title(PltPrm)

% Scatter plot of Modes and mode spacings 
subplot(1,3,3);
%scroll through classes
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    Mdf=[];
    Mddf=[];
    for jh=1:length(tH);
        ttH=tH(jh);
        mf=[ttH.Modes.cf];
        mf=sort(mf);
        df=[];
        ff=[];
        for jmd=2:length(mf)-1;
            df=[df (mf(jmd)-mf(jmd-1)) (mf(jmd+1)-mf(jmd))];
            ff=[ff mf(jmd)*ones(1,2)];
        end
        Mddf=[Mdf df];
        Mdf=[Mdf ff];
    end
    % plot this class
    hp=plot(Mddf/1e3,Mdf/1e3,sprintf('%s',MrkMt(jj)));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    hold on
end; legend(Mt);
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
set(gca,'xlim',[0 2]);
%set(gca,'yscale','log')
xlabel('Mode Difference (kHz)')
ylabel('Mode Frequency (kHz)')
