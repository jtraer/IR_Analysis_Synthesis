function PltIRStts_MdFrqHst(Dh,PltPrm,V)
% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)

% Mode frequencies
%subplot(1,3,1);
% scroll through classes
for jj=1:length(V)
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    Mdf=[];
    MdR=[];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.OnPwr]];
    end
    % remove any high frequency peaks
    ndx=find(Mdf>20e3);
    Mdf(ndx)=[];
    MdR(ndx)=[];
    % compute histogram
    [hst,ff]=histcounts(Mdf,10);
    hst=[0 hst 0];
    ff=[ff ff(end)+(ff(2)-ff(1))];
    ff(find(ff==0))=20;
    % plot this class
    hp=plot(hst,ff/1e3,sprintf('%s-',V(jj).mrk));
    set(hp,'color',V(jj).cmp);
    hold on
end; 
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
set(gca,'ylim',[0.02 20]);
%set(gca,'yscale','log')
xlabel('No of Modes histogram')
ylabel('Mode Frequency (kHz)')

% Mode spacing frequencies
%subplot(1,3,2);
%%scroll through classes
%for jj=1:length(V)
%    % collate all IRs that have this particular label
%    tH=[]; 
%    for jh=1:length(Dh);
%        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
%    end
%    % specify the ordinates and abscissa
%    Mdf=[];
%    MdR=[];
%    for jh=1:length(tH);
%        ttH=tH(jh);
%        df=[];
%        for jmd=2:length(ttH.Modes);
%            df=[df ([ttH.Modes(jmd:end).cf]-ttH.Modes(jmd-1).cf)];
%        end
%        df=abs(df);
%        Mdf=[Mdf df];
%        MdR=[MdR [tH(jh).Modes.OnPwr]];
%    end
%    % compute histogram
%    [hst,ff]=histcounts(Mdf);
%    hst=[0 hst 0];
%    ff=[ff ff(end)+(ff(2)-ff(1))];
%    ff(find(ff==0))=20;
%    % plot this class
%    hp=plot(hst,ff/1e3,V(jj).mrk);
%    set(hp,'color',V(jj).cmp);
%    hold on
%end; 
%hold off
%axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%%set(gca,'yscale','log')
%xlabel('No of Modes histogram')
%ylabel('Mode Frequency Difference (kHz)')
%title(PltPrm)
%
%% Scatter plot of Modes and mode spacings 
%subplot(1,3,3);
%%scroll through classes
%for jj=1:length(V)
%    % collate all IRs that have this particular label
%    tH=[]; 
%    for jh=1:length(Dh);
%        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
%    end
%    % specify the ordinates and abscissa
%    Mdf=[];
%    Mddf=[];
%    for jh=1:length(tH);
%        ttH=tH(jh);
%        mf=[ttH.Modes.cf];
%        mf=sort(mf);
%        df=[];
%        ff=[];
%        for jmd=2:length(mf)-1;
%            df=[df (mf(jmd)-mf(jmd-1)) (mf(jmd+1)-mf(jmd))];
%            ff=[ff mf(jmd)*ones(1,2)];
%        end
%        Mddf=[Mdf df];
%        Mdf=[Mdf ff];
%    end
%    % plot this class
%    hp=plot(Mddf/1e3,Mdf/1e3,V(jj).mrk);
%    set(hp,'color',V(jj).cmp);
%    hold on
%end; 
%hold off
%axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'xlim',[0 2]);
%%set(gca,'yscale','log')
%xlabel('Mode Difference (kHz)')
%ylabel('Mode Frequency (kHz)')
