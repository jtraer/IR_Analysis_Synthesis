function PltIRStts_PkAmp(Dh,PltPrm,V);

% preallocate one data point for each class for the legend
MkLgnd(V)
    
plt=zeros(length(V),1);
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissaa
    xx(jj)=jj;
    plt(jj)=mean(20*log10([tH.MaxAmp]));
    err(jj)=std(20*log10([tH.MaxAmp]));
    % plot
    hp=errorbar(xx(jj),plt(jj),err(jj),V(jj).mrk); hold on
    set(hp,'linewidth',3,'markersize',6);
    set(hp,'color',V(jj).cmp);
    hp=text(xx(jj)+0.2,plt(jj),1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
end; hold off 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[0 1.2*xlm(2)]);
set(gca,'ylim',[0 1.2*ylm(2)]);
xlabel('IR category')
ylabel('Peak Amplitude (dB)')
title(PltPrm)
