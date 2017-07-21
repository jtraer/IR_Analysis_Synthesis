function PltIRStts_Spcs(Dh,PltPrm,V)
% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)

% scroll through classes and plot 1 of each class
for jj=1:length(V)
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    for jh=1:length(tH);
        if jj==1; Spcff=H.Spcff; end
        plt=abs(tH(jh).spc);
        plt(find(plt<max(plt)*1e-5))=max(plt)*1e-5;
        plt=20*log10(abs(plt));
        plt=interp1(tH(jh).Spcff,plt,Spcff);
        mplt(:,jh)=plt;
    end
    plt=median(mplt,2);
    err=std(mplt,[],2);
    hp=plot(plt,Spcff/1e3,'-');
    set(hp,'color',V(jj).cmp);
    hold on
end; 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
xlabel('Power (dB)')
ylabel('Frequency (kHz)')
title(PltPrm)