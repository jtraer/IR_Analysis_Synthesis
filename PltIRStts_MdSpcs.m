function PltIRStts_MdSpcs(Dh,PltPrm,V)
% preallocate one data point for each class for the legend
MkLgnd(V)

% scroll through classes and plot 1 of each class
for jj=1:length(V)
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    for jh=1:length(tH);
        plt=abs(tH(jh).spc);
        plt(find(plt<max(plt)*1e-5))=max(plt)*1e-5;
        plt=20*log10(abs(plt));
        plt=plt+80*jj;
        hp=plot(plt,tH(jh).Spcff/1e3,'-');
        set(hp,'color',V(jj).cmp);
        hold on
    end
end; 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
xlabel('Power (dB)')
ylabel('Frequency (kHz)')
title(PltPrm)
