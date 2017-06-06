function PltIRStts_MdSpcs(H,PltPrm,Mt,MrkMt,cmpMt)

% scroll through classes and plot 1 of each class
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    % plot this class
    for jh=1;
        plt=abs(tH(jh).spc);
        plt(find(plt<max(plt)*1e-5))=max(plt)*1e-5;
        plt=20*log10(abs(plt));
        plt=plt+80*jj;
        hp=plot(plt,tH(jh).Spcff/1e3,sprintf('-',MrkMt(jj)));
        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
        hold on
    end
end; legend(Mt);
% now that legend is made we cycle through and plot all the Spectra
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    % plot this class
    for jh=1:length(tH);
        plt=abs(tH(jh).spc);
        plt(find(plt<max(plt)*1e-5))=max(plt)*1e-5;
        plt=20*log10(abs(plt));
        plt=plt+80*jj;
        hp=plot(plt,tH(jh).Spcff/1e3,sprintf('-',MrkMt(jj)));
        set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
        hold on
    end
end; 
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
xlabel('Power (dB)')
ylabel('Frequency (kHz)')
title(PltPrm)
