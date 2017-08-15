function PltIRSpc(H,C)

subplot(2,1,1)
plot(10*log10(abs(H.spc)),H.Spcff/1e3,'b');
hold on
if ~isempty(C)
    for jc=1:length(C);
        plot(10*log10(abs(C(jc).spc)),C(jc).Spcff/1e3,'k:');
    end
end
hold off
axis tight
set(gca,'yscale','log')
set(gca,'ylim',[min(H.ff) max(H.ff)]/1e3);
ylabel('Frequency (kHz)')
xlabel('Power (db)');
title(sprintf('Total IR - IR is %2.2fms start-to-finish',length(H.h)/H.fs*1e3));

Nsc=length(H.Attck);
for js=1:Nsc;
    subplot(3,Nsc,2*Nsc+js)
    plot(10*log10(abs(H.Attck(js).Spc)),H.Attck(js).ff/1e3,'b');
    axis tight
    set(gca,'yscale','log')
    set(gca,'ylim',[min(H.ff) max(H.ff)]/1e3);
    title(sprintf('%dms',ceil(H.Attck(js).T*1e3)));
end

