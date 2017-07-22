function PltIRSpc(H,C)

subplot(2,1,1)
plot(20*log10(abs(H.spc)),H.Spcff/1e3,'b');
hold on
if ~isempty(C)
    for jc=1:length(C);
        plot(20*log10(abs(C(jc).spc)),C(1).Spcff/1e3,'k:');
    end
end
hold off
axis tight
set(gca,'yscale','log')
ylabel('Frequency (kHz)')
xlabel('Power (db)');

Nsc=length(H.Attck);
for js=1:Nsc;
    subplot(3,Nsc,2*Nsc+js)
    plot(20*log10(abs(H.Attck(js).RwSpc)),H.Attck(js).ff/1e3,'b');
    axis tight
    set(gca,'yscale','log')
    set(gca,'ylim',[min(H.ff) max(H.ff)]/1e3);
    title(sprintf('%dms',ceil(H.Attck(js).T*1e3)));
end

