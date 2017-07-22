function PltCAL_Spc(C)
    
lcnt=0;
for jc=1:length(C);
    plot(20*log10(C(jc).Attck(3).RwSpc),C(jc).Attck(3).ff/1e3); hold on
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('%s: ch=%d',C(jc).Name,C(jc).Channel);
end
hold off
set(gca,'yscale','log')
legend(lgnd)
axis tight;
xlabel('Power (dB)');
ylabel('Frequency (kHz)');
OtPth=GtPthStm(GtPthStm(C(1).Path));
title([OtPth ': CAL Spectra']);
