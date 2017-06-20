function PltCAL_Spc(D,V)
    
lcnt=0;
for jd=1:length(D);
    plot(20*log10(D(jd).Attck(3).RwSpc),D(jd).Attck(3).ff/1e3); hold on
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Direct: ch=%d',D(jd).Channel);
end
for jv=1:length(V)
    plot(20*log10(V(jv).Attck(3).RwSpc),V(jv).Attck(3).ff/1e3,'--');
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Omnidirectional: ch=%d',V(jv).Channel);
end
hold off
set(gca,'yscale','log')
legend(lgnd)
axis tight;
xlabel('Power (dB)');
ylabel('Frequency (kHz)');
OtPth=GtPthStm(GtPthStm(D(1).Path));
title([OtPth ': CAL Spectra']);
