function PltCAL_Spc(D,V)
    
lcnt=0;
for jd=1:length(D);
    plot(D(jd).spc,D(jd).Spcff); hold on
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Direct: ch=%d',D(jd).Channel);
end
for jv=1:length(V)
    plot(V(jv).spc,V(jv).Spcff,'--');
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Omnidirectional: ch=%d',V(jv).Channel);
end
hold off
set(gca,'yscale','log')
legend(lgnd)
xlabel('Power (dB)');
ylabel('Frequency (kHz)');
OtPth=GtPthStm(GtPthStm(D(1).Path));
title([OtPth ': CAL Spectra']);
