function PltCAL_RT60(D,V)
    
lcnt=0;
for jd=1:length(D);
    plot(D(jd).RT60,D(jd).ff); hold on
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Direct: ch=%d',D(jd).Channel);
end
for jv=1:length(V)
    plot(V(jv).RT60,V(jv).ff,'--');
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Omnidirectional: ch=%d',V(jv).Channel);
end
hold off
set(gca,'yscale','log')
legend(lgnd)
xlabel('RT60 (s)');
ylabel('Frequency (kHz)');
OtPth=GtPthStm(GtPthStm(D(1).Path));
title([OtPth ': CAL DRR']);
