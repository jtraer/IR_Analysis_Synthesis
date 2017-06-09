function PltCAL_DRR(D,V)
    
lcnt=0;
for jd=1:length(D);
    plot(D(jd).DRR,D(jd).ff); hold on
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Direct: ch=%d',D(jd).Channel);
end
for jv=1:length(V)
    plot(V(jv).DRR,V(jv).ff,'--');
    lcnt=lcnt+1; lgnd{lcnt}=sprintf('Omnidirectional: ch=%d',V(jv).Channel);
end
hold off
set(gca,'yscale','log')
legend(lgnd)
xlabel('DRR (dB)');
ylabel('Frequency (kHz)');
