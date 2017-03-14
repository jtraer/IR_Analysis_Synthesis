function PltIRts(H)

%*** => compress the time series for plotting
h=sign(H.nh).*abs(H.nh).^(0.6);
%*** Plot
plot([1:length(H.nh)]/H.fs,h);
%*** => plot scale lines
hold on
for jln=1:3
    plot(([2 length(H.nh)]/H.fs),10^(-jln*0.6)*ones(1,2),'k:');
    plot(([2 length(H.nh)]/H.fs),-10^(-jln*0.6)*ones(1,2),'k:');
end
xlabel('Time (s)');
ylabel('Waveform amplitude (compressed)')
set(gca,'xscale','log')
title([H.Path ': Denoised IR'])
