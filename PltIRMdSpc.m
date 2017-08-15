function PltIRMdSpc(H)

% plot the spectrum
nft=2^ceil(log2(length(H.h)));
spc=fft(H.h,nft);
spc=spc(1:end/2);
Spcff=[1:nft/2]*H.fs/nft;

plot(10*log10(abs(spc)),Spcff/1e3);
hold on;
for jm=1:length(H.Modes)
    [~,ndx]=min(abs(Spcff-H.Modes(jm).cf));
    plot(10*log10([min(abs(spc)) abs(spc(ndx))]),Spcff(ndx)/1e3*ones(1,2),'r-');
end

xlabel('Power (db)');
ylabel('Frequency (kHz)')
set(gca,'yscale','log')
title([H.Path ': Power Spectrum and identified modes'])
