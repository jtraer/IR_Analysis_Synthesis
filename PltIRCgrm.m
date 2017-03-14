function PltIRCgrm(h,H);

Npts=length(h);
[fltbnk,ff,erbff]=make_erb_cos_filters(3*Npts,H.fs,length(H.ff),H.off(1),H.off(end));
Cgrm=generate_subbands([zeros(Npts,1); h; zeros(Npts,1)].',fltbnk);
Cgrm=Cgrm(Npts+[1:Npts],:).'; 
Cgrm=Cgrm(2:end-1,:);

plt=20*log10(abs(Cgrm)); 
pcolor([1:length(H.h)]/H.fs,H.ff/1e3,plt);
axis xy; shading flat
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title([H.Path ': Cochleagram']);
set(gca,'clim',max(max(plt))+[-80 0]);
set(gca,'yscale','log');
colorbar
