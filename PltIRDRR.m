function PltIRDRR(H,C);

hp=plot(H.DRR,H.ff/1e3);
pclr=get(hp,'color');
set(gca,'yscale','log');
hold on;
hp=plot(H.RT60+H.DRR_std/2,H.ff/1e3,':'); set(hp,'color',pclr)
hp=plot(H.RT60-H.DRR_std/2,H.ff/1e3,':'); set(hp,'color',pclr)
%** => plot de-noised recording
hold on
if ~isempty(C)
    plot(C(1).DRR,H.ff/1e3,'k:');
    plot(C(2).DRR,H.ff/1e3,'r:');
    plot(raw_aa,H.ff/1e3,'b:');
    plot(H.DRR(H.BdBndsFlg),H.ff(H.BdBndsFlg)/1e3,'k+');
end
hold off;
xlabel('DRR (s)');
ylabel('Frequency (kHz)')
title([H.Path ': DRR'])
