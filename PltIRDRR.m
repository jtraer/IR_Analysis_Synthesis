function PltIRDRR(H,C,raw_aa);

hp=plot(H.DRR,H.ff/1e3);
pclr=get(hp,'color');
set(gca,'yscale','log');
hold on;
hp=plot(H.RT60+H.DRR_std/2,H.ff/1e3,':'); set(hp,'color',pclr)
hp=plot(H.RT60-H.DRR_std/2,H.ff/1e3,':'); set(hp,'color',pclr)
%** => plot de-noised recording
hold on
if ~isempty(C)
    D=C.Direct;
    V=C.Omni;
    tD=D(find([D.Channel]==H.Channel));
    tV=V(find([V.Channel]==H.Channel));
    plot(tD.DRR,H.ff/1e3,'k:');
    plot(tV.DRR,H.ff/1e3,'r:');
    plot(raw_aa,H.ff/1e3,'b:');
    plot(H.DRR(H.BdBndsFlg),H.ff(H.BdBndsFlg)/1e3,'k+');
end
hold off;
xlabel('DRR (s)');
ylabel('Frequency (kHz)')
title([H.Path ': DRR'])
