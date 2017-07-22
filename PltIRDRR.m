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
    vndx_t=strcmp({C.Name},'Omni');
    vndx=find(vndx_t==1);
    V=C(vndx);
    V=V(find([V.Channel]==H.Channel));
    plot(V.DRR,H.ff/1e3,'k:');
    plot(H.DRR(H.BdBndsFlg),H.ff(H.BdBndsFlg)/1e3,'k+');
end
hold off;
xlabel('DRR (s)');
ylabel('Frequency (kHz)')
title([H.Path ': DRR'])
