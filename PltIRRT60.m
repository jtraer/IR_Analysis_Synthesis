function PltIRRT60(H,C)

hp=plot(H.RT60,H.ff/1e3);
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on;
hp=plot(H.RT60+H.RT60_std/2,H.ff/1e3,':'); set(hp,'color',pclr)
hp=plot(H.RT60-H.RT60_std/2,H.ff/1e3,':'); set(hp,'color',pclr)
%plot(60./H.Ns.bbt,H.Ns.ff/1e3,'-.');
%plot(60./H.Md.bbt,H.Md.ff/1e3,'d-');
set(gca,'yscale','log','xscale','log');
hold on
if ~isempty(C)
    vndx_t=strcmp({C.Name},'Omni');
    vndx=find(vndx_t==1);
    V=C(vndx);
    V=V(find([V.Channel]==H.Channel));
    plot(V.RT60,H.ff/1e3,'k:');
    plot(H.RT60(H.BdBndsFlg),H.ff(H.BdBndsFlg)/1e3,'k+');
end
hold off;
xlabel('RT60 (s)');
ylabel('Frequency (kHz)')
title([H.Path ': RT60'])
