function PltIRStts_MdRT60(H,PltPrm,Mt,MrkMt,cmpMt)

% scroll through classes
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    Mdf=[];
    MdR=[];
    for jh=1:length(tH);
        P=[tH(jh).Modes.MnPwr];
        [~,ndx]=sort(P,'descend');
        if length(ndx)>3;
            ndx=ndx(1:3);
        end
        Mdf=[Mdf [tH(jh).Modes(ndx).cf]];
        MdR=[MdR [tH(jh).Modes(ndx).RT60]];
    end
    % plot this class
    hp=plot(MdR,Mdf/1e3,sprintf('%s',MrkMt(jj)));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    hold on
end; legend(Mt);
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Mode RT60 (s)')
ylabel('Mode Frequency (kHz)')
title(PltPrm)
