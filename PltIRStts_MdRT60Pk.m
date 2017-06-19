function PltIRStts_MdRT60(Dh,PltPrm,V)

% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)

% scroll through classes
for jj=1:length(V)
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
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
    hp=plot(MdR,Mdf/1e3,V(jj).mrk);
    set(hp,'color',V(jj).cmp);
    hold on
end; 
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Mode RT60 (s)')
ylabel('Mode Frequency (kHz)')
title(PltPrm)
