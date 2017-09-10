function PltIRStts_MdOnPwr(Dh,PltPrm,V)
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
    frqlms=[1e6 0];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.OnPwr]];
        frqlms(1)=min([frqlms(1) [tH(jh).Modes.cf]]);
        frqlms(2)=max([frqlms(2) [tH(jh).Modes.cf]]);
    end
    ff=linspace(frqlms(1),frqlms(2),1e1);
    [Mdf,srt]=sort(Mdf);
    MdR=MdR(srt);
    ndx=find(diff(Mdf)==0);
    Mdf(ndx+1)=Mdf(ndx+1)+1e-3*rand(size(ndx));
    RR=interp1(Mdf,MdR,ff);
    % plot this class
    hp=plot(MdR,Mdf/1e3,V(jj).mrk);
    set(hp,'color',V(jj).cmp);
    hold on
    % plot the mean
    hp=plot(mean(MdR),mean(Mdf)/1e3,V(jj).mrk);
    set(hp,'linewidth',3,'markersize',6);
    set(hp,'color',V(jj).cmp);

    %hp=plot(RR,ff/1e3,'-');
    %set(hp,'color',V(jj).cmp);
    %set(hp,'linewidth',1,'markersize',6);
end; 
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
xlabel('Mode Onset power (dB)')
ylabel('Mode Frequency (kHz)')
title(PltPrm)
