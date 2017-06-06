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
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.RT60]];
    end
    % plot this class
    hp=plot(MdR,Mdf/1e3,sprintf('%s',MrkMt(jj)));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    hold on
end; legend(Mt);
% now that legend is made we cycle through and plot the individaul data points
for jj=1:length(Mt)
    % compile a substructure of just this class
    jmt=Mt(jj);
    tH=[];
    for jh=1:length(H);
        eval(sprintf('if strcmp(H(jh).%s,Mt(jj)); tH=[tH H(jh)]; end',PltPrm));
    end
    % interpolate
    Mdf=[];
    MdR=[];
    frqlms=[1e6 0];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.RT60]];
        frqlms(1)=min([frqlms(1) [tH(jh).Modes.cf]]);
        frqlms(2)=max([frqlms(2) [tH(jh).Modes.cf]]);
    end
    ff=linspace(frqlms(1),frqlms(2),1e1);
    [Mdf,srt]=sort(Mdf);
    MdR=MdR(srt);
    ndx=find(diff(Mdf)==0);
    Mdf(ndx+1)=Mdf(ndx+1)+1e-3*rand(size(ndx));
    RR=interp1(Mdf,MdR,ff);
    % plot
    hp=plot(RR,ff/1e3,sprintf('-'));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
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
