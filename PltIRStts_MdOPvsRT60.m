function PltIRStts_MdOPvsRT60(H,PltPrm,Mt,MrkMt,cmpMt)

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
    MdP=[];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.RT60]];
        MdP=[MdP [tH(jh).Modes.OnPwr]];
    end
    % plot this class
    hp=plot(MdR,MdP,sprintf('%s',MrkMt(jj)));
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
    MdP=[];
    MdR=[];
    OPlms=[1e6 0];
    for jh=1:length(tH);
        MdP=[MdP [tH(jh).Modes.OnPwr]];
        MdR=[MdR [tH(jh).Modes.RT60]];
        OPlms(1)=min([OPlms(1) [tH(jh).Modes.OnPwr]]);
        OPlms(2)=max([OPlms(2) [tH(jh).Modes.OnPwr]]);
    end
    PP=linspace(OPlms(1),OPlms(2),5);
    [MdP,srt]=sort(MdP);
    MdR=MdR(srt);
    ndx=find(diff(MdP)==0);
    MdP(ndx+1)=MdP(ndx+1)+1e-3*rand(size(ndx));
    RR=interp1(MdP,MdR,PP);
    % plot
    hp=plot(RR,PP,sprintf('-'));
    set(hp,'color',cmpMt(ceil(length(cmpMt)*jj/length(Mt)),:));
    hold on
end; 
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
i%set(gca,'xscale','log')
xlabel('Mode RT60 (s)')
ylabel('Mode Onset Power (dB)')
title(PltPrm)
