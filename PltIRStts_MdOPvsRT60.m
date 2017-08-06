function PltIRStts_MdOPvsRT60(Dh,PltPrm,V)
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
    MdP=[];
    OPlms=[1e6 0];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.RT60]];
        MdP=[MdP [tH(jh).Modes.OnPwr]];
        OPlms(1)=min([OPlms(1) [tH(jh).Modes.OnPwr]]);
        OPlms(2)=max([OPlms(2) [tH(jh).Modes.OnPwr]]);
    end
    PP=linspace(OPlms(1),OPlms(2),5);
    [MdP,srt]=sort(MdP);
    MdR=MdR(srt);
    ndx=find(diff(MdP)==0);
    MdP(ndx+1)=MdP(ndx+1)+1e-3*rand(size(ndx));
    % sometimes OnPwr is NaN (which is not good)
    bndx=unique([find(isnan(MdP)) find(isnan(MdR))]);
    MdP(bndx)=[];
    MdR(bndx)=[];
    RR=interp1(MdP,MdR,PP);
    % plot this class
    hp=plot(MdR,MdP,V(jj).mrk);
    set(hp,'color',V(jj).cmp);
    hold on
    hp=plot(RR,PP,'-');
    set(hp,'color',V(jj).cmp);
end; 
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[(1-0.2*sign(xlm(1)))*xlm(1) (1+0.2*sign(xlm(2)))*xlm(2)]);
%set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Mode RT60 (s)')
ylabel('Mode Onset Power (dB)')
title(PltPrm)
