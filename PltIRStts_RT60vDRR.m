function PltIRStts_RT60_vDRR(Dh,PltPrm,V);

% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)
    
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    ff=H.ff/1e3; mplt=zeros(length(ff),length(tH));
    for jh=1:length(tH);
        mplt(:,jh)=tH(jh).RT60;
        mplt2(:,jh)=tH(jh).DRR;
    end
    plt=median(mplt,2);
    err=std(mplt,[],2);
    plt2=median(mplt2,2);
    err2=std(mplt2,[],2);
    % plot
    hp=plot(plt,plt2,[V(jj).mrk '']); hold on
    set(hp,'linewidth',1,'markersize',6);
    set(hp,'color',V(jj).cmp);
    [mx,mxndx]=max(plt);
    hp=text(mx+0.1,ff(mxndx),1.001,V(jj).name);
    set(hp,'color',V(jj).cmp);
    for jp=1:length(plt)
        hp=plot(plt(jp)+err*[-1 1],plt2(jp)*ones(1,2),'-'); hold on
        set(hp,'color',V(jj).cmp);
        hp=plot(plt(jp)*ones(1,2),plt2(jp)+err2*[-1 1],'-'); hold on
        set(hp,'color',V(jj).cmp);
    end
end; hold off 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xscale','log');
set(gca,'xlim',[0.5*xlm(1) 1.2*xlm(2)]);
%set(gca,'ylim',[50 20e3]/1e3);
xlabel('RT60 (s)')
ylabel('DRR (dB)')
title(PltPrm)
