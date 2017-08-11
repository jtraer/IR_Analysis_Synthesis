function PltIRStts_mRT60vDRR(Dh,PltPrm,V);

% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)
    
% find calibration recording
%C=[]; D=[];
%for jh=1:length(Dh);
%    if strcmp(Dh(jh).Meta.Env.Class,'CAL');
%        load(sprintf('%s/%s',Dh(jh).PthStm,Dh(jh).name));
%        C=[C; H];
%        if strcmp(H.Meta.App.PolarAngle_fromTop,'90')&&strcmp(H.Meta.App.AzimuthalAngle_fromFront,'0')
%            D=H;
%        end
%    end
%end

for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    mplt=zeros(1,length(tH));
    for jh=1:length(tH);
        mplt(:,jh)=median(tH(jh).RT60);
        mplt2(:,jh)=median(tH(jh).DRR);
    end
    plt=mean(mplt,2);
    err=std(mplt,[],2);
    plt2=mean(mplt2,2);
    err2=std(mplt2,[],2);
    % plot
    hp=plot(plt,plt2,[V(jj).mrk '-']); hold on
    set(hp,'linewidth',1,'markersize',6);
    set(hp,'color',V(jj).cmp);
    hp=plot(plt+err*[-1 1],plt2*ones(1,2),'-');
    set(hp,'color',V(jj).cmp);
    hp=plot(plt*ones(1,2),plt2+err2*[-1 1],'-');
    set(hp,'color',V(jj).cmp);
    hp=text(plt,plt2,1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
end; hold off 
set(gca,'xscale','log');
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%set(gca,'xlim',[0.5*xlm(1) 1.2*xlm(2)]);
ylabel('DRR (dB)')
xlabel('RT60 (s)')
title(PltPrm)
