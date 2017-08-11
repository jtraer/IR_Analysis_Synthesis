function PltIRStts_Tgs(Dh,PltPrm,V);

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
        % if a calibration IR exists remove the speaker spectrum
        %if ~isempty(D);
        %    D_DRR=D.DRR;
        %    fprintf('Removing speaker DRR\n')
        %else 
        %    Dspc=zeros(size(tH(jh).DRR));
        %end
        mplt(:,jh)=tH(jh).Tgs*1e3;
    end
    plt=mean(mplt,2);
    err=std(mplt,[],2);
    % plot
    hp=errorbar(jj,plt,err,[V(jj).mrk '-']); hold on
    set(hp,'linewidth',1,'markersize',6);
    set(hp,'color',V(jj).cmp);
    hp=text(jj+0.1,plt+2,1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
end; hold off 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%set(gca,'xlim',[0.5*xlm(1) 1.2*xlm(2)]);
xlabel('Class')
ylabel('T_{G} (ms)')
title(PltPrm)
