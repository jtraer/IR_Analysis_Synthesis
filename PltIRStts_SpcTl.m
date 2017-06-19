function PltIRStts_SpcTl(Dh,PltPrm,V);

% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)

% find calibration recording
C=[]; D=[];
for jh=1:length(Dh);
    if strcmp(Dh(jh).Meta.Env.Class,'CAL');
        load(sprintf('%s/%s',Dh(jh).PthStm,Dh(jh).name));
        C=[C; H];
        if strcmp(H.Meta.App.PolarAngle_fromTop,'90')&&strcmp(H.Meta.App.AzimuthalAngle_fromFront,'0')
            D=H;
        end
    end
end
    
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    ff=H.Attck(6).ff/1e3; mplt=zeros(length(ff),length(tH));
    for jh=1:length(tH);
        % if a calibration IR exists remove the speaker spectrum
        if ~isempty(D);
            Dspc=20*log10(D.Attck(6).Spc);
            fprintf('Removing speaker spectrum\n')
        else 
            Dspc=zeros(size(tH(jh).Attck(6).Spc));
        end
        mplt(:,jh)=20*log10(tH(jh).Attck(6).Spc)-Dspc;
    end
    plt=median(mplt,2);
    err=std(mplt,[],2);
    % smooth for plotting purposes
    plt=medfilt1(plt,16);
    err=medfilt1(err,16);
    % plot
    hp=plot(plt,ff,['-']); hold on
    set(hp,'linewidth',1,'markersize',6);
    set(hp,'color',V(jj).cmp);
    [mx,mxndx]=max(plt);
    hp=text(mx+0.1,ff(mxndx),1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
    hp=plot(plt+err,ff,':'); hold on
    set(hp,'color',V(jj).cmp);
    hp=plot(plt-err,ff,':'); hold on
    set(hp,'color',V(jj).cmp);
end; hold off 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xlim',[0.5*xlm(1) 1.2*xlm(2)]);
set(gca,'yscale','log');
set(gca,'ylim',[20 20e3]/1e3);
xlabel('Power (dB)')
ylabel('Frequency (kHz)')
title(PltPrm)
