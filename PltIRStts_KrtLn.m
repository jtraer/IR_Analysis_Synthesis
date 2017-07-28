function PltIRStts_KrtLn(Dh,PltPrm,V);

% preallocate one data point for each class for the legend
MkLgnd(V,Dh,PltPrm)
    
for jj=1:length(V);
    % collate all IRs that have this particular label
    tH=[]; 
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); load(''%s/%s''); tH=[tH H]; end;',PltPrm,Dh(jh).PthStm,Dh(jh).name));
    end
    % specify the ordinates and abscissa
    for jh=1:length(tH);
        Npts(jh)=length(tH(jh).krt);
    end
    Npts=max(Npts);
    tt=[1:Npts]/H.fs;
    krt=zeros(Npts,length(tH));
    for jh=1:length(tH);
        jkrt=tH(jh).krt;
        mxndx=prctile(jkrt,99);
        mxndx=find(jkrt>mxndx);
        jkrt=jkrt(mxndx(1):end);
        krt(1:length(jkrt),jh)=jkrt;
        krt(length(jkrt):end,jh)=median(jkrt);
    end
    plt=median(krt,2);
    err=std(krt,[],2);
    % downsample
    ntt=linspace(min(tt),max(tt),max([max(tt)*1e2 10]));
    plt=interp1(tt,plt,ntt);
    err=interp1(tt,err,ntt);
    tt=ntt;
    % plot
    hp=plot(tt,plt,['-']); hold on
    set(hp,'linewidth',1,'markersize',6);
    set(hp,'color',V(jj).cmp);
    [mx,mxndx]=max(plt);
    hp=text(tt(mxndx)+0.1,mx,1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
    hp=plot(tt,plt+err,':'); hold on
    set(hp,'color',V(jj).cmp);
    hp=plot(tt,plt-err,':'); hold on
    set(hp,'color',V(jj).cmp);
end; 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
%set(gca,'xscale','log');
%set(gca,'xlim',[2e-3 0.2]);
set(gca,'yscale','log');ylm=get(gca,'ylim');
set(gca,'ylim',[1 ylm(2)]);
% add distance marks 
for jd=[1 3 10 30];
    plot(jd/340*[1 1],ylm,'k:')
    text(jd/340,2*ylm(1),sprintf('%dm',jd));
end
xlabel('time (s)')
ylabel('Kurtosis')
title(PltPrm)
