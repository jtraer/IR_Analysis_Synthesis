function PltIRStts_smKrt(Dh,PltPrm,V);

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
    smtt=[1/H.fs:0.001:(max(tt)-0.01)];
    smkrt=zeros(length(smtt),length(tH));
    krt=zeros(length(smtt),length(tH));
    bnwd=ceil(0.01*tH(1).fs);

    for jh=1:length(tH);
        jkrt=tH(jh).krt;
        mxndx=prctile(jkrt,99);
        mxndx=find(jkrt>mxndx);
        jkrt=jkrt(mxndx(1):end);
        for jb=1:length(smtt)
            [~,ndx]=min(abs(tt-smtt(jb)));
            if ndx<length(jkrt)-bnwd;
                smkrt(jb,jh)=length(find(jkrt(ndx+[0:bnwd])>3.4))/bnwd; 
            else 
                smkrt(jb,jh)=0;
            end
        end
    end
    plt=median(smkrt,2);
    err=std(smkrt,[],2);
    % plot
    hp=plot(smtt,plt,['-']); hold on
    set(hp,'linewidth',1,'markersize',6);
    set(hp,'color',V(jj).cmp);
    [mx,mxndx]=max(plt);
    hp=text(smtt(mxndx(1))+0.1,mx(1),1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
    hp=plot(smtt,plt+err,':'); hold on
    set(hp,'color',V(jj).cmp);
    hp=plot(smtt,plt-err,':'); hold on
    set(hp,'color',V(jj).cmp);
end; 
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'xscale','log');
set(gca,'xlim',[0.5*xlm(1) 1.2*xlm(2)]);
%set(gca,'yscale','log');ylm=get(gca,'ylim'); 
% add distance marks 
for jd=[1 3 10 30];
    plot(jd/340*[1 1],ylm,'k:')
    text(jd/340,2*ylm(1),sprintf('%dm',jd)); 
end
xlabel('time (s)')
ylabel('Fraction of 10ms windows with Early Reflections (Kurt>3.4)')
title(PltPrm)
