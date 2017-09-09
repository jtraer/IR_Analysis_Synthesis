function PltIRStts_MdJntDts(Dh,PltPrm,V)
% preallocate one data point for each class for the legend
%subplot(1,2,1); MkLgnd(V,Dh,PltPrm); hold on

% Joint distribution of Mode frequencie and RT60s
Npts=50;
JntDst=zeros(Npts,Npts,2);
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
    MdOP=[];
    for jh=1:length(tH);
        Mdf=[Mdf [tH(jh).Modes.cf]];
        MdR=[MdR [tH(jh).Modes.RT60]];
        MdOP=[MdOP [tH(jh).Modes.OnPwr]];
    end
    % remove any high frequency peaks (This is a HACK)
    ndx=find(Mdf>10e3);
    Mdf(ndx)=[];
    MdR(ndx)=[];
    MdOP(ndx)=[];
    % compute histogram
    JDff=e2freq(linspace(freq2e(tH(1).ff(1)),freq2e(10e3),Npts));
    JDR=(logspace(-2,0.5,Npts));
    JDOP=(linspace(-100,20,Npts));
    for jmd=1:length(Mdf);
        [~,fndx]=min(abs(JDff-Mdf(jmd)));
        [~,Rndx]=min(abs(JDR-MdR(jmd)));
        [~,OPndx]=min(abs(JDOP-MdOP(jmd)));
        JntDst(fndx,Rndx,1)=JntDst(fndx,Rndx,1)+1;
        JntDst(fndx,Rndx,2)=JntDst(fndx,OPndx,2)+1;
    end
    for jp=1:2;
        JntDst(:,:,jp)=medfilt2(JntDst(:,:,jp),[2 2]);
    end
    % find the maximum
    tmp=JntDst(:,:,1);
    [~,ndx]=max(tmp(:));
    [fndx,Rndx]=ind2sub([Npts,Npts],ndx);
    tmp=JntDst(:,:,2);
    [~,ndx]=max(tmp(:));
    [fndx2,OPndx]=ind2sub([Npts,Npts],ndx);
    % plot this class
    subplot(1,2,1)
    contour(JDR,JDff/1e3,JntDst(:,:,1),[0.2:0.2:1]*max(max(JntDst(:,:,1))),'linecolor',V(jj).cmp)
    hold on;
    hp=plot(JDR(Rndx),JDff(fndx)/1e3,V(jj).mrk);
    set(hp,'color',V(jj).cmp);
    hp=text(JDR(Rndx),JDff(fndx)/1e3,1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);
    
    subplot(1,2,2)
    contour(JDOP,JDff/1e3,JntDst(:,:,2),[0.2:0.2:1]*max(max(JntDst(:,:,1))),'linecolor',V(jj).cmp)
    hold on;
    hp=plot(JDOP(OPndx),JDff(fndx2)/1e3,V(jj).mrk);
    set(hp,'color',V(jj).cmp);
    hp=text(JDOP(OPndx),JDff(fndx2)/1e3,1.001,V(jj).name); 
    set(hp,'color',V(jj).cmp);    
end;
subplot(1,2,1)
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Mode RT60 (s)')
ylabel('Mode Frequency (kHz)')
subplot(1,2,2)
hold off
axis tight; xlm=get(gca,'xlim'); ylm=get(gca,'ylim');
set(gca,'yscale','log')
xlabel('Mode Onset Power (dB)')

