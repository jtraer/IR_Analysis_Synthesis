function Y=CnvExcFrc(H,P,TT,Tp,Pth,Nm)
x=P.h;
%%%%%%% CHECK IF WE NEED THIS - IT'S UGLY!!!!
if length(x)>101;
    x(end+[-500:0])=x(end+[-500:0]).*linspace(1,0,501).';
end
x=[x; zeros(ceil(H.fs/2),1)];

% make the necessary output folder
unix(sprintf('mkdir -p %s',Pth))

Y(1).y=x;
Y(1).Tp='IR';
Y(1).T=0;

% scroll through Time entries
cnt=length(Y);
for jT=TT; cnt=cnt+1;
    % create time vector
    tt=[0:1/H.fs:jT/1e3];
    % choose a Force type
    if strcmp(Tp,'Impulse')
        F=sin(pi*tt/max(tt));
    elseif strcmp(Tp,'Gauss')
        F=randn(size(tt));
    end
    F=F(:);
    % perfrom convolution
    y=RIRcnv(x,H.fs,F,H.fs,1);
    % save
    Y(cnt).y=y;
    Y(cnt).Tp=Tp;
    Y(cnt).T=jT;
end

% write
close all
fcnt=0;
for jy=1:length(Y)
    y=Y(jy).y;
    if jy<=2;
        Mx=max(abs(y));
        y=y/Mx;
    else
        y=y/Mx*(Y(2).T/Y(jy).T);
    end
    y=y*(1-1e-6);
    y=[zeros(ceil(H.fs/5),1); y];
    audiowrite(sprintf('%s/%s_%s_%04d.wav',Pth,Nm,Y(jy).Tp,round(Y(jy).T*10)),y,H.fs,'BitsPerSample',24);

    %** => plot time Series
    fcnt=fcnt+1; figure(fcnt)
    %*** => compress the time series for plotting
    h=sign(y).*abs(y).^(0.6);
    %*** Plot
    plot([1:length(y)]/H.fs,y);
    %*** => plot scale lines
    hold on
    for jln=1:3
        plot(([2 length(y)]/H.fs),10^(-jln*0.6)*ones(1,2),'k:');
        plot(([2 length(y)]/H.fs),-10^(-jln*0.6)*ones(1,2),'k:');
    end
    xlabel('Time (s)');
    ylabel('Waveform amplitude (compressed)')
    set(gca,'xscale','log')
    title(sprintf('%s - %04d',Y(jy).Tp,Y(jy).T));
    %set(gca,'fontsize',fntsz);
    saveas(gcf,sprintf('%s/%s_%s_%04d_IR',Pth,Nm,Y(jy).Tp,round(Y(jy).T*10)),'jpg') 
    
    %*** Plot Cgram
    fcnt=fcnt+1; figure(fcnt)
    Npts=length(y);
    [fltbnk,ff,erbff]=make_erb_cos_filters(3*Npts,H.fs,length(H.ff),H.ff(1),H.ff(end));
    C_s=generate_subbands([zeros(Npts,1); y; zeros(Npts,1)].',fltbnk);
    C_s=C_s(Npts+[1:Npts],:).'; 
    C_s=C_s(2:end-1,:);

    plt=20*log10(abs(C_s)); 
    pcolor([1:Npts]/H.fs,H.ff/1e3,plt);
    axis xy; shading flat
    xlabel('Time (s)');
    ylabel('Frequency (kHz)');
    title([H.Path ': Synth IR']);
    set(gca,'clim',max(max(plt))+[-80 0]);
    set(gca,'xlim',[0 Npts/H.fs])
    set(gca,'yscale','log');
    colorbar
    eval(sprintf('colormap(othercolor(''%s'',64));',P.cmp));
    %set(gca,'fontsize',fntsz);
    saveas(gcf,sprintf('%s/%s_%s_%04d_Cgrm',Pth,Nm,Y(jy).Tp,round(Y(jy).T*10)),'jpg') 
end
