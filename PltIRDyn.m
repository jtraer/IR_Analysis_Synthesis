function PltIRDyn(H);

%*** => compress the time series for plotting
h=H.nh;
hdot=diff(h)/H.fs;
hdd=diff(hdot)/H.fs;
h=sign(h).*abs(h).^(0.6);
hdot=sign(hdot).*abs(hdot).^(0.6);
hdd=sign(hdd).*abs(hdd).^(0.6);
h=h(1:end-2);
hdot=hdot(1:end-1);
hdot=hdot/max(abs(hdot));
hdd=hdd/max(abs(hdd));
Npts=10;
xx=linspace(min(h),max(h),Npts);
yy=linspace(min(hdot),max(hdot),Npts);
yy2=linspace(min(hdd),max(hdd),Npts);
Cnt=zeros(Npts,Npts);
Cnt2=zeros(Npts,Npts);
X=zeros(Npts,Npts);
Y=zeros(Npts,Npts);
Y2=zeros(Npts,Npts);
dX=zeros(Npts,Npts);
dY=zeros(Npts,Npts);
dY2=zeros(Npts,Npts);
T=zeros(Npts,Npts);
T2=zeros(Npts,Npts);
%measure phase
for jplt=1:(length(h)-1)
    [~,xndx]=min(abs(xx-h(jplt)));
    [~,yndx]=min(abs(yy-hdot(jplt)));
    [~,y2ndx]=min(abs(yy2-hdd(jplt)));
    % count
    Cnt(xndx,yndx)=Cnt(xndx,yndx)+1;
    Cnt2(xndx,y2ndx)=Cnt2(xndx,y2ndx)+1;
    cnt=Cnt(xndx,yndx);
    x=h(jplt);
    y=hdot(jplt);
    y2=hdd(jplt);
    dx=diff(h(jplt+[0:1]));  
    dy=diff(hdot(jplt+[0:1]));  
    dy2=diff(hdd(jplt+[0:1]));  
    t=jplt;
    X(xndx,yndx)=X(xndx,yndx)*(cnt-1)/cnt+x/cnt;
    Y(xndx,yndx)=Y(xndx,yndx)*(cnt-1)/cnt+y/cnt;
    Y2(xndx,y2ndx)=Y2(xndx,y2ndx)*(cnt-1)/cnt+y2/cnt;
    dX(xndx,yndx)=dX(xndx,yndx)*(cnt-1)/cnt+dx/cnt;
    dY(xndx,yndx)=dY(xndx,yndx)*(cnt-1)/cnt+dy/cnt;
    dY2(xndx,y2ndx)=dY2(xndx,y2ndx)*(cnt-1)/cnt+dy2/cnt;
    T(xndx,yndx)=T(xndx,yndx)*(cnt-1)/cnt+t/cnt;
    T2(xndx,y2ndx)=T(xndx,y2ndx)*(cnt-1)/cnt+t/cnt;
end
%*** Plot phase
ndx1=find(Cnt>0);
ndx2=find(Cnt2>0);
cmp=colormap(othercolor('Blues9',128));
cmp=flipud(cmp);
subplot(2,1,1);
for jplt=ndx1.';
    hp=quiver(X(jplt),Y(jplt),dX(jplt),dY(jplt)); hold on
    set(hp,'color',cmp(ceil(T(jplt)/length(h)*length(cmp)),:));
    set(hp,'linewidth',1);
    drawnow
end
subplot(2,1,2);
for jplt=ndx2.';
    hp=quiver(X(jplt),Y2(jplt),dX(jplt),dY2(jplt)); hold on
    set(hp,'color',cmp(ceil(T2(jplt)/length(h)*length(cmp)),:));
    set(hp,'linewidth',1);
    drawnow
end
subplot(2,1,1);
xlabel('Amplitude (compressed)');
ylabel('Velocity (compressed)')
subplot(2,1,2);
xlabel('Amplitude (compressed)');
ylabel('Acceleration (compressed)')
%set(gca,'xscale','log')
title([H.Path ': Denoised IR dynamics'])
