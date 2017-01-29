T=5;

if ~isempty(instrfind)
    fclose(instrfind);
    delete(instrfind);
end

serialPort = '/dev/tty.usbmodem1411'; 
s = serial(serialPort,'BaudRate',9600);
s = serial(serialPort,'BaudRate',115200);
fopen(s); 

y=zeros(1,44100*T);
tm=y;

tic
cnt=0;
while toc<T; cnt=cnt+1;
    tmp=fscanf(s, '%d');
    if ~isempty(tmp);
        y(cnt) = tmp;
        tm(cnt)=toc;
    end
end
fclose(s)

y=y(1:cnt);
tm=tm(1:cnt);

dtt=linspace(1e-4,1e-2,1e4);
hst=hist(diff(tm),dtt);

figure; plot(tm,y)
figure; plot(1./dtt,hst)

