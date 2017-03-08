% Takes two inputs that can be a time series or filenames and convolves
% them together efficiently

% try the following for good sounds
%path(path,'/Users/jtraer/LabBook/Workshop/Proposals/NRSA/SoundExamples/DrySounds/Music/Drums/')
%path(path,'/Users/jtraer/LabBook/Workshop/Proposals/NRSA/SoundExamples/AIR/')

%function [y,Cgrm,tt,ff]=RIRcnv(s,h,fs,dns)
function [y,tt,fs]=RIRcnv(x,fx,h,fh,dns)

% if inputs are filenames load them
if ischar(x)==1;
    [s,fs,Nbtx]=wavread(x);
else
    fs=fx;
    s=x;
end
if ischar(h)==1;
    [h,fh,Nbth]=wavread(h);
end    
% make sure the right size
if size(s,2)>size(s,1); s=s.'; end
if size(h,2)>size(h,1); s=s.'; end
if size(s,2)>1; s=s(:,1); end
if size(h,2)>1; h=h(:,1); end

% resample if need be
if fs~=fh; 
    disp('warning: sampling frequencies differ'); 
    disp(['sound fs=' num2str(fs) 'Hz']); 
    disp(['AIR fs=' num2str(fh) 'Hz']); 
    % interpolate the lower sampling frequency to a higher
    if fs<fh; 
        h=resample(h,fs,fh); fh=fs;
        disp('AIR has been resampled')
    else
        s=resample(s,fh,fs); fs=fh;
        disp('warning: sound has been resampled')
    end       
end
%downsample
if dns==[]; dns=1; end
s=decimate(s,dns); fs=fs/dns;
h=decimate(h,dns); fh=fh/dns;
% clean and zeropad
Npd=round(2*fs);
I=find(s~=0); s=s(min(I):max(I));
I=find(h~=0); h=h(min(I):max(I));
s=[s; zeros(Npd,1)];
h=[h; zeros(Npd,1)];


%%
tic
% convolve
Ns=length(s);
Nh=length(h);
Nft=2^(ceil(log(max(Ns,Nh))/log(2)));
S=fft(s,Nft);
H=fft(h,Nft);
X=S.*H;
y=ifft(X);
jndx=0; prc=1;
while prc>1e-8
    jndx=jndx+round(fs/32);
    prc=sum(abs(y(jndx:end)))/sum(abs(y));
end    
y=y(1:jndx);
Ny=length(y);
y=y./(abs(max(y)))*0.95;
tt=0:1/fs:(Ny-1)/fs;
tcnv=toc;
disp(['convolution performed: CPU time=' num2str(tcnv) 's'])
disp('try:')
disp('sound(y,fy)')

%% cochlear domain
% tic
% [Cgrm,ff,fltbnk]=Snd2Cgrm(y,fs,40,400);
% tch=toc;
% disp(['Cochleagrams computed: CPU time=' num2str(tch) 's'])
