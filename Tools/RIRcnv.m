%* Takes two input itime series and convolves them together efficiently
%* Inputs:
%** x, fs   : time series 1 (i.e. a source) and its sampling rate
%** h, fh   : time series 2 (i.e. an impulse repsonse filter) and its sampling rate 
%** dns     : a downsampling factor, if operation speed is more important than audio quality

%* Outputs:
%** y, fy   : time series of the convolution and its sampling rate


function [y,fy]=RIRcnv(s,fs,h,fh,dns)

%* make sure the input time series are column vectors - column index refers to channel
if size(s,2)>size(s,1); 
    s=s.'; 
end
if size(h,2)>size(h,1); 
    h=h.'; 
end
%* select only the first channel 
if size(s,2)>1; 
    s=s(:,1); 
    fprintf('Warning: only one source channel will be convolved\n'); 
end
if size(h,2)>1; 
    h=h(:,1); 
    fprintf('Warning: only one source channel will be convolved\n'); 
end

%* resample if need be
if fs~=fh; 
    fprintf('warning: sampling frequencies differ\n'); 
    fprintf('sound fs=%3.0f Hz\n',fs); 
    fprintf('IR fh=%3.0f Hz\n',fh); 
    % interpolate the lower sampling frequency to a higher
    if fs<fh; 
        h=resample(h,fs,fh); fh=fs;
        fprintf('IR has been resampled\n')
    else
        s=resample(s,fh,fs); fs=fh;
        fprintf('sound has been resampled\n')
    end       
end
%* downsample, if necessary
if dns==[]; dns=1; end
s=decimate(s,dns); fs=fs/dns;
h=decimate(h,dns); fh=fh/dns;
%* save the new sampling frequency
fy=fs;
%* clean and zeropad both signals - to prevent wrap around artifacts
I=find(s~=0); s=s(min(I):max(I)); 
I=find(h~=0); h=h(min(I):max(I)); Nh=length(h);
s=[zeros(size(s)); s; zeros(size(s))];
h=[zeros(size(h)); h; zeros(size(h))];
%* create a delta function for comparison
d=zeros(size(h));
d(Nh+1)=1;

%* convolve
tic
%** first find teh smallest power of 2 to give a signal longer than BOTH inputs
Ns=length(s);
Nh=length(h);
Nft=2^(ceil(log(max(Ns,Nh))/log(2)));
%** Fourier transform
S=fft(s,Nft);
H=fft(h,Nft);
D=fft(d,Nft);
%** multiply and invert back to time domain
X=S.*H;
y=ifft(X); 
X=S.*D;
yd=ifft(X);
% check if the convolution is not centered and if so realign
if (abs(y(1))>1e-8||abs(y(end))>1e-8)
    y=fftshift(y);
end
% trim inaudible regions from the start and tail
I=find(abs(y)>max(abs(y))*1e-4);
y=y(min(I):max(I));
% rescale the amplitude of the convolution
y=y*max(abs(yd))/max(abs(y));
tcnv=toc;
fprintf('convolution performed: CPU time=%2.3fs\n',tcnv);
