function h=RmvSpkTrnsFn(H,C);

h=H.h;
fs=H.fs;
for jc=1:length(C);
    Nf(jc)=length(C(jc).spc);
end
Nf(find(Nf==0))=[]; % remove the omnidirectional IR
[Nf,ndx]=min(Nf);
cff=C(ndx(1)).Spcff;
cSpc=zeros(Nf,1);
for jc=1:length(C);
    spc=interp1([C(jc).Spcff(:); C(jc).fs],[C(jc).spc(:); 0],cff(:));
    C(jc).Path
    size(cSpc)
    size(spc)
    cSpc=cSpc+spc(:)/length(C);
end

% the IRs have been filtered to remove all signal above and below the frequency bands of note (specified in Sb_fs).  Thus we cannot use these frequencies in filtering for they will explode for trivial reasons (i.e. they have been filtered out of the Calibration IRs). We will filter these away at the end anyways
% enforce the lower frequency values of the Calibration IR to have moderate values (such that they won't explode). 
lndx=find(cff<min(C(1).ff)); 
cSpc(lndx)=cSpc(max(lndx)+1).*(10.^linspace(2,0,length(lndx)));
% and the upper frequencies
undx=find(cff>max(C(1).ff)); 
cSpc(undx)=cSpc(min(undx)-1).*(10.^linspace(0,2,length(undx)));

% divide the Measured IR by the Calibration IR in the frequency domain
Npd=length(h);
tmp=[zeros(Npd,1); h; zeros(Npd,1)];
nft=2^max(ceil(log2([length(tmp) 2*length(cSpc)])));
tff=[1:nft]/nft*fs;
NH=fft(h,nft);
NH=NH(1:nft/2);
T=interp1([cff(:); tff(end)],[cSpc(:); cSpc(end) ],tff(1:nft/2));
T2=medfilt1(T,ceil(length(T)/(4*length(C(1).ff)))); % smooth to avoid deep notches in the calibration which will create audible artifacts in the output 
T=max([T; T2]);
fNH=NH./abs(sqrt(T(:)))*mean(abs(T)); % the sqrt is because we compute the spectrum with pwelch which computes power (i.e. amplitude squared)
h_cal=ifft([fNH; 0; flipud(conj(fNH(2:end)))]);
h_cal=h_cal(1:Npd);
h=h_cal;
