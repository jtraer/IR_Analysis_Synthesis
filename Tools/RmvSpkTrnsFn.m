function h=RmvSpkTrnsFn(h,cSpc,cff);
    
Npd=length(h);
tmp=[zeros(Npd,1); h; zeros(Npd,1)];
nft=2^max(ceil(log2([length(tmp) length(cSpc)])));
tff=[1:nft]/nft*H.fs;
NH=fft(h,nft);
NH=NH(1:nft/2);
T=interp1([0 cff tff(end)],[cSpc(1); cSpc; cSpc(end) ],tff(1:nft/2));
fNH=NH./abs(T(:))*mean(abs(T));
h_cal=ifft([fNH; 0; flipud(conj(fNH(2:end)))]);
h_cal=h_cal(1:Npd);
h=h_cal;

