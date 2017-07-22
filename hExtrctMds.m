function M=hExtrctMds(H,MxNft);

% compute a spectrogram at a rangle of scales
Nft=2^ceil(log2(length(H.h))); % This is the longest possible spectrum we can compute
MxNft=min([Nft MxNft]);
Mxff=[1:MxNft/2]*H.fs/MxNft;
nft=256;
fftcnt=0;
while nft<=MxNft; fftcnt=fftcnt+1;
    fprintf('Nft=%d\n',nft);
    spc=fft(H.h,nft);
    spc=abs(spc(1:nft/2));
    ff=[1:nft/2]*H.fs/nft;
    spc(1)=min(abs(spc));
    spc(end)=min(abs(spc));

    Mds=zeros(MxNft/2,1);

    % Scroll through bins
    for jbn=2:length(ff)-1;
%        figure(101);
%        plot(20*log10(abs(spc)),ff); hold on;
        % fit a peak to this
        %** Initialize a peak
        options=optimset('display','off');
        
        Pk0=[ff(jbn) spc(jbn) ff(2)-ff(1)];
        %if fftcnt==1;
        %    [mx,mxndx]=max(spc);
        %    Pk0=[ff(mxndx) mx ff(2)-ff(1)];
        %else
        %    Pk0=Pk;
        %end
        Pkmn=[ff(1) min(spc) ff(2)-ff(1)];
        Pkmx=[ff(end) 2*max(spc) 1e3];
        Pk=fminsearchbnd(@(Pk) FtPk2Spc(spc,ff,Pk),Pk0,Pkmn,Pkmx,options);
        %Pk=fminbnd(@(Pk) FtPk2Spc(spc,ff,f0,spc(jbn),Pk),ff(2)-ff(1),(ff(end)-ff(1)),options);
        %** record peak
        [err,Pkk]=FtPk2Spc(spc,ff,Pk);
        FV=1-rms(spc-Pkk)/rms(spc);
        [~,ndx]=min(abs(Mxff-Pk(1)));
        bnwd=ceil(Pk(3)/(Mxff(2)-Mxff(1)));
        mn=max(ndx-bnwd,1);
        mx=min(ndx+bnwd,length(Mxff));
        Mds([mn:mx])=Mds([mn:mx])+FV*sin(linspace(0,pi,length([mn:mx]))).';

%        plot(20*log10(Pkk+min(abs(spc))),ff,'r'); 
%        drawnow; 
%        hold off

    end
    BMds(:,fftcnt)=Mds;
%    figure(102);
%    imagesc(20*log10(BMds+1e-24)); colorbar
%    axis xy;
%    drawnow; 

    nft=nft*2;

end

% Modes
Mds=sum(BMds,2);
%** find turning points
dMds=diff(Mds);
undx=find(dMds>0);
dndx=find(dMds<0);
%** scroll through and count modes
mcnt=0;
mn1ndx=1;
ndx=2; dr=1;
while ndx<=length(dMds);
    if ndx==length(dMds);
        dMds(ndx)=-dr*abs(dMds(ndx));
    end
    if sign(dMds(ndx))==-dr;
        tpndx=ndx;
        if dr==1;
            pkndx=tpndx;
        else
            mn2ndx=tpndx;
            %** we have a mode
            mcnt=mcnt+1;
            M(mcnt).cf=Mxff(pkndx);
            bnndx=min([(mn2ndx-pkndx) (pkndx-mn1ndx)]);
            M(mcnt).bw=min([abs(ff(mn2ndx)-ff(pkndx)) abs(ff(pkndx)-ff(mn1ndx))]);
            M(mcnt).FV=sum(Mds(pkndx)+[-1 1]*bnndx);
            mn1ndx=mn2ndx;
        end
        dr=-dr;
    end
    ndx=ndx+1;
end

figure(103);
plot(Mxff,20*log10(Mds+1e-24)); hold on
plot(Mxff(undx),20*log10(Mds(undx)),'o')
plot(Mxff(dndx),20*log10(Mds(dndx)),'o')
for jm=1:length(M);
    plot(M(jm).cf+[-1 1]*M(jm).bw,20*log10(M(jm).FV)*ones(1,2),'k+-');
end

