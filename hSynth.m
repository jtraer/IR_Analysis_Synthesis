function H=hSynth(H,ftp);
set(0,'DefaultFigureVisible','off');
%* == hSynth.m i.e. Make synthetics ==
%* Takes a structure of IR properties (i.e. H, output by hPrp.m), and makes a family of synthetics

rcnt=0;
%rcnt=rcnt+1; R(rcnt).Name='R_SpcCnst';    R(rcnt).cmp='BuPu3';
%rcnt=rcnt+1; R(rcnt).Name='R_SpcRot';    R(rcnt).cmp='BuPu3';
%rcnt=rcnt+1; R(rcnt).Name='R_Lng2';    R(rcnt).cmp='BuPu3';
%rcnt=rcnt+1; R(rcnt).Name='R_Lng3';    R(rcnt).cmp='BuPu3';
%rcnt=rcnt+1; R(rcnt).Name='R_Lng4';    R(rcnt).cmp='BuPu3';
%rcnt=rcnt+1; R(rcnt).Name='R_Lng5';    R(rcnt).cmp='BuPu3';
rcnt=rcnt+1; R(rcnt).Name='R_Shrt';    R(rcnt).cmp='BuPu3';
%rcnt=rcnt+1; R(rcnt).Name='R_NoMod';      R(rcnt).cmp='YlGn3';
%rcnt=rcnt+1; R(rcnt).Name='R_Lin';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_MrSym';        R(rcnt).cmp='BuPu9';
%rcnt=rcnt+1; R(rcnt).Name='R_LssSym';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_HiSpc';        R(rcnt).cmp='BuPu9';
%rcnt=rcnt+1; R(rcnt).Name='R_LoSpc';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_Sn';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_Tp10';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_Tp5';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_RmvMds';        R(rcnt).cmp='Reds3';
%rcnt=rcnt+1; R(rcnt).Name='R_AddMds';        R(rcnt).cmp='Reds3';

scnt=0;
scnt=scnt+1; S(scnt).Name='S_Exp';        S(scnt).cmp='Greens9';
scnt=scnt+1; S(scnt).Name='S_Mds';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_SbMds';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_MvMds';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_MdsSpcCnst';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_MdsRndAmp';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_MdsRndDcy';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_FtMds';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_LngMds2';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_LngMds3';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_LngMds4';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_LwMds';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_HiMds';        S(scnt).cmp='Greens9';
%scnt=scnt+1; S(scnt).Name='S_LinEnd';     S(scnt).cmp='Reds9';
%scnt=scnt+1; S(scnt).Name='S_LinStrt';    S(scnt).cmp='Oranges9';
%scnt=scnt+1; S(scnt).Name='S_SpcCnst';    S(scnt).cmp='BuPu9';

fntsz=15;

%* == Compute the cochleagram ==
%** zeropad to avoid edge effects
Npts=length(H.h);
Nbnds=length(H.ff);
[fltbnk,ff,erbff]=make_erb_cos_filters(3*Npts,H.fs,Nbnds,H.ff(1),H.ff(end));
Cgrm=generate_subbands([zeros(Npts,1); H.nh; zeros(Npts,1)].',fltbnk);
Cgrm=Cgrm(Npts+[1:Npts],:).'; 
%** Remove the extreme bands
Cgrm=Cgrm([2:(end-1)],:);
H.ff=ff([2:(end-1)]);

%** rank the bands in terms of power
CgrmP=sum(Cgrm.^2,2);
[~,Pordr]=sort(CgrmP,'descend');

tt=[1:length(H.nh)]/H.fs;
%* == Scroll through cochlear channels ==
for jbn=1:Nbnds; 
    % Extract the subband
    tmp=Cgrm(jbn,:); 
    % take envelope
    tmp2=abs(hilbert([zeros(1,Npts) tmp zeros(1,Npts)]));  
    tmp2=tmp2(Npts+[1:Npts]);
    % compute a new altered subband
    for jr=1:length(R)
        R(jr).Name
        %** Spectrally constant
        if strcmp('R_SpcCnst',R(jr).Name);
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(median(H.RT60));
            ImpDcy=DesDcy-Dcy;
            tst=tmp.*10.^(-ImpDcy/20*tt);  
            if abs(tst(end))<abs(tmp(end))
                R(jr).Cgrm(jbn,:)=tst;
            else
                R(jr).Cgrm(jbn,:)=tmp;
            end
        %** Spectral rotation
        elseif strcmp('R_SpcRot',R(jr).Name)
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(H.RT60(end-jbn+1));
            ImpDcy=DesDcy-Dcy;
            tst=tmp.*10.^(-ImpDcy/20*tt);  
            if abs(tst(end))<abs(tmp(end))
                R(jr).Cgrm(jbn,:)=tst;
            else
                R(jr).Cgrm(jbn,:)=tmp;
            end
        %** Slower decay (Longer IR)
        elseif strcmp('R_Lng2',R(jr).Name)
            fct=2;
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(fct*H.RT60(jbn));
            ImpDcy=DesDcy-Dcy;
            tmp2=zeros(ceil(fct+1)*Npts,1);
            tmp3=tmp.*(10.^(Dcy/20*tt)); 
            tmp2(1:Npts)=tmp3;
            tmp2(floor(Npts/2):Npts)=tmp2(floor(Npts/2):Npts).*cos(linspace(0,pi/2,length(tmp2(floor(Npts/2):Npts))).');
            for jj=1:ceil(fct)*2;
                tmp2(floor(Npts*jj/2)-1+[1:Npts])=tmp2(floor(Npts*jj/2)-1+[1:Npts])+tmp3(:).*sin(linspace(0,pi,Npts)).';
            end
            % test to make sure this is flat and not growing
            pwr=20*log10(abs(tmp2));
            pwr(find(pwr<median(pwr)))=median(pwr);
            p=polyfit([1:length(tmp2)].'/H.fs,pwr,1);
            tmp2=tmp2.*(10.^(-p(1)*[1:length(tmp2)].'/H.fs/20));
            % add desired decay
            tst=tmp2.*10.^(-DesDcy/20*[1:length(tmp2)].'/H.fs);  
            R(jr).Cgrm(jbn,:)=tst;
        elseif strcmp('R_Lng3',R(jr).Name)
            fct=3;
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(fct*H.RT60(jbn));
            ImpDcy=DesDcy-Dcy;
            tmp2=zeros(ceil(fct+1)*Npts,1);
            tmp3=tmp.*(10.^(Dcy/20*tt)); 
            tmp2(1:Npts)=tmp3;
            tmp2(floor(Npts/2):Npts)=tmp2(floor(Npts/2):Npts).*cos(linspace(0,pi/2,length(tmp2(floor(Npts/2):Npts))).');
            for jj=1:ceil(fct)*2;
                tmp2(floor(Npts*jj/2)-1+[1:Npts])=tmp2(floor(Npts*jj/2)-1+[1:Npts])+tmp3(:).*sin(linspace(0,pi,Npts)).';
            end
            % test to make sure this is flat and not growing
            pwr=20*log10(abs(tmp2));
            pwr(find(pwr<median(pwr)))=median(pwr);
            p=polyfit([1:length(tmp2)].'/H.fs,pwr,1);
            tmp2=tmp2.*(10.^(-p(1)*[1:length(tmp2)].'/H.fs/20));
            % add desired decay
            tst=tmp2.*10.^(-DesDcy/20*[1:length(tmp2)].'/H.fs);  
            R(jr).Cgrm(jbn,:)=tst;
        elseif strcmp('R_Lng4',R(jr).Name)
            fct=4;
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(fct*H.RT60(jbn));
            ImpDcy=DesDcy-Dcy;
            tmp2=zeros(ceil(fct+1)*Npts,1);
            tmp3=tmp.*(10.^(Dcy/20*tt)); 
            tmp2(1:Npts)=tmp3;
            tmp2(floor(Npts/2):Npts)=tmp2(floor(Npts/2):Npts).*cos(linspace(0,pi/2,length(tmp2(floor(Npts/2):Npts))).');
            for jj=1:ceil(fct)*2;
                tmp2(floor(Npts*jj/2)-1+[1:Npts])=tmp2(floor(Npts*jj/2)-1+[1:Npts])+tmp3(:).*sin(linspace(0,pi,Npts)).';
            end
            % test to make sure this is flat and not growing
            pwr=20*log10(abs(tmp2));
            pwr(find(pwr<median(pwr)))=median(pwr);
            p=polyfit([1:length(tmp2)].'/H.fs,pwr,1);
            tmp2=tmp2.*(10.^(-p(1)*[1:length(tmp2)].'/H.fs/20));
            % add desired decay
            tst=tmp2.*10.^(-DesDcy/20*[1:length(tmp2)].'/H.fs);  
            R(jr).Cgrm(jbn,:)=tst;
        elseif strcmp('R_Lng5',R(jr).Name)
            fct=5;
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(fct*H.RT60(jbn));
            ImpDcy=DesDcy-Dcy;
            tmp2=zeros(ceil(fct+1)*Npts,1);
            tmp3=tmp.*(10.^(Dcy/20*tt)); 
            tmp2(1:Npts)=tmp3;
            tmp2(floor(Npts/2):Npts)=tmp2(floor(Npts/2):Npts).*cos(linspace(0,pi/2,length(tmp2(floor(Npts/2):Npts))).');
            for jj=1:ceil(fct)*2;
                tmp2(floor(Npts*jj/2)-1+[1:Npts])=tmp2(floor(Npts*jj/2)-1+[1:Npts])+tmp3(:).*sin(linspace(0,pi,Npts)).';
            end
            % test to make sure this is flat and not growing
            pwr=20*log10(abs(tmp2));
            pwr(find(pwr<median(pwr)))=median(pwr);
            p=polyfit([1:length(tmp2)].'/H.fs,pwr,1);
            tmp2=tmp2.*(10.^(-p(1)*[1:length(tmp2)].'/H.fs/20));
            % add desired decay
            tst=tmp2.*10.^(-DesDcy/20*[1:length(tmp2)].'/H.fs);  
            R(jr).Cgrm(jbn,:)=tst;
        %** Fatser decay (shorter IR)
        elseif strcmp('R_Shrt',R(jr).Name)
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(H.RT60(jbn)/2);
            ImpDcy=DesDcy-Dcy;
            tst=tmp.*10.^(-ImpDcy/20*tt);  
            R(jr).Cgrm(jbn,:)=tst;
        %**Make more symmetircal 
        elseif strcmp('R_MrSym',R(jr).Name)
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(H.RT60(jbn)*1*(jbn/(Nbnds+1)));
            ImpDcy=DesDcy-Dcy;
            tst=tmp.*10.^(-ImpDcy/20*tt);  
            R(jr).Cgrm(jbn,:)=tst;
        %**Make less symmetircal 
        elseif strcmp('R_LssSym',R(jr).Name)
            Dcy=60/(H.RT60(jbn));
            DesDcy=60/(H.RT60(jbn)*1*((Nbnds-jbn)/(Nbnds+1)));
            ImpDcy=DesDcy-Dcy;
            tst=tmp.*10.^(-ImpDcy/20*tt);  
            R(jr).Cgrm(jbn,:)=tst;
        %** replace bins with sinusoids
        elseif strcmp('R_Sn',R(jr).Name)
            if jbn==1;
                h=zeros(size(tt));
            end
            Dcy=60/(H.RT60(jbn));
            x=sin(2*pi*H.ff(jbn)*tt);
            x=x*10.^(-H.DRR(jbn)/20);
            x=x.*(10.^(-Dcy/20*tt));
            h=h+x;
            R(jr).Cgrm(jbn,:)=tmp;
        %** zero out smallest bins 
        elseif strcmp('R_Tp10',R(jr).Name)
            ndx=find(jbn==Pordr);
            if ndx>10;
                R(jr).Cgrm(jbn,:)=zeros(size(tmp));
            else
                R(jr).Cgrm(jbn,:)=tmp;
            end
        elseif strcmp('R_Tp5',R(jr).Name)
            ndx=find(jbn==Pordr);
            if ndx>5;
                R(jr).Cgrm(jbn,:)=zeros(size(tmp));
            else
                R(jr).Cgrm(jbn,:)=tmp;
            end
        elseif strcmp('R_RmvMds',R(jr).Name)
            R(jr).Cgrm(jbn,:)=tmp;
        elseif strcmp('R_AddMds',R(jr).Name)
            R(jr).Cgrm(jbn,:)=tmp;
        %** Increase onset power in high frequencies
        elseif strcmp('R_HiSpc',R(jr).Name)
            R(jr).Cgrm(jbn,:)=tmp*10^((jbn-Nbnds/2)/Nbnds*30/20);
        %** Increase onset power in low frequencies
        elseif strcmp('R_LoSpc',R(jr).Name)
            R(jr).Cgrm(jbn,:)=tmp*10^(-(jbn-Nbnds/2)/Nbnds*30/20);
        %** Remove Modulations
        elseif strcmp('R_NoMod',R(jr).Name)
            Dcy=20*log10(abs(tmp2));
            DesDcy=H.DRR(jbn)-60/H.RT60(jbn)*tt;
            R(jr).Cgrm(jbn,:)=tmp.*10.^((DesDcy-Dcy)/20);
        %** impose linear decay
        elseif strcmp('R_Lin',R(jr).Name)
            Dcy=20*log10(abs(tmp));
            Nln=ceil(H.RT60(jbn)*(60+H.DRR(jbn))/60*H.fs);
            if Nln>Npts; Nln=Npts; 
            elseif Nln<1; Nln=1; 
            end
            DesDcy=linspace(10^(H.DRR(jbn)/20),10^(-60/20),Nln);    
            DesDcy=20*log10(DesDcy);
            DesDcy=[DesDcy -100*ones(1,Npts-Nln)];
            R(jr).Cgrm(jbn,:)=tmp.*10.^((DesDcy-Dcy)/20);
        end
    end
end

%keyboard

% resynthesize the remixed IRs
for jr=1:length(R)
    rmxCgrm=R(jr).Cgrm;
    jNpts=size(R(jr).Cgrm,2);
    [fltbnk,ff,erbff]=make_erb_cos_filters(3*jNpts,H.fs,Nbnds,H.ff(1),H.ff(end));
    rmxCgrm=[zeros(1,size(rmxCgrm,2)); rmxCgrm; zeros(1,size(rmxCgrm,2))];
    rmxh=collapse_subbands([zeros(size(rmxCgrm)) rmxCgrm zeros(size(rmxCgrm))].',fltbnk);
    rmxh=rmxh(jNpts+[1:jNpts]);
    rmxh=rmxh/max(abs(rmxh));
    R(jr).h=rmxh;
    if strcmp('R_Sn',R(jr).Name)
        R(jr).h=h(:);
    end
    if strcmp('R_RmvMds',R(jr).Name)
        h=R(jr).h;
        for jm=1:length(H.Modes);
            if H.Modes(jm).cf-H.Modes(jm).Wd/2>0;
                [b,a]=iirnotch(H.Modes(jm).cf/(H.fs/2),H.Modes(jm).Wd/2/(H.fs/2));
                h=filtfilt(b,a,h);
            end
        end
        R(jr).h=h;
    end
    if strcmp('R_AddMds',R(jr).Name)
        h=R(jr).h;
        for jm=1:length(H.Modes);
            md=sin(2*pi*H.Modes(jm).cf*(1+randn(1)/10)*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/H.Modes(jm).RT60/20));
            h=h+md(:);
        end
        R(jr).h=h;
    end
end

%* Make synthetics
for js=1:length(S)
    if strcmp('S_Exp',S(js).Name)
        [h_s,C_s]=hMakeSynth(H.DRR,60./H.RT60,H.ff([1 end]),Npts,H.fs,-70,50,'Exp'); 
    elseif strcmp('S_LinEnd',S(js).Name)  
        [h_s,C_s]=hMakeSynth(H.DRR,60./H.RT60,H.ff([1 end]),Npts,H.fs,-70,50,'Lin2'); 
    elseif strcmp('S_LinStrt',S(js).Name)  
        [h_s,C_s]=hMakeSynth(H.DRR,60./H.RT60,H.ff([1 end]),Npts,H.fs,-70,50,'Lin3'); 
    elseif strcmp('S_SpcCnst',S(js).Name)  
        [h_s,C_s]=hMakeSynth(median(H.DRR)*ones(size(H.DRR)),60./median(H.RT60)*ones(size(H.RT60)),H.ff([1 end]),Npts,H.fs,-70,50,'Exp'); 
    elseif strcmp('S_Mds',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes);
            md=sin(2*pi*H.Modes(jm).cf*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/H.Modes(jm).RT60/20));
            h_s=h_s+md(:);
            sPth=sprintf('%s/Synth%03dBnds',H.Path,Nbnds); unix(sprintf('mkdir -p %s',sPth));
            md=[zeros(ceil(H.fs/5),1); md(:)];
            audiowrite(sprintf('%s/Md%02d_%04d.wav',sPth,cnt,round(H.Modes(jm).cf)),md./max(abs(md))*0.99,H.fs,'BitsPerSample',24);
            ths=[zeros(ceil(H.fs/5),1); h_s(:)];
            audiowrite(sprintf('%s/cMd%02d_%04d.wav',sPth,cnt,round(H.Modes(jm).cf)),ths./max(abs(ths))*0.99,H.fs,'BitsPerSample',24);
        end
    elseif strcmp('S_SbMds',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        I=rand(size(srt));
        srt(find(I<0.5))=[];
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes);
            md=sin(2*pi*H.Modes(jm).cf*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/H.Modes(jm).RT60/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_MvMds',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes);
            I=ceil(rand(1)*length(H.Modes));
            f=H.Modes(I).cf;
            md=sin(2*pi*H.Modes(I).cf*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/H.Modes(jm).RT60/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_MdsSpcCnst',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes);
            md=sin(2*pi*H.Modes(jm).cf*tt);
            md=md*10^(median([H.Modes.OnPwr])/20);
            md=md.*(10.^(-tt*60/median([H.Modes.RT60])/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_MdsRndAmp',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes);
            md=sin(2*pi*H.Modes(jm).cf*tt); 
            md=md*10^((20*rand(1)+H.Modes(jm).OnPwr)/20); 
            md=md.*(10.^(-tt*60/H.Modes(jm).RT60/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_MdsRndDcy',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes);
            md=sin(2*pi*H.Modes(jm).cf*tt);
            md=md*10^((H.Modes(jm).OnPwr)/20);
            md=md.*(10.^(-tt*60/H.Modes(jm).RT60/(1+max([0.8 randn(1)/5]))/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_FtMds',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes)
            for jj=1:5;
                md=sin(2*pi*(H.Modes(jm).cf*(1+randn(1)/5))*tt);
                md=md*10^(H.Modes(jm).OnPwr*(1+randn(1)/5)/20);
                md=md.*(10.^(-tt*60/(H.Modes(jm).RT60*(1+randn(1)/5))/20));
                h_s=h_s+md(:);
            end
        end
    elseif strcmp('S_LngMds2',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes)
            md=sin(2*pi*(H.Modes(jm).cf)*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/(H.Modes(jm).RT60*2)/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_LngMds3',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes)
            md=sin(2*pi*(H.Modes(jm).cf)*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/(H.Modes(jm).RT60*3)/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_LngMds4',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes)
            md=sin(2*pi*(H.Modes(jm).cf)*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/(H.Modes(jm).RT60*4)/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_LwMds',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes)
            md=sin(2*pi*(H.Modes(jm).cf/2)*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/(H.Modes(jm).RT60)/20));
            h_s=h_s+md(:);
        end
    elseif strcmp('S_HiMds',S(js).Name)
        h_s=zeros(Npts,1);
        tt=[1:Npts]/H.fs;
        % scroll through modes
        [~,srt]=sort([H.Modes.MnPwr],'descend');
        cnt=0;
        for jm=srt; cnt=cnt+1; %1:length(H.Modes)
            md=sin(2*pi*(H.Modes(jm).cf*2)*tt);
            md=md*10^(H.Modes(jm).OnPwr/20);
            md=md.*(10.^(-tt*60/(H.Modes(jm).RT60)/20));
            h_s=h_s+md(:);
        end
    end
    S(js).h=h_s(:);
    S(js).Cgrm=C_s;
end

%* == Plot ==
P1.h=H.nh;
P1.Name='RW';
P1.cmp='Blues9';
P1.Cgrm=Cgrm;

sPth=sprintf('%s/Synth%03dBnds',H.Path,Nbnds);
unix(sprintf('mkdir -p %s',sPth));

P=[orderfields(P1) orderfields(R) orderfields(S)];
for jp=1:length(P)
    Nm=P(jp).Name;
    y=[zeros(ceil(H.fs/5),1); P(jp).h; zeros(ceil(H.fs/5),1)];
    y=y/max(abs(y))*(1-1e-24);
    audiowrite(sprintf('%s/%s.wav',sPth,Nm),y,H.fs,'BitsPerSample',24);
end

%    CnvExcFrc(H,P(jp),linspace(0.1,10,3),'Impulse',sprintf('%s/ExcFrc%03dBnds',H.Path,Nbnds),Nm);
%    CnvExcFrc(H,P(jp),linspace(0.1,10,3),'Gauss',sprintf('%s/ExcFrc%03dBnds',H.Path,Nbnds),Nm);
