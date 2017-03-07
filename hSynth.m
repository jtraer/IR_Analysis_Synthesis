function H=hSynth(H,ftp);
%* == hSynth.m i.e. Make synthetics ==
%* Takes a structure of IR properties (i.e. H, output by hPrp.m), and makes a family of synthetics

fprintf('We should add something here to manimpualte modes');

set(0,'DefaultFigureVisible','off');

rcnt=0;
rcnt=rcnt+1; R(rcnt).Name='R_SpcCnst';    R(rcnt).cmp='BuPu3';
rcnt=rcnt+1; R(rcnt).Name='R_NoMod';      R(rcnt).cmp='YlGn3';
rcnt=rcnt+1; R(rcnt).Name='R_Lin';        R(rcnt).cmp='Reds3';

scnt=0;
scnt=scnt+1; S(scnt).Name='S_Exp';        S(scnt).cmp='Greens9';
scnt=scnt+1; S(scnt).Name='S_LinEnd';     S(scnt).cmp='Reds9';
scnt=scnt+1; S(scnt).Name='S_LinStrt';    S(scnt).cmp='Oranges9';
scnt=scnt+1; S(scnt).Name='S_SpcCnst';    S(scnt).cmp='BuPu9';

fntsz=15;

%* == Compute kurtosis in 10ms windows ==

%* == Compute the cochleagram ==
%** zeropad to avoid edge effects
Npts=length(H.nh);
Nbnds=length(H.ff);
[fltbnk,ff,erbff]=make_erb_cos_filters(3*Npts,H.fs,Nbnds,H.ff(1),H.ff(end));
Cgrm=generate_subbands([zeros(Npts,1); H.nh; zeros(Npts,1)].',fltbnk);
Cgrm=Cgrm(Npts+[1:Npts],:).'; 
%** Remove the extreme bands
Cgrm=Cgrm([2:(end-1)],:);
H.ff=ff([2:(end-1)]);

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
        %** Remove Modulations
        elseif strcmp('R_NoMod',R(jr).Name)
            Dcy=20*log10(abs(tmp2));
            DesDcy=H.DRR(jbn)-60/H.RT60(jbn)*tt;
            R(jr).Cgrm(jbn,:)=tmp.*10.^((DesDcy-Dcy)/20);
        %** impose linear decay
        elseif strcmp('R_Lin',R(jr).Name)
            Dcy=20*log10(abs(tmp2));
            Nln=ceil(H.RT60(jbn)*(60+H.DRR(jbn))/60*H.fs);
            if Nln>Npts; Nln=Npts; end
            DesDcy=linspace(10^(H.DRR(jbn)/20),10^(-60/20),Nln);    
            DesDcy=20*log10(DesDcy);
            DesDcy=[DesDcy -100*ones(1,Npts-Nln)];
            R(jr).Cgrm(jbn,:)=tmp.*10.^((DesDcy-Dcy)/20);
        end
    end
end

% resynthesize the remixed IRs
for jr=1:length(R)
    rmxCgrm=R(jr).Cgrm;
    rmxCgrm=[zeros(1,size(rmxCgrm,2)); rmxCgrm; zeros(1,size(rmxCgrm,2))];
    rmxh=collapse_subbands([zeros(size(rmxCgrm)) rmxCgrm zeros(size(rmxCgrm))].',fltbnk);
    rmxh=rmxh(Npts+[1:Npts]);
    rmxh=rmxh/max(abs(rmxh));
    R(jr).h=rmxh;
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
    end
    S(js).h=h_s(:);
    S(js).Cgrm=C_s;
end

%* == Plot ==
P1.h=H.nh;
P1.Name='RW';
P1.cmp='Blues9';
P1.Cgrm=Cgrm;

P=[orderfields(P1) orderfields(R) orderfields(S)];
for jp=1:length(P)
    Nm=P(jp).Name;
    CnvExcFrc(H,P(jp),linspace(0.1,10,3),'Impulse',sprintf('%s/ExcFrc%03dBnds',H.Path,Nbnds),Nm);
    CnvExcFrc(H,P(jp),linspace(0.1,10,3),'Gauss',sprintf('%s/ExcFrc%03dBnds',H.Path,Nbnds),Nm);
end 
