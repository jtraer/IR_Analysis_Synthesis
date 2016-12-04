function H=hExtrct(rc,G,Pth);
%* == hExtrt.m i.e. IR extract ==
%* Takes a waveform of recorded audio and a structure G of a golay sequence and properties and extracts the environmental IR.
%** If no path is specified make it an empty string
if nargin<3; Pth=''; end

%* Match lengths of audio
%** Find numberof completed Golay cycles and truncate accordingly
Nrps=floor(length(rc)/G.Ng/2);
rc=rc(1:Nrps*G.Ng*2);

%* Split recorded audio into sections
rcnt=0;
for jrp=1:Nrps
	sc=rc((jrp-1)*2*G.Ng+[1:(2*G.Ng)]);
	
	%** Truncate large peaks in recorded audio
    %*** IRs recorded in absence of noise have kurtosis ~3 (empirically determined).  Higher values are likely due to extraneous noise -- which can be more sparse.
    cnt=0; 
    while kurtosis(sc)>5; 
        cnt=cnt+1;
        %*** => find peaks (probably due to noise)
        thrsh=prctile(abs(sc),99.99-1e-4*cnt); % gradually increase the number of points to be shrunk
        I=find(abs(sc)>thrsh); 
        %*** => decrease magnitude of peaks ito the median value so that the noise is reduced 
        sc(I)=sign(sc(I)).*rand(size(I))*prctile(abs(sc),50);  
    end
    %** save sections to reconstitute a de-noised time-series
    rc_dn((jrp-1)*2*G.Ng+[1:2*G.Ng])=sc;

	%* Split into complementary sections
	scA=sc(1:G.Ng);
	scB=sc(G.Ng+[1:G.Ng]);

	%* Extract IR snapshots
    %** Fourier transform
    A=fft(G.a(:),G.Ng);    
    Ap=fft(scA(:),G.Ng);
    B=fft(G.b(:),G.Ng);    
    Bp=fft(scB(:),G.Ng);
	%** Cross correlate (if the recorded signals [Ap,Bp] are conjugated the IR is flipped backwards in time)
    AA=Ap.*conj(A);
    BB=Bp.*conj(B);
    aa=ifft(AA);
    bb=ifft(BB);
    %** Combine sequences to eliminate self-noise
    h_snp(:,jrp)=aa+bb;
end
fprintf('%s: raw IR extracted.  Saving...\n',Pth)

%* Average over snapshots and measure variability
H.h=mean(h_snp,2);
H.h_var=std(h_snp,[],2);
H.h_snps=h_snp;
H.fs=G.fs;

%* Estimate the region that likely contains signal (beginning and end are usually noise)
%** Plot and query user for the start and end times of the IR
figure; 
plot(20*log10(abs(H.h))); 
hold on;
plot(20*log10(abs(H.h+sign(H.h).*H.h_var)),':'); 
plot(20*log10(abs(H.h-sign(H.h).*H.h_var)),':'); 
title('Manually extract start and end of IR');
drawnow;
M=GtMtDt(sprintf('%s/IR_StartEnd.txt',Pth),{'Start_index';'End_index'});
eval(sprintf('xlm(1)=%s;',M.Start_index))
eval(sprintf('xlm(2)=%s;',M.End_index))
%** trim out just the IR from the background noise
H.h=H.h([xlm(1):xlm(2)]);
H.h_var=H.h_var([xlm(1):xlm(2)]);
H.h_snps=H.h_snps([xlm(1):xlm(2)],:);

%* Save metadata
H.DateCreated=date;
H.No_of_snapshots=Nrps;
H.Golay_Code=G.Name;

%* Plot and save (if a path was given)
if length(Pth)>0;
    close all
    %** => plot recording
    figure(1)
    plot([1:length(rc)]/G.fs,rc);
    %** => plot de-noised recording
    hold on
    plot([1:length(rc)]/G.fs,rc_dn);
    xlabel('Time (s)');
    ylabel('Waveform amplitude')
    title([Pth ': Recording'])
    %saveas(gcf,sprintf('%s/Recording',Pth),'fig'); 
    saveas(gcf,sprintf('%s/Recording',Pth),'jpg');

    %** => plot IR
    figure(2)
    %*** => compress the time series for plotting
    h=sign(H.h).*abs(H.h).^(0.6);
    h_var1=H.h+sign(H.h).*H.h_var;
    h_var1=sign(h_var1).*abs(h_var1).^(0.6); 
    h_var2=H.h-sign(H.h).*H.h_var;
    h_var2=sign(h_var2).*abs(h_var2).^(0.6); 
    %*** Plot
    plot([1:length(H.h)]/H.fs,h);
    %*** => plot variablility
    hold on
    plot([1:length(H.h)]/H.fs,h_var1,':');
    plot([1:length(H.h)]/H.fs,h_var2,':');
    xlabel('Time (s)');
    ylabel('Waveform amplitude (compressed)')
    set(gca,'xscale','log')
    title([Pth ': Raw IR and variability'])
    %saveas(gcf,sprintf('%s/Raw_IR',Pth),'fig');
    saveas(gcf,sprintf('%s/Raw_IR',Pth),'jpg');
   
    %** => plot snapshots
    figure(3)
    %*** => scroll through snaphsots
    for jsnp=1:Nrps
        h=H.h_snps(:,jsnp);
        %*** => compress the time series for plotting
        h=sign(h).*abs(h).^(0.6);
        %*** Plot
        plot([1:length(h)]/H.fs,h);
        hold on
    end; hold off
    xlabel('Time (s)');
    ylabel('Waveform amplitude (compressed)')
    set(gca,'xscale','log')
    title([Pth ': Raw IR snapshots'])
    %saveas(gcf,sprintf('%s/Raw_IR_Snapshots',Pth),'fig');
    saveas(gcf,sprintf('%s/Raw_IR_Snapshots',Pth),'jpg');
end

