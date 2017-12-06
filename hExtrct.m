function H=hExtrct(rc,fr,G,Pth);
%* == hExtrt.m i.e. IR extract ==
%* Takes a waveform of recorded audio and a structure G of a golay sequence and properties and extracts the environmental IR.
%** If no path is specified make it an empty string
if nargin<4; Pth=''; end

set(0,'DefaultFigureVisible','off');

G
isempty(G)

if ~isempty(G)
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
      cnt=0; fprintf('Rep %d: Removing large peaks...',jrp)
      while kurtosis(sc)>5; 
          cnt=cnt+1;
          %*** => find peaks (probably due to noise)
          thrsh=prctile(abs(sc),99.99-1e-4*cnt); % gradually increase the number of points to be shrunk
          I=find(abs(sc)>thrsh); 
          %*** => decrease magnitude of peaks ito the median value so that the noise is reduced 
          sc(I)=sign(sc(I)).*rand(size(I))*prctile(abs(sc),50);  
      end; fprintf('done\n')
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
else
  % if G is empty then the raw recoridng is an IR
  h_snp=rc(:);
  Nrps=1;
  rc_dn=[];
end

%* Average over snapshots and measure variability
H.h=mean(h_snp,2);
H.h_var=std(h_snp,[],2);
H.h_snps=h_snp;
H.fs=fr;

%* Estimate the region that likely contains signal (beginning and end are usually noise)
%** Plot and query user for the start and end times of the IR
if exist(sprintf('%s/IR_StartEnd.txt',Pth))~=2; 
    set(0,'DefaultFigureVisible','on');
end
figure;
subplot(3,1,1)
plot(20*log10(abs(H.h))); 
hold on;
plot(20*log10(abs(H.h+sign(H.h).*H.h_var)),':'); 
plot(20*log10(abs(H.h-sign(H.h).*H.h_var)),':');
axis tight
ylm=get(gca,'ylim');
for jj=1:4;
    x=(jj/340)*H.fs;
    plot([1 1]*x,ylm,'k--'); % 1m increments
    text(x,ylm(1)+0.1*(ylm(2)-ylm(1)),sprintf('%dm path',jj));
end
title('Manually extract start and end of IR');
set(gca,'xlim',[0 1000]);
subplot(3,1,2)
plot(20*log10(abs(H.h))); 
hold on;
plot(20*log10(abs(H.h+sign(H.h).*H.h_var)),':'); 
plot(20*log10(abs(H.h-sign(H.h).*H.h_var)),':'); 
axis tight
ylm=get(gca,'ylim');
for jj=1:3;
    x=(10*jj/340)*H.fs;
    plot([1 1]*x,ylm,'k--'); % 1m increments
    text(x,ylm(1)+0.1*(ylm(2)-ylm(1)),sprintf('%dm path',10*jj));
end
set(gca,'xlim',[0 5000]);
subplot(3,1,3)
plot(20*log10(abs(H.h))); 
hold on;
plot(20*log10(abs(H.h+sign(H.h).*H.h_var)),':'); 
plot(20*log10(abs(H.h-sign(H.h).*H.h_var)),':'); 
axis tight
ylm=get(gca,'ylim');
for jj=1:4;
    x=(20*jj/340)*H.fs;
    plot([1 1]*x,ylm,'k--'); % 1m increments
    text(x,ylm(1)+0.1*(ylm(2)-ylm(1)),sprintf('%dm path',20*jj));
end
set(gca,'xlim',[0 25000]);
drawnow; 
M=GtMtDt(sprintf('%s/IR_StartEnd.txt',Pth),{'Start_index';'End_index'});
eval(sprintf('xlm(1)=%s;',M.Start_index))
eval(sprintf('xlm(2)=%s;',M.End_index))
%** trim out just the IR from the background noise
H.h=H.h([xlm(1):xlm(2)]);
H.h_var=H.h_var([xlm(1):xlm(2)]);
H.h_snps=H.h_snps([xlm(1):xlm(2)],:);

%* Measure spectra
nft=2^ceil(log2(length(H.h)));
spc=fft(H.h,nft);
H.spc=spc(1:end/2);
H.Spcff=[1:nft/2]*H.fs/nft;

%* Manually identify Mode frequencies
%if exist(sprintf('%s/IR_Modes.txt',Pth))~=2; 
%    fid=fopen(sprintf('%s/IR_Modes.txt',Pth),'w');
%    set(0,'DefaultFigureVisible','on');
%    figure; 
%    plot(20*log10(abs(H.spc)),H.Spcff);
%    axis tight; ylm=get(gca,'ylim');
%    set(gca,'ylim',[-ylm(2)/5 ylm(2)]);
%    hold on; plot([min(20*log10(abs(H.spc))) max(20*log10(abs(H.spc)))],zeros(1,2),'k--');
%    title('Manually extract modes: click under dashed line when all are selected');
%    drawnow; 
%    Mf=1;
%    mcnt=0;
%    while Mf>0; mcnt=mcnt+1;
%        Mf=ginput(1);
%        Mf=Mf(2);
%        if Mf>0;
%            [~,ndx]=min(abs(H.Spcff-Mf));
%            if ndx>50;
%                [~,ndx2]=max(abs(H.spc(ndx+[-50:50])))
%                Mf=H.Spcff(ndx+ndx2-51);
%            else 
%                [~,ndx2]=max(abs(H.spc([1:(ndx+50)])))
%                Mf=H.Spcff(ndx2);
%            end
%            Mdf(mcnt)=Mf;
%            [~,ndx]=min(abs(H.Spcff-Mf));
%            plot(20*log10(abs(H.spc(ndx))),H.Spcff(ndx),'ro');
%            drawnow
%            fprintf(fid,'Md%02d\t%f\n',mcnt,round(Mf));
%        end
%    end
%    fclose(fid)
%end
%[MdNdx,Mdf]=textread(sprintf('%s/IR_Modes.txt',Pth),['%s %f']);
%for jj=1:length(Mdf);
%    M=GtMtDt(sprintf('%s/IR_Modes.txt',Pth),{sprintf('Md%02d',jj)});
%    eval(sprintf('Modes(jj).cf=str2num(M.Md%02d);',jj,jj));
%end
%if isempty(Mdf);
%    Modes=[];
%end
%H.Modes=Modes;

%* Save metadata
H.DateCreated=date;
H.No_of_snapshots=Nrps;
if isempty(G);
  G.Name='';
  G.fs=fr;
end
H.Golay_Code=G.Name;

%* Plot and save (if a path was given)
set(0,'DefaultFigureVisible','off');
if length(Pth)>0;
    close all
    %** => plot recording
    figure(1)
    plot([1:length(rc)]/G.fs,rc);
    %** => plot de-noised recording
    hold on
    if ~isempty(rc_dn)
      plot([1:length(rc)]/G.fs,rc_dn);
    end 
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
    %*** => plot markers for trimming
    if length(h)>=100
        plot([10 30 100]/H.fs,h([10 30 100]),'d');
    else
        plot([10 30]/H.fs,h([10 30]),'d');
    end
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
    
    %** => plot spectrum
    figure(4)
    %*** => scroll through snaphsots
    plot(20*log10(abs(H.spc)),H.Spcff/1e3);
    hold on;
    xlabel('Power (db)');
    ylabel('Frequency (kHz)')
    set(gca,'yscale','log')
    title([Pth ': Raw IR spectra'])
    %saveas(gcf,sprintf('%s/Raw_IR_Snapshots',Pth),'fig');
    saveas(gcf,sprintf('%s/Raw_IR_Spc',Pth),'jpg');
end

