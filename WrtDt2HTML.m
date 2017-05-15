%* == function WrtDt2HTML.m : Takes stimuli/IRs stored in a structure BH and makes a set of plots to be arranged in an interactive html  
%** Written by James Traer - jtraer@mit.edu

function WrtDt2HTML(BH,fNm,PltPrms)

%* == Inputs
%** BH  : a structure of audio files 
%** fNm : a filename to save the data under
fntsz=15;
tmtpth='/Users/jtraer/LabBook/Workshop/PerceptualExperiments/Sandbox/timit'
if exist(tmtpth)==7
    Ds=GetTimit(tmtpth,'shuffle');
    Ds=Ds(1);
    for js=1:length(Ds);
        [TMT(js).s,TMT(js).fs]=audioread(Ds(js).name);
        audiowrite('tmt.wav',TMT(js).s,TMT(js).fs);
    end
else
    [TMT(1).s,TMT(1).fs]=audioread(tmt.wav);
end

%* == Preamble : make sure we have the requisite paths and search for flags
cfl=mfilename('fullpath'); %eval(sprintf('! ./Subroutines/ArchiveCode.sh %s.m',cfl));   
[stts,out]=unix(sprintf('grep -n HACK %s.m',cfl)); 
if length(out)>length(cfl)+153; 
    fprintf(2,'\n\nWARNING: check these hacks\n%s\n',out(length(cfl)+153:end)); 
    keyboard; 
end

%* == Log : Just started today (3rd-Oct)
%** TODO: we need better audio (the IR convolved with speech)
%** TODO: when debugged search and replace variable names for more generic ones

%* ==== Crunch ====

%** get the path to which we save files
sndx=[]; for jtmp=1:length(fNm); if strcmp(fNm(jtmp),'/'); sndx=[sndx jtmp]; end; end
if isempty(sndx)
    OtPth='';
else
    OtPth=fNm(1:sndx);
end

%** get the field names of the structure
flds=fieldnames(BH);
%** Open the JSON file for writing
fid2=fopen(sprintf('%s.json',fNm),'w')
fprintf(fid2,'var stimList = [\n')
%** => start loop over IRs (jIR)
for jIR=[1:length(BH)]
    %** =>  make a folder to save images and audio to
    FldrNm=sprintf('%s/h_%s_%dBnds',OtPth,BH(jIR).Name,length(BH(jIR).ff)-2);
    eval(sprintf('! mkdir -p %s',FldrNm));

    %** => add a comma to separate indices in the JSON
    if jIR~=1;
        fprintf(fid2,',\n');
    end
    %** => open up a new element in the JSON
    fprintf(fid2,'{\n"blah":\t"blah"'); % the blahs put in an empty field so that hereafter all fields can be written in the same format

    %** => Plot special values 
    %*** => This should be fixed when we want to generalize this code to other stimuli
    %*** => T0
    fprintf(fid2,',\n"T0":\t%2.3f',median(BH(jIR).RT60));
    %*** => path to audio
    %**** TODO update this to rescale volumes!!!
    Dh=dir(sprintf('%s/h_denoised_%03d.wav',BH(jIR).Path,length(BH(jIR).ff)));
    eval(sprintf('! cp %s/%s %s/h.wav',BH(jIR).Path,Dh(1).name,FldrNm));
    fprintf(fid2,',\n"sound":\t"%s/h.wav"',FldrNm);
    %*** => copy to a folder of just audio
    eval(sprintf('! cp %s/%s IRMAudio/Audio/%s.wav',BH(jIR).Path,Dh(1).name,BH(jIR).Name));

    %*** => convolve with a TIMIT sentence
    if length(Ds)>0;
        [y,fy]=RIRcnv(TMT(1).s,TMT(1).fs,BH(jIR).nh,BH(jIR).fs,1);
        audiowrite(sprintf('%s/tmt1.wav',FldrNm),y,fy);
        fprintf(fid2,',\n"speech":\t"%s/tmt1.wav"',FldrNm);
    end
    %%*** => write an image of the time series
    unix(sprintf('sips -s format png %s/IR.jpg --out %s/ts.png',BH(jIR).Path,FldrNm));
    fprintf(fid2,',\n"TimeSeries":\t"%s/ts.png"',FldrNm);

    %*** => write an image of the synthetic time series
    %figure;
    %tt=[1:length(BH(jIR).h_Eco)]/BH(jIR).fs;
    %plot(tt*1e3,sign(BH(jIR).h_Eco).*abs(BH(jIR).h_Eco).^(0.6));
    %hold on
    %plot(tt([2 end])*1e3,10^(-1*0.6)*ones(1,2),'k:');
    %plot(tt([2 end])*1e3,-10^(-1*0.6)*ones(1,2),'k:');
    %plot(tt([2 end])*1e3,10^(-2*0.6)*ones(1,2),'k:');
    %plot(tt([2 end])*1e3,-10^(-2*0.6)*ones(1,2),'k:');
    %plot(tt([2 end])*1e3,10^(-3*0.6)*ones(1,2),'k:');
    %plot(tt([2 end])*1e3,-10^(-3*0.6)*ones(1,2),'k:');
    %set(gca,'xlim',[-20 1.2e3],'xscale','log');
    %xlabel('Time (ms)'); ylabel('Amplitude (compressed)')
    %set(gca,'fontsize',fntsz);
    %drawnow
    %saveas(gcf,sprintf('%s/ts_Eco',FldrNm),'epsc');
    %saveas(gcf,sprintf('%s/ts_Eco.png',FldrNm));
    %fprintf(fid2,',\n"Eco_TimeSeries":\t"%s/ts_Eco.png"',FldrNm);
    %%**** => save the audio
    %audiowrite(sprintf('%s/h_Eco.wav',FldrNm),BH(jIR).h_Eco,BH(jIR).fs,'BitsPerSample',24);
    %fprintf(fid2,',\n"sound_Eco":\t"%s/h_Eco.wav"',FldrNm);
    %close all

    %%*** => plot C-gram
    unix(sprintf('sips -s format png %s/Cgram.jpg --out %s/Cgrm.png',BH(jIR).Path,FldrNm));
    fprintf(fid2,',\n"Cgrm":\t"%s/Cgrm.png"',FldrNm);
    
    %%*** => plot synthetic C-gram
    %figure;
    %h=BH(jIR).h_Eco.'; fs=BH(jIR).fs;
    %Npts=length(h);
    %[fltbnk,bff,erbff]=make_erb_cos_filters(3*Npts,fs,Nbnds,BH(jIR).ff(1),BH(jIR).ff(end));
    %Cgrm=generate_subbands([zeros(Npts,1); h; zeros(Npts,1)].',fltbnk);
    %Cgrm=Cgrm(Npts+[1:Npts],:).'; fprintf('Cgram made\n')
    %imagesc(BH(jIR).tt,[1:Nbnds+2],20*log10(abs(hilbert(Cgrm.').')));
    %axis xy
    %colorbar
    %set(gca,'clim',[-60 0])
    %cmp=othercolor('BuGn9',256);
    %colormap(cmp)
    %set(gca,'ytick',[2 Nbnds+2],'yticklabel',sprintf('%2.2f\n',BH(jIR).ff([2 Nbnds-1])/1e3));
    %set(gca,'xlim',[0 0.2],'ylim',[1 Nbnds+2]);
    %xlabel('Time (s)');
    %ylabel('Frequency (kHz)')
    %set(gca,'fontsize',fntsz);
    %drawnow
    %saveas(gcf,sprintf('%s/Cgrm_Eco',FldrNm),'epsc');
    %saveas(gcf,sprintf('%s/Cgrm_Eco.png',FldrNm));
    %fprintf(fid2,',\n"Cgrm_Eco":\t"%s/Cgrm_Eco.png"',FldrNm);
    %close all
    
    %%*** => plot RT60
    unix(sprintf('sips -s format png %s/RT60.jpg --out %s/RT60.png',BH(jIR).Path,FldrNm));
    fprintf(fid2,',\n"RT60":\t"%s/RT60.png"',FldrNm);
    
    %%*** => plot spectrum
    unix(sprintf('sips -s format png %s/IR_AttckSpc.jpg --out %s/Spc.png',BH(jIR).Path,FldrNm));
    fprintf(fid2,',\n"Spc":\t"%s/Spc.png"',FldrNm);
    
    %** ==> loop over fields in structure and write them to the JSON file (jfld)
    for jf=1:length(flds); 
        fld=flds{jf};
        vl=getfield(BH(jIR),fld);
        %*** ==> if it's a string write it
        if ischar(vl);
            fprintf(fid2,',\n"%s":\t"%s"',fld,vl);
        else
            if isstruct(vl)==0
                %*** ==> if it's a single number write it
                if length(vl)==1
                    fprintf(fid2,',\n"%s":\t%3.4f',fld,vl);
                end
            end
        end
    end
    %** ==< end-loop over fields in structure 

    %** => when all images are written close out the json structure
    fprintf(fid2,'\n}');
    %** => close all figures
    close all
end % for jIR
%** =< end-loop over folders of IRs (jIR)

% end the json file
fprintf(fid2,'\n]');
fclose(fid2)

%* == write the HTML file ==
%** Delete current lines
[~,LnNdx]=unix('sed -n ''/<img src="IRMAudio/='' IR_Data_Summary.html');
LnNdx=str2num(LnNdx);
for jln=1:length(LnNdx);
    unix(sprintf('sed -i.bak -e ''%dd'' IR_Data_Summary.html',LnNdx(1)));
end
%** Write new ones
[~,LnNdx]=unix('sed -n ''/<div id="Stats">/='' IR_Data_Summary.html');
LnNdx=str2num(LnNdx);
for jPlt=1:length(PltPrms);
    Dplt=dir(sprintf('IRMAudio/%s/*.png',PltPrms{jPlt}));
    for jp=1:length(Dplt);
        unix(sprintf('awk ''NR==%d{print "    <img src=\\"IRMAudio/%s/%s\\">"}7'' IR_Data_Summary.html >tmp.html',LnNdx+1,PltPrms{jPlt},Dplt(jp).name)); 
        unix('mv tmp.html IR_Data_Summary.html')
    end
end

%* == TODO: Save this code to a summary file
%eval(sprintf('! grep "%%\\*" %s.m > tmp.org',cfl))
%eval('! sed ''s/^ *//g'' < tmp.org > tmp2.org');  % remove whitespace
%eval('! sed ''s/^[ \\t]+//g'' < tmp2.org > tmp.org');  % remove tabs too
%eval(sprintf('! sed ''s/^.//'' tmp.org > %s.org',cfl))   % remove '%' so emacs can read the indenting
