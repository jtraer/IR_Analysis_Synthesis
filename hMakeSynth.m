% This function generates synthetic IRs from the following input parameters
% -- aalph  : the subband dependent Wet-Dry ratios
% -- bbt    : the subband dependent decay rates 
% -- Flmz   : the high and low frequency limits
% -- Nt     : the lenghth of the output synthetic in samples
% -- fs     : the sampling frequency
% -- dbthrsh: the cutoff treshold above which we normalize the IRs to have the same power
% -- ERlm   : the limit (in ms) before which we consider signal early reeflections and hence we do not consider peaks in this region as contributing to the magnitude of the IR distortion
% -- DcyMd  : the mode of decay

% possible decay modes include
%%% -- Exp  : matches exponentials as in HCgrm2Snd2.m
%%% -- Lin1 : matches a linear with the same starting value and same RT60 as
%%% the exponential (ie. too much energy)
%%% -- Lin2 : matches the same RT60 and energy as the exponential (ie alpha is very small)
%%% -- Lin3 : matches the same alpha and area (ie. RT60 is now very small)
%%% -- TTxx :  matches an exponential but then truncates the tail to remove xx% energy

% === Written by James Traer, jtraer@mit.edu, March 2015

function [ns,C]=hMakeSynth(aalph,bbt,Flmz,Nt,fs,dbthrsh,ERlm,DcyMd); 

% switch debugger on or off
dbgFLG=0; fcnt=0;
% generate a time vector and get the index at the Eearl-reflection/diffuse tail boundary
tt=[1:Nt]/fs;
[~,ERndx]=min(abs(tt-ERlm/1e3));
% generate pink noise
Nbnds=length(aalph)-2;
[fltbnk,bff,erbff]=make_erb_cos_filters(3*Nt,fs,Nbnds,Flmz(1),Flmz(2));
nz=randn(Nt,1);
NZgrm=generate_subbands([zeros(Nt,1); nz; zeros(Nt,1)].',fltbnk);
NZgrm=NZgrm(Nt+[1:Nt],:).'; 
%ensure the power is pink (in envelope space!!)  and 0dB
nSbRMS=rms(abs(hilbert(NZgrm.'))); 
NZgrm=NZgrm./(nSbRMS.'*ones(1,Nt));
C=NZgrm;
% DEBUGGER: plot the Noise Cochleagram 
if dbgFLG==1; fcnt=fcnt+1; figure(500+fcnt); imagesc(20*log10(abs(C))); axis xy; colormap('gray'); colorbar; title('MakeSynthIR.m: The baseline "pink" noise'); end
% scroll through subbands and impose the required decay envelope
for jf=1:Nbnds+2
    % extract the required subband
    tmp1=C(jf,:); 
    % compute the time at which this band becomes inaudible
    tj=(-dbthrsh+aalph(jf))/bbt(jf); if isinf(tj)==1; tj=max(tt); end
    % iteratively enforce the canonical exponential decay until the DRR is measured as desired
    err=1e6; desired_a=aalph(jf); a_we_impose=aalph(jf);
    tmpExp=tmp1.*10.^((a_we_impose-bbt(jf)*tt)/20); 
    while abs(err)>1e10; 
        % enforce decay
        tmpExp=tmp1.*10.^((a_we_impose-bbt(jf)*tt)/20); 
        % measure the properties
        tmp=[zeros(size(tmpExp)) tmpExp zeros(size(tmpExp))];
        tmp=(abs(hilbert(tmp)));
        tmp=tmp([length(tmpExp)+[1:length(tmpExp)]]);
        [Pft,~,~,~]=FtPlyDcy(tmp,tt,1,1);
        measured_a=Pft(2); 
        % measure the discrepancy between desired and measured
        err=measured_a-desired_a; fprintf('a_we_impose=%2.2f, err=%2.2f\n',a_we_impose,err)
        % implement a correction
        a_we_impose=a_we_impose-err;
    end
    % find the peaks in the diffuse tail (t>ERlm) at which this this subband is audible
    tmp2=20*log10(abs(tmpExp)); I=find(tmp2>dbthrsh); 
    I(find(I<ERndx))=[]; 
    % sum the energy contained in these peaks as a measure of the subband audible power
    AExp=sum(-dbthrsh+tmp2(I))/fs;

    % define the subband threshold we will use for matching
    thrshj=dbthrsh;
    % if the start of the reverb is below the audible threshold then this subband is likely inaudible.  Anything we will do is thus arbitrary and so why not set an arbitrary matching threshold
    if (AExp<=0||a_we_impose<=dbthrsh); 
        % redefine the subband threshold
        thrshj=min([max(tmp2(ERndx:end)) a_we_impose-ERlm/1e3*bbt(jf)])-10; 
        I=find(tmp2>thrshj);
        I(find(I<ERndx))=[];
        % redefine the area to be matched
        AExp=sum(-thrshj+tmp2(I))/fs;
        % and redefine the time to be considered
        tj=(-thrshj+a_we_impose)/bbt(jf);
    end

    % DEBUGGER: plot the channel decays
    if dbgFLG==1; if rem(jf,ceil(Nbnds/4))==0; fcnt=fcnt+1; figure(500+fcnt); plot(tt,20*log10(abs(tmpExp))); hold on; plot(tt,(aalph(jf)-bbt(jf)*tt),'r--'); plot(ones(2,1)*tt(ERndx),[thrshj 0],'m:'); title(sprintf('MakeSynthIR.m: Subband %d/%d',jf,(Nbnds+2))); set(gca,'ylim',[thrshj 0]); hold off; end; end
        
    % if the Decay Mode is Exponential we use the value already defined
    if strcmp(DcyMd,'Exp')==1   
        tmp1=tmpExp;
    % linear decay with the same starting amplitude and decay time (thus more energy)
    elseif strcmp(DcyMd,'Lin1')==1   
        [~,ndx]=min(abs(tt-tj));
        % compute the adjustment
        tmp1(1:ndx)=tmp1(1:ndx).*(10^(a_we_impose/20)*(1-[0:1:(ndx-1)]/ndx)); tmp1(ndx+1:end)=0;
        % record error
        %eerr(jf)=max([abs(tj-ot)/tj abs(a_we_impose-oa)/a_we_impose]);
        if dbgFLG==1; figure(601); subplot(2,1,1); plot(tt,20*log10(abs(tmp1)),'k-.'); hold on; plot(tt,aalph(jf)-bbt(jf)*tt,'c--');  hold off; axis([0 2 -80 0]); drawnow; end

    % linear decay with matched energy and decay time (thus lower starting amplitude)
    elseif strcmp(DcyMd,'Lin2')==1  
        % find the index beyond which the IR is inaudible
        %[~,ndx]=min(abs(tt-tj)); 
        I=find(20*log10(abs(tmpExp))>thrshj); ndx=max(I);
        % If the index to be matched is 0 just set the subband and don't iterate
        if ndx<=ERndx; tmp1=zeros(size(tmp1));
            a2=a_we_impose;
            tmphs=tmp1;
            tmp1(1:ndx)=linspace(10^(a2/20),10^(thrshj/20),ndx).*tmphs(1:ndx); tmp1(ndx+1:end)=0;
        %...otherwise compute something interesting
        else
            % set the DRR to the ecological value and set upper and lower limits for the DRR value
            a2=a_we_impose; ulm=thrshj+2*(a2-thrshj)+2; llm=thrshj;
            % save the subband time structure
            tmphs=tmp1;
            % iterate through values of subband DRR (a2) until the error is small
            err=1e6; itcnt=0;
            while (ulm-llm)>1e-1; itcnt=itcnt+1; 
                % define a linearly decaying time series
                tmp1(1:ndx)=linspace(10^(a2/20),10^(thrshj/20),ndx).*tmphs(1:ndx); tmp1(ndx+1:end)=0;
                % transfer to log-domain
                tmp2=20*log10(abs(tmp1((1:ndx)))); 
                % compute the audible power
                I=find(tmp2>thrshj);
                I(find(I<ERndx))=[];
                A2=sum(-thrshj+tmp2(I))/fs;
                % compute the audible power difference between this and an exponential  made with these parameters
                err=(A2-AExp)/AExp;
                % respecify the upper and lower limits to the DRR (a2) accordingly
                if err>0; ulm=a2; else llm=a2; end; 
                % in case we are trapped in a local minimum start again
                if (rem(itcnt,1e2)==0); ulm=0; llm=thrshj; end
                % recalculate a2 midway between the limits
                a2=(llm+ulm-2*thrshj)/2+thrshj;

                if itcnt>1e2; keyboard; end

                % DEBUGGER
                if dbgFLG==1; fprintf('Lin2: subband %d/%d: aalph=%2.2f; a2=%2.2f, AExp=%2.2f; A2=%2.2f; Aerr=%2.4f\n',jf,Nbnds,aalph(jf),a2,AExp,A2,err); end
            end
        end
    % linear decay that matches energy and starting amplitude (thus a shorter RT60)  
    elseif strcmp(DcyMd,'Lin3')==1
        t2=tj; ulm=length(tmp1)/fs; llm=ERlm/1e3;
        % save the subband time structure
        tmphs=tmp1;
        % iterate through values of t2 until the error is as small as can be
        err=1e6; itcnt=0;
        while abs(ulm-llm)>1/fs; itcnt=itcnt+1; 
            % define a linearly decaying time series
            [~,ndx]=min(abs(tt-t2)); 
            tmp1(1:ndx)=linspace(10^(a_we_impose/20),10^(thrshj/20),ndx).*tmphs(1:ndx); tmp1(ndx+1:end)=0;
            % transfer to log-domain
            tmp2=20*log10(abs(tmp1((1:ndx)))); 
            % compute the audible power
            I=find(tmp2>thrshj); 
            I(find(I<ERndx))=[]; 
            A2=sum(-thrshj+tmp2(I))/fs;
            % compute the audible power difference between this and an exponential  made with these parameters
            err=(A2-AExp)/AExp;
            % respecify the upper and lower limits to the DRR (a2) accordingly
            if err>0; ulm=t2; else llm=t2; end; 
            % recalculate t2 midway between the limits
            t2=(llm+ulm)/2;
            % DEBUGGER
            if dbgFLG==1; fprintf('Lin3: subband %d/%d: aalph=%2.2f; t2=%2.2f, AExp=%2.2f; A2=%2.2f; Aerr=%2.4f\n',jf,Nbnds,aalph(jf),t2,AExp,A2,err); end
        end
        %t2=tj; ulm=2*tj; llm=0;
        %tmphs=tmp1;
        %while abs(ulm-llm)>1/fs;
        %    [~,ndx]=min(abs(tt-t2));
        %    tmp1(1:ndx)=linspace(10^(a_we_impose/20),10^(thrshj/20),ndx); tmp1(ndx+1:end)=0;
        %    tmp2=20*log10(abs(tmp1((1:ndx)))); I=find(tmp2>thrshj);
        %    A2=sum(-thrshj+tmp2(I))/fs;
        %    err=A2-AExp;
        %    if err>0; ulm=t2; else llm=t2; end; 
        %    t2=(llm+ulm)/2;
        %    tmp1(1:ndx)=tmp1(1:ndx).*tmphs(1:ndx);
        %end
        %if dbgFLG==1; figure(601); subplot(2,1,1); plot(tt,20*log10(abs(tmp1)),'k-.'); hold on; plot(tt(1:ndx),tmp2,'r--'); plot(tt,aalph(jf)-bbt(jf)*tt,'c--');  hold off; axis([0 2 -80 0]); drawnow; end   

    % linear decay initialized by equality but then scaled downwards
    % uniformly in amplitude until the area is equal to that of the
    % exponential
    %elseif strcmp(DcyMd,'Lin4')==1
    %    [~,ndx]=min(abs(tt-tj));
    %    tmphs=tmp1;
    %    eps=10^(aalph(jf)/100); err=1e6; ulm=dbthrsh; llm=0;
    %    while abs(err)>Aj*1e-3;
    %        tmp1(1:ndx)=10^(aalph(jf)/20-eps/20)*(1-[0:1:(ndx-1)]/ndx); tmp1(ndx+1:end)=0;
    %        tmp2=20*log10(tmp1(1:ndx)); I=find(tmp2>-dbthrsh);      
    %        A2=sum(dbthrsh+tmp2(I))/fs;
    %        err=A2-Aj;
    %        if err>0; llm=eps; else ulm=eps; end; if 20*log10(abs(eps))>=dbthrsh; err=0; end
    %        eps=(llm+ulm)/2; if abs(ulm-llm)<1/fs; err=0; end
    %        tmp1(1:ndx)=tmp1(1:ndx).*tmphs(1:ndx);
    %    end
    %    if DBsw==1; figure(601); subplot(2,1,1); plot(tt,20*log10(abs(tmp1)),'k-.'); hold on; plot(tt(1:ndx),tmp2,'r--'); plot(tt,aalph(jf)-bbt(jf)*tt,'c--');  hold off; axis([0 2 -80 0]); drawnow; end   


    % % synthesize an IR with energy that oscillates sinusoidally about
    % % exponential decay
    % elseif strcmp(DcyMd(1:3),'Sin')==1  
    %     if strcmp(DcyMd(4),'P')==1; s=1; elseif strcmp(DcyMd(4),'M')==1; s=-1; end
    %     % number of cycles
    %     Nc=str2double(DcyMd(5:end));
    %     % compute amplitude of the sinusoid
    %     % -- period
    %     w=2*pi/(tj/Nc);
    %     % -- amplitude such that at a quarter cycle there has been no drop
    %     a=bbt(jf)*(tj/Nc)/4;
    %     % impose
    %     tmp1=tmp1.*10.^((aalph(jf)-bbt(jf)*tt+s*a*sin(w*tt))/20); 
%   %      tmp1=tmp1.*(10^(aalph(jf)/20)*exp(-bbt(jf)*tt+s*a*sin(w*tt)));

    % % a guassian window with symmetric attacks and decays and the direct
    % % arrival in the middle
    % elseif strcmp(DcyMd,'Gauss')==1        
    %    % combute the gradient of decay
    %    b=3/tj^2*(3*log(10)+log(aalph(jf)));
    %    % and the additive factor
    %    a=-3*log(10)+b*tj^2/4;
    %    % impose
    %    [~,ndx]=min(abs(tt-tj));
    %    tmp1(end/2-floor(ndx/2)+[1:ndx])=exp(a-b*(tt(1:ndx)-tj/2).^2);
    %    tmp1(end/2)=1;

    %% a symmetric window that exponentially rises and then falls symmetrically with the direct arroval in the middle    
    %elseif strcmp(DcyMd,'SymmExp')==1
    %    [~,ndx]=min(abs(tt-tj));
    %    ndx2=floor(length(tmp1)/2)-1;
    %    tmp1(end/2+[1:ndx2])=aalph(jf)*exp(-2*bbt(jf)*tt(1:ndx2));
    %    tmp1(end/2+[-ndx2:-1])=aalph(jf)*exp(-2*bbt(jf)*fliplr(tt(1:ndx2)));
    %    tmp1(end/2)=1;

    %% a box car normalized to have the same energy
    %elseif strcmp(DcyMd,'Bxcr')==1
    %    tbx=0.75*tj;
    %    a=tj*(dbthrsh+aalph(jf))/2;
    %    a=10^((a-dbthrsh)/20)/tbx;
    %    [~,ndx]=min(abs(tt-tbx));
    %    tmp1([1:ndx])=a*tmp1([1:ndx]).*ones(1,ndx);
	%tmp1([ndx+1:end])=0;
    %    tmp1(1)=1;

    % % multiple exponential decays of different gradients at different times  
    % elseif strcmp(DcyMd(1:3),'Mlt')==1
    %     fct=2;
    %     [~,ndx]=min(abs(tt-tj));
    %     Nsg=str2double(DcyMd(7:end));
    %     dtj=floor(ndx/Nsg); 
    %     Alph=aalph(jf);
    %     for jsg=1:(Nsg-1)
    %         sndx=(jsg-1)*dtj; if sndx==0; sndx=1; end
    %         tmp1(sndx+[1:dtj])=Alph*exp(-bbt(jf)*fct^(Nsg/2-jsg)*tt(1:dtj).');
    %         Alph=tmp1(sndx+dtj);
    %     end
    %     tmp1(sndx+dtj+1:end)=Alph*exp(-bbt(jf)*fct^(-Nsg/2)*tt(1:length(tmp1(sndx+dtj+1:end))).');    
    end  
    % add term to cochleagram
    C(jf,:)=tmp1(1:Nt);
    % DEBUGGER: plot the channel decays
    if dbgFLG==1; if rem(jf,ceil(Nbnds/4))==0; fcnt=fcnt+1; figure(500+fcnt); plot(tt,20*log10(abs(tmp1))); hold on; plot(tt,20*log10(abs(tmpExp)),'c--'); plot(tt,(aalph(jf)-bbt(jf)*tt),'r--'); plot(ones(2,1)*tt(ERndx),[thrshj 0],'m:'); title(sprintf('MakeSynthIR.m: Subband %d/%d for a %s IR',jf,(Nbnds+2),DcyMd)); set(gca,'ylim',[thrshj 0]); hold off; end; end
end
% DEBUGGER: plot the cochleagram
if dbgFLG==1; fcnt=fcnt+1; figure(500+fcnt); imagesc(20*log10(abs(C))); axis xy; colormap('gray'); colorbar; set(gca,'clim',[dbthrsh max(aalph)]); set(gca,'xlim',[0 2*max(ceil(((60+aalph)./bbt)*fs))]); title('MakeSynthIR.m: Cochleagram after the envelopes applied'); end

% zeropad
C=[zeros(size(C)) C zeros(size(C))];
% transfer back to time domain
ns=collapse_subbands(C.',fltbnk);
ns=ns(Nt+[1:Nt]); ns=ns.'; 
% final touches
ns(1)=1; ns=ns/max(abs(ns));
C=C(:,Nt+[1:Nt]);

%keyboard

% plot the time series
if dbgFLG==1; fcnt=fcnt+1; figure(500+fcnt); plot(tt,ns); set(gca,'xlim',[0 2*max(ceil(((60+aalph)./bbt)))]); title(sprintf('MakeSynthIR.m: Time series for %s IR',DcyMd)); end
