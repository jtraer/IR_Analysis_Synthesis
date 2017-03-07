% FtPlyDcy.m : Fit Polynomial Decay : Fits a two-segment model to an input time series in log-space under
% the assumption that the beginning section is decaying and
% the later section is measurement noise and hence is flat
% This computes the parameters of polynomial decay (in log space), amplitude, RT60 and Fraction of Variance Explained by the model  

%Inputs:
%-- x       : the time series
%-- tt      : a vector of time values (must be the same size as x)
%-- NPly    : the degree of the polynomial to be fit
%-- Qswtch  : if set to 1, all time consuming steps are neglected to speed up processing time.  This is used if this file will need to be used many times in a for loop to get approximate answers
%-- Init0   : (optional) a vector of initialization variables to be used in the
%               LMS fitting process (starting value [dB],decay rate
%               [dB/s],Inflection point where data meets noise [s])

% James Traer - jtraer@mit.edu - May 2014

function [bt,NsFlr,Test,FVE,Init]=FtPlyDcy(x,tt,NPly,Qswtch,Init0); 
if Qswtch==0; DBsw=1; itnm=1e4; else; DBsw=0; itnm=1e2; end
if nargin<5; Init0=[]; end

% Iteratively smooth all oscillations to yield a monotonic function 
% after each smoothing save a time series and a downsampled time vector
smx=x; smt=tt; xmns=ones(1,2*length(x)); cnt=0;   
while (isempty(xmns)~=1); cnt=cnt+1;
    % retain each version for model fitting 
    smrc{cnt}=smx;  ttrc{cnt}=smt; Lrc(cnt)=length(smt);
    % find neighboring peaks in the time series
    dx=sign(diff(smx));   
    dx=sign(diff(dx)); dx(1)=0;
    xmns=find(dx==1); % these are the local minima
    omn=1;  smx2=zeros(size(xmns));    smt2=zeros(size(xmns));
    if gpuDeviceCount~=0; smx2=gpuArray(smx2); smt2=gpuArray(smt2); end 
    % scroll through the peaks in the time series and median filter 
    % construct a new time series with each oscillatory cycle replaced by a single point
    for jmn=1:length(xmns)
        smx2(jmn)=median(smx(gather(omn):gather(xmns(jmn))));
        smt2(jmn)=mean(smt(gather(omn):gather(xmns(jmn))));
        omn=xmns(jmn);
    end
    % as long as there are peaks
    if isempty(xmns)~=1
        % pad the start and end
        smx2=[(smx(1)+smx2(1))/2 smx2 0]; smt2=[smt(1) smt2 smt(end)];
        % make sure we don't have any repeated time samples
        [smt2,I]=unique(smt2); smx2=smx2(I);
        % interpolate back to the original sampling - this will be used in thenext cycle
        [~,endx]=min(abs(tt-max(smt2)));
        smt=tt(1:(gather(endx)-1));
        smx=interp1(smt2,smx2,smt,'linear'); 
        % smooth
        %smx2=medfilt2(smx,[1 ceil(xmns(1)/8)],'symmetric');
    else % if there are no oscillations the function is monotonic
        smx=smx(isnan(smx)==0);
        smt=smt(1:length(smx));
    end
end
hold off

% iterate through the different smoothed versions of the time series and fit an exponential fit model to each
ftcnt=0;
for jft=1:length(smrc); 
	% cycle through various estimates of the noise floor
	for jprc=[10:10:90]
		% crudely compute the Noise floor (just a guess really)
		NsFlr0=prctile(20*log10(abs(smrc{jft})),jprc);
		% compute the index at which the signal falls less than the noise floor
		Ndx0=find(20*log10(abs(smrc{jft}))<=NsFlr0); Ndx0=min(Ndx0);
		if isempty(Ndx0); Ndx0=length(smrc{jft}); end; 
        % and the time of this index
        Tinf0=ttrc{jft}(Ndx0);
        % only proceed if the noise floor removes most late fluctuations
        Itl=find(tt>Tinf0);
        PkTl=find(20*log10(abs(x(Itl)))>NsFlr0);
        length(PkTl)/length(Itl);
        if length(PkTl)/length(Itl)<1; 
            % fit a linear decay
            ws = warning('off','all');  % Turn off warning
            Ft=polyfit(ttrc{jft}(1:gather(Ndx0)),20*log10(abs(smrc{jft}(1:gather(Ndx0)))),NPly);
            warning(ws)  % Turn it back on.
            % save the initializing parameters  
            fInit=[Tinf0 Ft]; 
            % use a least-mean squares fitting algorithm and save the LMS error to
            % that fit from the unsmoothed data (x)
            ftcnt=ftcnt+1; rcndx(ftcnt)=jft; 
            txx=20*log10(abs(smrc{jft}(1:Lrc(jft))));
            t3=ttrc{jft}(1:Lrc(jft));	
            t3=t3(~isnan(txx)); txx=txx(~isnan(txx));
            mnvl=[t3(2) -inf*ones(1,length(fInit)-1)]; mnvl(end-1)=-5e3;
            mxvl=[t3(end) inf*ones(1,length(fInit)-1)]; mxvl(end-1)=-1;
            best_fit{ftcnt}=fminsearchbnd(@(fit) PlyNsFlrFt(gather(txx),t3,fit), [gather(fInit)],mnvl,mxvl,optimset('MaxFunEvals',itnm,'MaxIter',itnm,'Display','off'));
            % now compute the actual error to the raw (UNSMOOTHED) data
            err(ftcnt)=PlyNsFlrFt(20*log10(abs(x)),tt,best_fit{ftcnt});
            % and if a starting value has been recommended try with that too
            if isempty(Init0)==0
                fprintf('Double check that FtPltDcy.m is doing what you want.  This seems dubious -- I thought we had migrated everything to an arbitrary polynomial fitter but here we have a linear model.  It is in an if statement that is probably defunct but if yous ee this message printed... it is being used\n.')
                ftcnt=ftcnt+1; rcndx(ftcnt)=jft; 
                best_fit{ftcnt}=fminsearchbnd(@(fit) TwoLnFt(txx,t3,fit), [Init0], [-inf 0 t3(2)], [0 inf t3(end)],optimset('MaxFunEvals',itnm,'MaxIter',itnm,'Display','off')); 
                %Enndx(ftcnt)=EndNdx; 
                err(ftcnt)=asymp_error3(20*log10(abs(x)),tt,best_fit{ftcnt});
            end
        end
    end
end
% find the smallest LMS fit of all the smoothed versions
[~,mndx]=min(err);
best_of_best_fits=best_fit{mndx};
Lbst=Lrc(rcndx(mndx));

% Extract the paramaters to output
% inflection point
Test=best_of_best_fits(1);   if Test<tt(2); Test=tt(2); end; 
[~,infndx]=min(abs(tt-Test));
% decay parameters
bt=best_of_best_fits(2:end);
%Compute the noisefloor
NsFlr=polyval(bt,Test);

% isolate the section before the inflection point
tt3=tt(1:infndx); xx3=x(1:infndx);
% compute the Fraction of Variance Explained by linear vs. non-linear fits
%FVE=GtGdnssLnFt(xx3,tt3);
dbx=20*log10(abs(xx3));
FVE=1-err(mndx)/rms(dbx-mean(dbx)); 
%FVE=1-rms(dbx-(bt(2)*tt3+bt(2)))/rms(dbx-mean(dbx)); 
%FVE=1-err(mndx)/rms(dbx); 

% ====== Plot ======
if DBsw==1;
     figure(157); 
     % raw
     plot(tt,20*log10(abs(x))); hold on; 
     % smoothed
     smx=smrc{rcndx(mndx)}; smt=ttrc{rcndx(mndx)};
     plot(smt,20*log10(abs(smx)),'c')
     % area considered
     plot(tt(Lrc(rcndx(mndx)))*ones(1,2),[min(20*log10(abs(x))) max(20*log10(abs(x)))],'k:');
     % raw fit
     hp=plot([min(tt) max(tt)],ones(1,2)*NsFlr,'-.');
     set(hp,'color',[1-((length(bt)-2)/26) 0 0]);
     hp=plot([tt(1:infndx)],polyval(bt,tt(1:infndx)),'--');
     set(hp,'color',[1-((length(bt)-2)/26) 0 0]);
     hp=plot([Test],polyval(bt,Test),'rs');
     set(hp,'color',[1-((length(bt)-2)/26) 0 0]);
     axis tight; 
     set(gca,'xlim',[0 3*Test+1e-24],'ylim',[NsFlr-2 max([(bt(end)+5) (NsFlr+2)])]);
     set(gca,'xtick',[Test],'xticklabel',[num2str(round(Test*1e3))],'ytick',bt(end),'yticklabel',[num2str(round(bt(end)))])
     set(get(gca,'title'),'string','');
     title(sprintf('b=%2.1f,FVE=%d',bt(1),round(FVE*100)))
     drawnow;
end
