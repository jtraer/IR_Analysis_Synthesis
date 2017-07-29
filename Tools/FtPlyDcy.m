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

function [bt,NsFlr,Test,FVE,Init]=FtPlyDcy(x,tt,NPly,Init0); 
if nargin<4; Init0=[]; end

% Iteratively smooth and save smoothed time series and downsampled time vectors
smx=x; smt=tt; cnt=1;   
smrc{1}=smx; ttrc{1}=tt;
while length(smx)>10; cnt=cnt+1;
    % downsample
    smt2=smt(1:2:end);
    smx=interp1(smt,smx,smt2);
    smt=smt2;
    % retain each version for model fitting 
    smrc{cnt}=smx;  ttrc{cnt}=smt; 
end

% to save time make sure smrc is not too large (though we still want the first and the last)
tmp=ceil(length(smrc)/10);
smrc2{1}=smrc{1};
ttrc2{1}=ttrc{1};
for jj=1:(floor(length(smrc)/tmp));
    smrc2{1+jj}=smrc{jj*tmp};
    ttrc2{1+jj}=ttrc{jj*tmp};
end
smrc2{length(smrc2)+1}=smrc{end};
ttrc2{length(ttrc2)+1}=ttrc{end};
smrc=smrc2;
ttrc=ttrc2;

% iterate through the different smoothed versions of the time series and fit an exponential fit model to each
ftcnt=0;
for jft=1:length(smrc); 
	% cycle through various estimates of the noise floor
	for jprc=[1 10:10:90 99]
        ftcnt=ftcnt+1; 
		% crudely compute the Noise floor (just a guess really)
		NsFlr0=prctile(20*log10(abs(smrc{jft})),jprc);
		% compute the index at which the signal falls less than the noise floor
		Ndx0=find(20*log10(abs(smrc{jft}))<=NsFlr0); Ndx0=min(Ndx0);
		if isempty(Ndx0); Ndx0=length(smrc{jft}); end; 
        % and the time of this index
        Tinf0=ttrc{jft}(Ndx0);
        % fit a linear decay
        ws = warning('off','all');  % Turn off warning
        Ft=polyfit(ttrc{jft}(1:gather(Ndx0)),20*log10(abs(smrc{jft}(1:gather(Ndx0)))),NPly);
        warning(ws)  % Turn it back on.
        % save the initializing parameters  
        fInit=[Tinf0 Ft]; 
        % use a least-mean squares fitting algorithm and save the LMS error to
        % that fit from the unsmoothed data (x)
        txx=20*log10(abs(smrc{jft}+1e-24));
        t3=ttrc{jft};	
        mnvl=[t3(2) -inf*ones(1,length(fInit)-1)]; mnvl(end-1)=-5e3;
        mxvl=[t3(end) inf*ones(1,length(fInit)-1)]; mxvl(end-1)=-1;
        best_fit{ftcnt}=fminsearchbnd(@(fit) PlyNsFlrFt(gather(txx),t3,fit), [gather(fInit)],mnvl,mxvl,optimset('MaxFunEvals',1e3,'MaxIter',1e3,'Display','off'));
        % now compute the actual error to the raw (UNSMOOTHED) data
        err(ftcnt)=PlyNsFlrFt(20*log10(abs(x)),tt,best_fit{ftcnt});
        % and if a starting value has been recommended try with that too
        if isempty(Init0)==0
            ftcnt=ftcnt+1; 
            best_fit{ftcnt}=fminsearchbnd(@(fit) PlyNsFlrFt(txx,t3,fit), [Init0], [-inf 0 t3(2)], [0 inf t3(end)],optimset('MaxFunEvals',itnm,'MaxIter',itnm,'Display','off')); 
            err(ftcnt)=asymp_error3(20*log10(abs(x)),tt,best_fit{ftcnt});
        end
    end
end
% find the smallest LMS fit of all the smoothed versions
[~,mndx]=min(err);
best_of_best_fits=best_fit{mndx};

% Extract the paramaters to output
% inflection point
Test=best_of_best_fits(1);   
if Test<tt(2); Test=tt(2); end; 
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
%figure(157); 
%% raw
%plot(tt,20*log10(abs(x))); hold on; 
%% smoothed
%smx=smrc{rcndx(mndx)}; smt=ttrc{rcndx(mndx)};
%plot(smt,20*log10(abs(smx)),'c')
%% area considered
%plot(tt(Lrc(rcndx(mndx)))*ones(1,2),[min(20*log10(abs(x))) max(20*log10(abs(x)))],'k:');
%% raw fit
%hp=plot([min(tt) max(tt)],ones(1,2)*NsFlr,'-.');
%set(hp,'color',[1-((length(bt)-2)/26) 0 0]);
%hp=plot([tt(1:infndx)],polyval(bt,tt(1:infndx)),'--');
%set(hp,'color',[1-((length(bt)-2)/26) 0 0]);
%hp=plot([Test],polyval(bt,Test),'rs');
%set(hp,'color',[1-((length(bt)-2)/26) 0 0]);
%axis tight; 
%set(gca,'xlim',[0 3*Test+1e-24],'ylim',[NsFlr-2 max([(bt(end)+5) (NsFlr+2)])]);
%set(gca,'xtick',[Test],'xticklabel',[num2str(round(Test*1e3))],'ytick',bt(end),'yticklabel',[num2str(round(bt(end)))])
%set(get(gca,'title'),'string','');
%title(sprintf('b=%2.1f,FVE=%d',bt(1),round(FVE*100)))
%drawnow;
