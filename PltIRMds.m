function PltIRMds(H,C)

subplot(2,3,1);
hp=plot([H.Modes.RT60],[H.Modes.cf]/1e3,'o');
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on;
for jj=1:length(H.Modes);
    hp=plot(H.Modes(jj).RT60*ones(1,2),(H.Modes(jj).cf+[-1 1]/2*H.Modes(jj).bw)/1e3,'-');
end
set(gca,'yscale','log','xscale','log');
set(gca,'xlim',[1e-2 3])
hold on
hold off;
axis tight;
xlabel('Mode RT60 (s)');
ylabel('Mode Frequency (kHz)')

subplot(2,3,2);
hp=plot([H.Modes.OnPwr],[H.Modes.cf]/1e3,'o');
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on;
for jj=1:length(H.Modes);
    hp=plot(H.Modes(jj).OnPwr*ones(1,2),(H.Modes(jj).cf+[-1 1]/2*H.Modes(jj).bw)/1e3,'-');
end
set(gca,'yscale','log');
set(gca,'xlim',[-60 20])
hold on
hold off;
axis tight;
xlabel('Mode Onset power (dB)');
ylabel('Mode Frequency (kHz)')
title([H.Path ': Modes '])

subplot(2,3,3);
hp=plot([H.Modes.RT60],[H.Modes.OnPwr],'o');
pclr=get(hp,'color');
set(hp,'linewidth',3);
set(gca,'xscale','log');
set(gca,'xlim',[0 3])
set(gca,'ylim',[-60 0])
hold on
hold off;
axis tight;
xlabel('Mode RT60 (s)');
ylabel('Onset Power (dB)')

subplot(2,2,3);
%[hst,hstx]=histcounts([H.Modes.cf]/1e3,10);
%hstx=hstx(1:end-1)+(hstx(2)-hstx(1))/2;
%hp=plot(hstx,hst);
histogram([H.Modes.cf]/1e3,length(H.Modes));
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on
hold off;
axis tight;
xlabel('Mode Frequencies (kHz)');
ylabel('No of modes')

subplot(2,2,4);
df=[];
for jj=2:length(H.Modes);
    df=[df ([H.Modes(jj:end).cf]-H.Modes(jj-1).cf)];
end
%[hst,hstx]=histcounts([H.Modes.cf]/1e3,10);
%hstx=hstx(1:end-1)+(hstx(2)-hstx(1))/2;
%hp=plot(hstx,hst);
histogram([df]/1e3,2*length(H.Modes));
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on
hold off;
axis tight;
xlabel('Mode Frequency differences (kHz)');
ylabel('No of modes')

