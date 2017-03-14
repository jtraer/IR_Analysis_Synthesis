function PltIRPrps(H);

ps=get(gcf,'position');
set(gcf,'position',[ps(1:2) ps(3)*3 ps(4)/2]);
%* Plot Time series
T0=median(H.RT60);
tlm=[-T0*0.1 T0]*1e3;
subplot(1,6,1)
plot([1:length(H.nh)]/H.fs*1e3,sign(H.nh).*abs(H.MaxAmp*H.nh).^(0.6));
set(gca,'xlim',tlm);
hold on;
text(0.1*T0,max(abs(H.MaxAmp*H.nh))^0.6,1.001,sprintf('Peak amp=%2.2f',H.MaxAmp));
% TODO print the maximim amplitude
hold off;
set(gca,'xlim',[0 600]);
xlabel('Time (ms)');
ylabel('Amp (compressed)');
title(sprintf('%s',H.Name));

%* Plot Kurtosis
subplot(1,6,2)
plot([1:length(H.nh)]/H.fs*1e3,H.krt);
set(gca,'xlim',tlm);
hold on;
plot([1 length(H.nh)]/H.fs*1e3,3*ones(1,2),'k--');
plot(H.Tgs*1e3,3,'k+')
hold off;
set(gca,'yscale','log')
set(gca,'xlim',[0 100]);
xlabel('Time (ms)');
ylabel('Kurtosis');
title(sprintf('Kurtosis'));

%* Plot ER spectrum
subplot(1,6,3)
plot(H.spcER,H.ff/1e3);
hold on;
hold off;
set(gca,'yscale','log')
set(gca,'xlim',[-10 10]);
xlabel('Power (dB)');
ylabel('Freq (kHz)');
title(sprintf('ER Spectrum'));

%* Plot tail spectrum
subplot(1,6,4)
plot(H.spcGs,H.ff/1e3);
hold on;
hold off;
set(gca,'yscale','log')
set(gca,'xlim',[-1 1]);
xlabel('Power (dB)');
ylabel('Freq (kHz)');
title(sprintf('Tail Spectrum'));

%* Plot DRR
subplot(1,6,5)
plot(-H.DRR,H.ff/1e3);
hold on;
hold off;
set(gca,'yscale','log')
set(gca,'xlim',[-15 15]);
xlabel('DRR (dB)');
ylabel('Freq (kHz)');
title(sprintf('DRR'));

%* Plot RT60
subplot(1,6,6)
plot(H.RT60,H.ff/1e3);
hold on;
hold off;
set(gca,'yscale','log');
set(gca,'xlim',[0 2]);
xlabel('RT60 (s)');
ylabel('Freq (kHz)');
title(sprintf('RT60'));
drawnow
