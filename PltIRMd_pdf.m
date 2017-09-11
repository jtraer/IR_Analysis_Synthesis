function PltIRMd_pdf(H,C)

ff=linspace(H.ff(1),H.ff(end),100);
Rtt=linspace(0,3*max([H.Modes.RT60]),100);
Opp=linspace(min([H.Modes.OnPwr])-20,max([H.Modes.OnPwr])+20,100);
Rtt=linspace(0,3*max([H.Modes.RT60]),100);
Opp=linspace(min([H.Modes.OnPwr])-20,max([H.Modes.OnPwr])+20,100);

subplot(1,3,1);
pdf=normpdf(ff,H.MdStts.pdf_f.mu,H.MdStts.pdf_f.sigma);
hp=plot(pdf,ff/1e3,'-');
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on;
set(gca,'yscale','log','xscale','log');
set(gca,'xlim',[1e-2 3])
hold on
hold off;
axis tight;
xlabel('Gaussian fit');
ylabel('Mode Frequency (kHz)')

subplot(1,3,2);
pdf=normpdf(Rtt,H.MdStts.pdf_RT60.mu,H.MdStts.pdf_RT60.sigma);
hp=plot(Opp,pdf,'-');
pclr=get(hp,'color');
set(hp,'linewidth',3);
hold on;
axis tight;
set(gca,'xscale','log');
set(gca,'xlim',[1e-2 3]);
xlabel('Mode RT60 (s)');
ylabel('Gaussian fit')
title([H.Path ': Modes '])

subplot(1,3,3);
pdf=normpdf(Opp,H.MdStts.pdf_OP.mu,H.MdStts.pdf_OP.sigma);
hp=plot(Opp,pdf,'-');
pclr=get(hp,'color');
set(hp,'linewidth',3);
axis tight;
set(gca,'xlim',[-80 40]);
xlabel('Mode Onset Power (db)');
ylabel('Gaussian fit')
