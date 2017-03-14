function PltIRKrt(H,C,VrKrt,ftp);

plot([1:length(H.krt)]/H.fs,H.krt);
hold on
plot([1:length(H.krt)]/H.fs,(3+2*VrKrt)*ones(size(H.krt)),'k--');
plot([1:length(H.krt)]/H.fs,(3-2*VrKrt)*ones(size(H.krt)),'k--');
plot(H.Tgs,3,'ko');
if ~isempty(C)
    plot([1:length(C(1).krt)]/C(1).fs,C(1).krt,'k:');
end
hold off;
set(gca,'yscale','log')
xlabel('Time (s)');
ylabel('kurtosis')
title([H.Path ': Kurtosis'])
