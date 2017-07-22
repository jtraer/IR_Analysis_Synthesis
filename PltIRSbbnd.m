function BdBndsFlg=PltIRSbbnd(H,V,jbn,tmp,tmp2,tmp3,sb_fs,Pft,Test,sdDRR,sdB);

BdBndsFlg=0;
% Plot
plot([1:length(tmp)]/H.fs,20*log10(abs(tmp))); 
hold on
plot([1:length(tmp2)]/H.fs,20*log10(abs(tmp2)),'c');
plot([0:(length(tmp3)-1)]/sb_fs,20*log10(abs(tmp3)),'g:');
plot([1:length(tmp)]/H.fs,Pft(2)+Pft(1)*[1:length(tmp)]/H.fs,'r--');
plot([1:length(tmp)]/H.fs,(Pft(2)+sdDRR/2)+(Pft(1)-sdB/2)*[1:length(tmp)]/H.fs,'r:');
plot([1:length(tmp)]/H.fs,(Pft(2)-sdDRR/2)+(Pft(1)+sdB/2)*[1:length(tmp)]/H.fs,'r:');
plot([1 length(tmp)]/H.fs,Pft(2)+Pft(1)*Test*ones(1,2),'k--');
if ~isempty(V)
    plot([1:length(tmp)]/H.fs,tmp3(1)-(60/V.RT60(jbn))*[1:length(tmp)]/H.fs,'m--');
    %** Check if recorded decay is less than or equal to the speaker-microphone IR
    if V.RT60(jbn)>(-60/Pft(1))*0.75;
        BdBndsFlg=1;
        text(0.5*Test,Pft(2),1.001,sprintf('Danger: Speaker-Microphone IR RT60 is %d%% of recorded',round(100*V.RT60(jbn)/(-60/Pft(1)))));
    end
end
hold off
title(sprintf('%s: Band %d',H.Path,jbn));
set(gca,'xlim',[0 3*Test]);
set(gca,'ylim',[Pft(2)+Pft(1)*Test-20 Pft(2)+20]);
xlabel('Time (s)')
ylabel('Power (dB)')
