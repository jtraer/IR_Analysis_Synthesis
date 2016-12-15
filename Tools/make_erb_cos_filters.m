%generates filters with cosine frequency response functions on an erb-transformed frequency axis 
%sr=48000; %sampling rate
%N=30; %number of channels
%low_lim = 100; %low cutoff of lowest band
%hi_lim = 20000; %high cutoff of highest band

function [filts,Hz_cutoffs,freqs] = make_erb_cos_filters(signal_length, sr, N, low_lim, hi_lim)

if rem(signal_length,2)==0 %even length
    nfreqs = signal_length/2;%does not include DC
    max_freq = sr/2;
    freqs = [0:max_freq/nfreqs:max_freq]; %go all the way to nyquist
else %odd length
    nfreqs = (signal_length-1)/2;
    max_freq = sr*(signal_length-1)/2/signal_length; %max freq is just under nyquist
    freqs = [0:max_freq/nfreqs:max_freq];
end   
cos_filts = zeros(nfreqs+1,N);

if hi_lim>sr/2
    hi_lim = max_freq;
end
%make cutoffs evenly spaced on an erb scale
cutoffs = e2freq([freq2e(low_lim) : (freq2e(hi_lim)-freq2e(low_lim))/(N+1) : freq2e(hi_lim)]);

for k=1:N
    l = cutoffs(k);
    h = cutoffs(k+2); %adjacent filters overlap by 50%
    l_ind = min(find(freqs>l));
    h_ind = max(find(freqs<h));
    avg = (freq2e(l)+freq2e(h))/2;
    rnge = (freq2e(h)-freq2e(l));
    cos_filts(l_ind:h_ind,k) = cos((freq2e( freqs(l_ind:h_ind) ) - avg)/rnge*pi); %map cutoffs to -pi/2, pi/2 interval
end

%add lowpass and highpass to get perfect reconstruction
filts = zeros(nfreqs+1,N+2);
filts(:,2:N+1) = cos_filts;
h_ind = max(find(freqs<cutoffs(2))); %lowpass filter goes up to peak of first cos filter
filts(1:h_ind,1) = sqrt(1 - filts(1:h_ind,2).^2);
l_ind = min(find(freqs>cutoffs(N+1))); %highpass filter goes down to peak of last cos filter
filts(l_ind:nfreqs+1,N+2) = sqrt(1 - filts(l_ind:nfreqs+1,N+1).^2);

Hz_cutoffs = cutoffs;

%subplot(2,1,1); plot(freqs,sum(filts.^2,2))
%subplot(2,1,2); semilogx(freqs,sum(filts.^2,2))
