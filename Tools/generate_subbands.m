function subbands = generate_subbands(sample_sound, filts)

if size(sample_sound,1)==1 %turn into column vector
    sample_sound = sample_sound';
end
N=size(filts,2)-2;
signal_length = length(sample_sound);
filt_length = size(filts,1);
fft_sample = fft(sample_sound);
if rem(signal_length,2)==0 %even length - 
    fft_filts = [filts' fliplr(filts(2:filt_length-1,:)')]'; %generate negative frequencies in right place; filters are column vectors
else %odd length
    fft_filts = [filts' fliplr(filts(2:filt_length,:)')]';
end
fft_subbands = fft_filts.*(fft_sample*ones(1,N+2));%multiply by array of column replicas of fft_sample
subbands = real(ifft(fft_subbands)); %ifft works on columns; imag part is small, probably discretization error?
