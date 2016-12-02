function H=hExtrct(rc,G);
%* Takes a waveform of recorded audio and a structure G of a golay sequence and properties and extracts the environmental IR.

%* Match lengths of audio
%** Find numberof completed Golay cycles and truncate accordingly
Nrps=floor(length(rc)/G.Lg);
rc=rc(1:Nrps*G.Lg);

%* Split recorded audio into sections
rcnt=0;
for jrp=1:Nrps
	sc=rc((jrp-1)*2*G.Lg+[1:(2*G.Lg)]);
	
	%** Truncate large peaks in recorded audio

	%* Split into sections
	scA=sc(1:G.Lg);
	scB=sc(G.Lg+[1:G.Lg]);

	%* Extract IR snapshots
             %[th,shft]=FreqDomainProcessing7(g,ab(:,jch),fs,NGly,SpkrSpc,Nbnds); 
	% Something with scA and G.a etc
end

%* Measure variability

%* Create a structure



