function Ds=GetTimit(sndpth,Shff)

fprintf('Reading source sentences... '); 
Ds=[];
for jd=1:8
	for jtree=1:2
		if jtree==1; trstrng='/test'; else trstrng='/train'; end
		Dtmp=dir([sndpth trstrng '/dr' num2str(jd) '/']);
		for jt=1:length(Dtmp)
			if (Dtmp(jt).isdir==1 && strcmp(Dtmp(jt).name(1),'.')~=1)
				D2tmp=dir([sndpth trstrng '/dr' num2str(jd) '/' Dtmp(jt).name '/*.wav']);
				jch=0;
				while jch<length(D2tmp) 
					jch=jch+1;
					if (strcmp(D2tmp(jch).name,'sa1.wav')==1 || strcmp(D2tmp(jch).name,'sa2.wav') || strcmp(D2tmp(jch).name(1:2),'sx')==1)
						D2tmp(jch)=[];
						jch=jch-1;
					else
						D2tmp(jch).name=[sndpth trstrng '/dr' num2str(jd) '/' Dtmp(jt).name '/' D2tmp(jch).name];
					end
				end
				Ds=[Ds; D2tmp];
			end
		end
	end
end

% shuffle if necessary
if strcmp(Shff,'shuffle')
    [~,srtndx]=sort(rand(size(Ds)));
    Ds=Ds(srtndx);
end
