function K=IRKrt(H,Tbn);
%* IRKrt.m - function to measure kurtosis of an IR (in structure H) in time bins (Tbn) in milliseconds
%* James Traer - jtraer@mit.edu - Jul-2017

%** Find the numbr of points in the bin and make sure it is an even number
Nbn=ceil(Tbn/1e3*H.fs); 
if Nbn>length(H.h)*2;
    Nbn=length(H.h)/4;
end
Nbn=Nbn+rem(Nbn,2); 
%** Repeat the first and last sections
tmp=[H.h(1:Nbn/2); H.h; H.h(end-Nbn/2+1:end)]; 
%** Scroll through data points
krt=zeros(length(H.h),1); 
stndx=0; Nstp=1; 
while stndx<(length(tmp)/Nstp-Nbn); stndx=stndx+1; 
    sc=tmp((stndx-1)*Nstp+[1:Nbn]); 
    krt(stndx,:)=kurtosis(sc); 
end; 
K.krt=krt;
%** Compute the expected variance of kurtosis for samples of Gaussian noise  
VrKrt=24*Nbn*(Nbn-1)^2/((Nbn-3)*(Nbn-2)*(Nbn+3)*(Nbn+5)); 
K.VrKrt=VrKrt;
%** Classify data points as "Sparse" or "Noise-like" 
Sndx=find(krt>3+2*VrKrt);  %Sparse
Nndx=find(krt<=3+2*VrKrt); %Noise-like
%** Compute the crossover to Gaussian statistics as the point at which there has been as many Gaussian points as sparse (this is a stable measure but it is also arbitrary and crude)
%*** Find the maximum (preumably this is near the first arrival)
[~,mxndx]=max(abs(H.h));
Sndx(find(Sndx<=mxndx))=[];
Nndx(find(Nndx<=mxndx))=[];
NGs=0; NER=1; cnt=mxndx; 
while (NGs<=NER&&cnt<length(krt)); cnt=cnt+1; 
    NGs=length(find(Nndx<=cnt)); 
    NER=length(find(Sndx<=cnt)); 
end; 
K.Tgs=cnt/H.fs;
K.Sndx=Sndx;
K.Nndx=Nndx;
