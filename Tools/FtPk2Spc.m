function [err,Pkk]=FtPk2Spc(spc,ff,Pk)

Pk_f=Pk(1);
Pk_a=Pk(2);
Pk_wd=Pk(3);

Pk=zeros(length(ff),1);
bnwd=ceil(Pk_wd/(ff(2)-ff(1)));
[~,Cbn]=min(abs(ff-Pk_f));
mn=max(Cbn-bnwd,1);
mx=min(Cbn+bnwd,length(ff));
Pk([mn:mx])=Pk_a*sin(linspace(0,pi,length([mn:mx])));

err=abs(spc)-Pk;
err=rms(err);
Pkk=Pk;
