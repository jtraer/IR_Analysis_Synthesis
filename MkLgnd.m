function MkLgnd(V,Dh,PltPrm);

for jj=1:length(V); hp=plot(-1,-1,V(jj).mrk); hold on
    set(hp,'color',V(jj).cmp);
    set(hp,'linewidth',3,'markersize',10);
    cnt=0;
    for jh=1:length(Dh);
        eval(sprintf('if strcmp(Dh(jh).%s,V(jj).name); cnt=cnt+1; end;',PltPrm));
    end
    V(jj).name=sprintf('%s: %dIRs',V(jj).name,cnt);
end; 
legend({V.name}); 
