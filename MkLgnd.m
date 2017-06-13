function MkLgnd(V);

for jj=1:length(V); hp=plot(-1,-1,V(jj).mrk); hold on
    set(hp,'color',V(jj).cmp);
    set(hp,'linewidth',3,'markersize',10);
end; 
legend({V.name}); 
