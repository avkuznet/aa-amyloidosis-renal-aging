function savefig_emf(f,name)
% Export a figure at the size used for the figures of this paper.
set(f,'Units','inches','Position',[1 1 5.2 3.9], ...
      'PaperUnits','inches','PaperPosition',[0 0 5.2 3.9]);
set(findall(f,'-property','FontSize'),'FontSize',11);
print(f,name,'-dmeta');
end
