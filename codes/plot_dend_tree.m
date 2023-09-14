function plot_dend_tree(G)
  figure
  set(gcf,'Visible','on')
  set(gcf,'units','normalized','position',[0 0 0.75 0.85])
  subplot(1,2,1)
  hold on
  pbaspect([1 1 1])
  ng=plot(G);
  xlabel('{\mu}m');
  ylabel('{\mu}m');
  ng.XData = G.Nodes.Coord(:,1);
  ng.YData = G.Nodes.Coord(:,2);
  ng.ZData = G.Nodes.Coord(:,3);
  
  ng.NodeLabel={};
  ng.LineWidth = G.Edges.width; 
  
  ng.MarkerSize = abs(G.Nodes.Size);
  ng.EdgeAlpha = 1.0; %0.25;
  ng.ArrowSize = 0;
  box on
  set(gca,'linewidth',3)
  
  subplot(1,2,2)
  M = adjacency(G);
  spy(M+M');
  title('Connectivity (Hines)');
  set(gca,'linewidth',3)
  set(gcf,'defaultAxesFontSize',30)
end

