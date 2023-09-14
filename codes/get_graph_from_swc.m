function G=get_graph_from_swc(filename)
%% this function takes a .swc file and generates a digraph of the geometry
% G = is a digraph with edges, nodes, and attached information such as
% radius, edgelength, subset type, and coordinates

use_radius= 1 ;
[~,id,pid,coord,r,subset]=readSWC(filename);
if ~isempty(r<0)
    r =abs(r);      % occasionally and .swc has a negative radius value
end

if (use_radius == 0)
    % use coordinate distance to compute radius
    crdist = (coord(:,1).^2 + coord(:,2).^2 + coord(:,3).^2).^0.5;
    s = sort(crdist);
    crdist(crdist==0) = s(2);
    markerSize = (1./(crdist))*max(r)*10;
else
    % radius information
    markerSize = r; %(r./max(r)).*2;
end

% s,t are for making the edges
t = id(2:end); s = pid(2:end);

% this creates a graph must be digraph! inorder to track parent to child
G = digraph(s,t);

%% assigns information to the nodes and edges of the graph
G.Nodes.Size=markerSize;
G.Nodes.type = subset;
G.Edges.width = markerSize(G.Edges.EndNodes(:,1));
G.Nodes.Coord=coord;

end