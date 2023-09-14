function M = make_diffusion_stencil(G)

num_nodes = length(G.Nodes.Coord(:,1));
M = zeros(num_nodes);
scf = 1e-6; %convert to meters, the coordinates in the .swc are in micrometers

for i=1:num_nodes
    pred = predecessors(G,i)';
    succ = successors(G,i)';
    nghb = [pred,succ];
    
    % first set the diagonal entries
    dl=G.Nodes.Coord(nghb,:)-G.Nodes.Coord(i,:); dl = dl.*scf;
    
    dl = dl.^2; dl = sum(dl,2); dl = dl.^0.5; total_dl = sum(dl);
    dl = dl.^(-1);
    M(i,i)  = -1*sum(dl);
    
    % now set the off diagonal entries
    for j=nghb
        dl=G.Nodes.Coord(j,:)-G.Nodes.Coord(i,:); dl = dl.*scf;
        dl = dl.^2; dl = sum(dl,2); dl = dl.^0.5; dl = dl.^(-1);
        M(i,j) = dl;
    end
    
    M(i,:) = M(i,:).*(2/total_dl);
end

%spy(M)

end