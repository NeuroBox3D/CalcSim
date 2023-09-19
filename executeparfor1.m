
morphfile = 'morphos/NMO_35893/ref2.swc';
num_syn = 1000;
n_inact=make_synapse_dist(morphfile,num_syn);

parfor i=0:length(n_inact)
tic
    calcium_simulationSBDF(morphfile,sprintf('vtk_%i',i),sprintf('data_%i',i),i,0);
toc
end
