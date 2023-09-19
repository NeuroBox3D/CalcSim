function n_inact=make_synapse_dist(filename,num_syn)
addpath(genpath([pwd, filesep, 'codes' ]));

dt = 10e-6;

if (nargin == 0)
    filename = 'morphos/beam64um/beamref5.swc';
end

G = get_graph_from_swc(filename);
n = length(G.Nodes.Size);

%num_syn = 1000;
ta =0.020; tb = 0.025;   % range of start times
synamps = [3.0,5.0].*1e-8;
a1 = synamps(1); a2 = synamps(2);    % range of amplitudes

synapses = make_random_synapses2(num_syn,n,ta,tb,a1,a2,dt);
save('synapses0.mat','synapses','-v7.3');

n_inact=num_syn-[900:-100:100,90:-10:10];
for i=1:length(n_inact)
    load(sprintf('synapses%i.mat',0),'synapses');
    synapses = inactivate_syn(synapses,n_inact(i));
    save(sprintf('synapses%i.mat',i),'synapses','-v7.3');
end
end

function synapses = make_random_synapses2(num_syn,n,ta,tb,a1,a2,dt)
t0 = ta:dt:tb;

for j=1:num_syn
    synapses(j).at_node = randi(n,1,1);
    synapses(j).start_time = t0(randi(length(t0),1,1));
    synapses(j).end_time = synapses(j).start_time + 0.01;
    synapses(j).amp = ((a2-a1).*rand(1,1)+a1);
end
end

function synapses = inactivate_syn(syn,n)
    active_ind = [];
    for i=1:length(syn)
        if (syn(i).amp ~= 0)
            active_ind = [active_ind i];
        end
    end
    
    for i=1:n
        syn(active_ind(i)).amp = 0;
    end
    synapses = syn;
end
