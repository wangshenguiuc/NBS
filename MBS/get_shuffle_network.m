function [ random_net ] = get_shuffle_network( adj_net,random_network_type )
%GET_SHUFFLE_random_net Summary of this function goes here
%   method=4: preserve degree for each node (strong randomness)
%   method=8: preserve degree distribution only (weak randomness)

ngene = size(adj_net,1);
adj_net(1:ngene+1:ngene*ngene) = 0;
if strcmp(random_network_type,'Random_network_preserve_degree_distribution')
    up_diag_adj_net = adj_net;
    for i=1:ngene
        up_diag_adj_net(i,i:ngene) = 0;
    end
    shuffle_gid = randperm(ngene);
    random_net = up_diag_adj_net(shuffle_gid,:);
    random_net = max(random_net,random_net');
elseif strcmp(random_network_type,'Random_network_preserve_degree')
    old_deg = sum(adj_net);
    random_net = adj_net;
    [vi,vj] = find(random_net);
    nedge = length(vi);
    nshuffle = nedge*10;
    
    for i=1:nshuffle        
        if mod(i,1000)==0
            [vi,vj] = find(random_net);
            nold_edge = nedge;
            nedge = length(vi);
            assert(nedge==nold_edge);
            fprintf('%f finished\n',i/nshuffle);
        end
        pair = randsample(nedge,2);
        x1 = vi(pair(1));
        x2 = vi(pair(2));
        y1 = vj(pair(1));
        y2 = vj(pair(2));
        if x1==y2 || x2==y1
            continue
        end
        if random_net(x1,y1)> 0 && random_net(x2,y2) > 0 && random_net(x1,y2) == 0 && random_net(x2,y1) == 0
            random_net(x1,y2) = random_net(x1,y1);
            random_net(x2,y1) = random_net(x2,y2);
            random_net(x1,y1) = 0;
            random_net(x2,y2) = 0;
            
            random_net(y2,x1) = random_net(y1,x1);
            random_net(y1,x2) = random_net(y2,x2);
            random_net(y1,x1) = 0;
            random_net(y2,x2) = 0;
        end
    end
    new_deg = sum(random_net);
    assert(isequal(old_deg,new_deg));
else
    error('wrong method');
end

end

