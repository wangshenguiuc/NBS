addpath 'evaluation'
addpath 'permute_matrix'
addpath 'util'
addpath('matlab_submine')
addpath('evaluation')
addpath('graph_search')
addpath('find_subnet')
addpath('permute_matrix')
addpath('util')
addpath('NBS')
addpath('Hotnet')
warning off
rng(8);
ngene = 8;
npat = 7;
adj_network = zeros(ngene,ngene); % gene by gene ME network
p2g_c = cell(1,2);
net_c = cell(1,2);
%% ME example
mat = zeros(npat,ngene);
mat(1,1) = 1;
mat(2,2) = 1;
mat(3,3) = 1;
mat(4,4) = 1;
%
% mat(7,5) = 1;
% mat(7,6) = 1;
mat(7,6) = 1;
mat(6,6) = 1;
mat(5,8) = 1;

p2g_c{1} = mat;
%% co occur example
mat = zeros(npat,ngene);
mat(1:2,1:4) = 1;
% mat([1,2,3],2) = 1;
% mat([1,2,3,4],3) = 1;
% mat([2,4],4) = 1;

% mat(7,5) = 1;
mat(7,6) = 1;
mat(6,6) = 1;
mat(5,8) = 1;


p2g_c{2} = mat;
adj_network = zeros(ngene,ngene);
adj_network(1,2) = 1;
adj_network(2,4) = 1;
adj_network(1,3) = 1;
adj_network(4,3) = 1;
adj_network(1,4) = 1;
adj_network(2,3) = 1;

adj_network(5,6) = 1;
adj_network(6,7) = 1;
adj_network(7,8) = 1;
adj_network(4,5) = 1;
adj_network(5,7) = 1;

adj_network = max(adj_network,adj_network');
net_c{1} = adj_network;

adj_network = zeros(ngene,ngene);
adj_network(1,2) = 1;
adj_network(2,4) = 1;
adj_network(1,3) = 1;
adj_network(4,3) = 1;
adj_network(1,4) = 1;
adj_network(2,3) = 1;

adj_network = max(adj_network,adj_network');
net_c{2} = adj_network;
pat_diff_all = cell(1,2);
for nid = 1:2
    adj_network = net_c{nid};
    for pid =1:2%length(p2g_c)
        p2g = p2g_c{pid};
        npat = size(p2g,1);
        % for i=1:5
        %         adj_network(i,3) = 1;
        % end
        % adj_network(1,2) = 1;
        % adj_network(3,8) = 1;
        
        
        pat_diff = run_diffusion(adj_network,0.3,p2g); % run random walk
        for i=1:npat
            pat_diff(i,:) = pat_diff(i,:)*sum(p2g(i,:));
        end
        p2g = p2g(npat:-1:1,:);
        pat_diff = pat_diff(npat:-1:1,:);
        pat_diff(isnan(pat_diff)) = 0;
        HMobj = HeatMap(adj_network);
        addTitle(HMobj, 'network');
        pat_diff_all{pid} = [pat_diff_all{pid},pat_diff];
        
    end
    
end
HMobj = HeatMap(pat_diff_all{1},'Colormap',redblue(100));
plot(HMobj);
addTitle(HMobj, 'after propogate');
print(['../output/figure1_simulated_data/ME',num2str(pid),'.pdf'],'-dpdf');

HMobj = HeatMap(pat_diff_all{2},'Colormap',redblue(100));
plot(HMobj);
addTitle(HMobj, 'before propogate');
print(['../output/figure1_simulated_data/cooccur',num2str(pid),'.pdf'],'-dpdf');

% HMobj = HeatMap(adj_network');
% addTitle(HMobj, 'network')
% pat_diff = run_diffusion(adj_network,0.5,p2g);

