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
network_file = '../data/network/InBio-Map_Entrez.sif';
% network_file = '..\data\network\FI_id.txt';
gene_name_mapping_file = '../data/util/gene_name_mapping.txt';
result_dir = '../output/subtype/';

dataset = 'TCGA';
if strcmp(dataset,'TCGA')
    survival_dir = '../data/survival/raw_survival/';
    filter_no_survival_pat = true;
    min_mutation = 10;
    ME_hop_l = 1:3;
elseif  strcmp(dataset,'ICGC')
    survival_dir = '../data/survival/ICGC/';
    filter_no_survival_pat = false;
    min_mutation = 10;
elseif  strcmp(dataset,'msk')
    survival_dir = '../data/survival/msk/';
    filter_no_survival_pat = false;
    min_mutation = 1;
    ME_hop_l = 0;
elseif  strcmp(dataset,'Xenograft')
    survival_dir = '../data/survival/Xenograft/';
    filter_no_survival_pat = false;
    min_mutation = 10;
else
    error('wrong dataset');
end

max_nclst = 6;
rst = 0.5;
weighted_net = false;
low_d_type=1;
input_mat_type=4;
clst_type=2;
dist_type=1;
dim = 10:10:50;



mutation_file = [survival_dir,cancer_type,'_mut.txt'];
survival_file = [survival_dir,cancer_type,'_surv.txt_clean']; %% remove patient with NA survival

[pat_info,pat_survival,adj_network,i2p,i2g,p2i,g2i] = ...
    MBS_read_data(mutation_file,network_file,survival_file,min_mutation,weighted_net,filter_no_survival_pat,false);

degree = sum(adj_network);
[npat,ngene] = size(pat_info);

if strcmp(cancer_type,'HNSC')
    HPV_status = read_TCGA_HPV('../data/tumor_stage/hpv_combo.txt',p2i);
    selected_pat = find(HPV_status==0);
else
    npat = size(pat_info,1);
    selected_pat = 1:npat;
end
[g2hgnc,i2hgnc,hgnc2i,hgnc2g] = map_gene_name(g2i,i2g,gene_name_mapping_file);
p2g = pat_info(selected_pat,:);
[npat,ngene] = size(p2g);

nneg_pat = length(selected_pat);
select_i2p = containers.Map(1:nneg_pat,values(i2p,num2cell(selected_pat)));

p2g = full(p2g);
mute_freq = sum(p2g);

topk = 100;

[ME_net,pat_diff] = MBS_read_ME_network(cancer_type,rst,pat_info,i2hgnc,topk);
[our_clust_pv,indClust,~,our_cluster_max,our_cluster_min]=...
    our_NBS_cluster_consensus(pat_info,{[pat_diff]},pat_survival,rst,...
    selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,'test');



[freq_clust_pv,indClust,~,freq_cluster_max,freq_cluster_min]=...
    our_NBS_cluster(pat_info,{[pat_info]},pat_survival,rst,...
    selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,'Mutation_profile');

pat_diff = run_diffusion(adj_network,rst,pat_info);
[diff_clust_pv,~,~,diff_cluster_max,diff_cluster_min]=...
    our_NBS_cluster(pat_info,{pat_diff},pat_survival,rst,...
    selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,'PPI_based');

[diff_clust_pv,~,~,diff_cluster_max,diff_cluster_min]=...
    our_NBS_cluster(pat_info,{pat_diff},pat_survival,rst,...
    selected_pat,i2p,degree,cancer_type,max_nclst,2,1,1,clst_type,dist_type,50,'our_NBS');


topk = 100;
[ME_net,pat_diff] = MBS_read_ME_network(cancer_type,rst,pat_info,i2hgnc,topk,0);

[our_clust_pv,indClust,~,our_cluster_max,our_cluster_min]=...
    our_NBS_cluster(pat_info,{[pat_diff]},pat_survival,rst,...
    selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,'ME_based');
topk = 100;

[ME_net,pat_diff] = MBS_read_ME_network(cancer_type,rst,pat_info,i2hgnc,topk);
[our_clust_pv,indClust,~,our_cluster_max,our_cluster_min]=...
    our_NBS_cluster(pat_info,{[pat_diff]},pat_survival,rst,...
    selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,'ME_Network_based');
