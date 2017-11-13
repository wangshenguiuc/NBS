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
survival_dir = '../data/survival/raw_survival/';
filter_no_survival_pat = true;
min_mutation = 10;
max_nclst = 6;
rst = 0.5;
weighted_net = false;
low_d_type=1;
input_mat_type=4;
clst_type=2;
dist_type=1;
dim = 10:10:50;
nperm = 10;
diff_mat_type = 3;

mutation_file = [survival_dir,cancer_type,'_mut.txt'];
survival_file = [survival_dir,cancer_type,'_surv.txt_clean']; %% remove patient with NA survival

[pat_info,pat_survival,adj_network,i2p,i2g,p2i,g2i] = ...
    MBS_read_data(mutation_file,network_file,survival_file,min_mutation,weighted_net,filter_no_survival_pat,false);

degree = sum(adj_network);

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
ME_file =  ['../inter_result/ME_network/',cancer_type,'_ME_',num2str(0),'.txt'];
[gid1,gname1,gid2,gname2,pv,odds] = textread(ME_file,'%d%s%d%s%f%f');
nnet_edge =length(gid1);
pv_left_mat = ones(ngene,ngene);
for ii=1:nnet_edge
    pv_left_mat(gid1(ii),gid2(ii)) = pv(ii);
end

topk = 100;

[e1,e2,~] = textread('../data/network/InBio-Map_Symbol.sif','%s%s%f','delimiter','\t');
fout = fopen('../data/network/InBio-Map_Symbol_edge.sif','w');
for i=1:length(e1)
    fprintf(fout,'%s\t%s\t%d\n',e1{i},e2{i},1);
end
fclose(fout);
our_chisq = zeros(nperm,max_nclst);
% network_list_l = {'../data/network/FI_hgnc.txt'};
network_list_l = {'../data/network/string_integrated.txt'};%,'../data/network/FI_hgnc.txt','../data/network/string_integrated.txt'};
for network_list_c = network_list_l
    network_list = char(network_list_c);
    net_name = strsplit(network_list,'/');
    net_name = char(net_name(end));
%     [e1,e2,~] = textread(network_list,'%s%s%f','delimiter','\t');
%     filter = isKey(hgnc2i,e1) & isKey(hgnc2i,e2);
%     e1 = e1(filter);
%     e2 = e2(filter);
%     eid1 = cell2mat(values(hgnc2i,e1));
%     eid2 = cell2mat(values(hgnc2i,e2));
%     new_net = sparse(eid1,eid2,1,ngene,ngene);
%     new_net = max(new_net,new_net');
%     new_pv_left_mat = max(pv_left_mat,1 - new_net);
%     
%     [o,v] = sort(new_pv_left_mat(:),'ascend');
%     pv_thres = o(topk);
%     all_pv_net = new_pv_left_mat <= min(0.99999,pv_thres);
%     pat_diff = run_diffusion(all_pv_net,rst,pat_info);
%     nan_pat = isnan(sum(pat_diff,2));
%     pat_diff(nan_pat,:) = 0;
    
    [ME_net,pat_diff] = MBS_read_ME_network([net_name,'/',cancer_type],rst,pat_info,i2hgnc,topk);
    
    [our_clust_pv,indClust,our_chisq,our_cluster_max,our_cluster_min]=...
        our_NBS_cluster(pat_info,{pat_diff},pat_survival,rst,...
        selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,...
        [net_name,'_ME']);
    
%     all_pv_net = new_net;
%     pat_diff = run_diffusion(all_pv_net,rst,pat_info);
%     nan_pat = isnan(sum(pat_diff,2));
%     pat_diff(nan_pat,:) = 0;    
%     [our_clust_pv,indClust,our_chisq,our_cluster_max,our_cluster_min]=...
%         our_NBS_cluster(pat_info,{pat_diff},pat_survival,rst,...
%         selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,...
%         net_name);
%     flog = fopen(['../inter_result/single_gene_analysis/',cancer_type,net_name],'w');
%     for nc = 2:max_nclst
%         fprintf(flog,'nc=%d\tpv=%f\n',nc,our_clust_pv(nc));
%     end
%     fclose(flog);
end


