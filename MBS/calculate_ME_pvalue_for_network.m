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
addpath('steiner_tree_util')
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
    if pv(ii)<0.5
        pv_left_mat(gid1(ii),gid2(ii)) = pv(ii);
    end
end

our_chisq = zeros(nperm,max_nclst);
% network_list_l = {'../data/network/FI_hgnc.txt'};

network_list_l = {'../data/network/PathwayCommons9.txt','../data/network/string_integrated.txt','../data/network/FI_hgnc.txt','../data/network/string_inBio_pathComm_FI.txt'};
for network_list_c = network_list_l
    network_list = char(network_list_c);
    net_name = strsplit(network_list,'/');
    net_name = char(net_name(end));
    [e1,e2,~] = textread(network_list,'%s%s%f','delimiter','\t');
    filter = isKey(hgnc2i,e1) & isKey(hgnc2i,e2);
    e1 = e1(filter);
    e2 = e2(filter);
    eid1 = cell2mat(values(hgnc2i,e1));
    eid2 = cell2mat(values(hgnc2i,e2));
    new_net = sparse(eid1,eid2,1,ngene,ngene);
    new_net = max(new_net,new_net');
    for nhop=1:3
        adj_network = new_net;
        for ct=2:nhop
            adj_network = adj_network*adj_network;
            adj_network = min(adj_network,1);
            adj_network = max(adj_network,adj_network');
        end
        adj_network = min(adj_network,1);
        
        
        s = sum(pat_info);
        high_freq_gene = find(s>=3);
        low_freq_gene = setdiff(1:ngene,high_freq_gene);
        adj_network(:,low_freq_gene)=0;
        adj_network(low_freq_gene,:)=0;        
        
        fprintf('cancer_type=%s ngene=%d\n',cancer_type,length(high_freq_gene));
        
        [pv_left_mat,odds_mat] = fast_local_ME_test(high_freq_gene,p2g,adj_network,1,cancer_type);
        pv_left_mat = full(pv_left_mat);
        
        [vi,vj] = find(pv_left_mat<0.5);
        nedge = length(vi);
        network_folder = ['../inter_result/ME_network/',net_name,'/'];
        if ~exist(network_folder)
            mkdir(network_folder);
        end
        fout_net = fopen([network_folder,cancer_type,'_ME_',num2str(nhop),'.txt'],'w');
        for i=1:nedge
            fprintf(fout_net,'%d\t%s\t%d\t%s\t%f\t%f\n',vi(i),i2hgnc(vi(i)),vj(i),i2hgnc(vj(i)),pv_left_mat(vi(i),vj(i)),odds_mat(vi(i),vj(i)));
        end
        fclose(fout_net);        
    end
    
end


