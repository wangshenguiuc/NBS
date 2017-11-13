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
dataset = 'TCGA';
if strcmp(dataset,'TCGA')
    survival_dir = '../data/survival/raw_survival/';
    filter_no_survival_pat = true;
    min_mutation = 10;
elseif  strcmp(dataset,'ICGC')
    survival_dir = '../data/survival/ICGC/';
    filter_no_survival_pat = false;
    min_mutation = 10;
elseif  strcmp(dataset,'msk')
    survival_dir = '../data/survival/msk/';
    filter_no_survival_pat = false;
    min_mutation = 1;
elseif  strcmp(dataset,'Xenograft')
    survival_dir = '../data/survival/Xenograft/';
    filter_no_survival_pat = false;
    min_mutation = 10;
else
    error('wrong dataset');
end

network_file = '../data/network/InBio-Map_Entrez.sif';
% network_file = '..\data\network\FI_id.txt';
gene_name_mapping_file = '../data/util/gene_name_mapping.txt';
result_dir = '../output/subtype/';
dataset = 'TCGA';
max_nclst = 6;
rst = 0.5;
weighted_net = false;
low_d_type=1;
input_mat_type=4;
clst_type=2;
dist_type=1;
dim = 10:10:50;
% query_cancer_type = 'batch';
nhop = 0;
files = dir([survival_dir '*.txt_clean'])';
for fi = 1:length(files)
    f = files(fi);
    if length(f.name)<5
        continue
    end
    fl = strsplit(f.name,'.');
    fl = fl{1};
    fl = strsplit(fl,'_');
    cancer_type = fl{1};

    if ~strcmp(query_cancer_type,'batch') && ~strcmp(cancer_type,query_cancer_type)
        continue
    end
    mutation_file = [survival_dir,cancer_type,'_mut.txt'];
    survival_file = [survival_dir,cancer_type,'_surv.txt_clean']; %% remove patient with NA survival
    
    [pat_info,pat_survival,adj_network,i2p,i2g,p2i,g2i] = ...
        MBS_read_data(mutation_file,network_file,survival_file,min_mutation,weighted_net,filter_no_survival_pat,false);
    [npat,ngene] = size(pat_info);
    if npat<100
        fprintf('%s not enough patient\n',cancer_type);
        continue
    end
    fprintf('cancer_type=%s=%s start compute adjnetwork\n',cancer_type,query_cancer_type);
    if nhop == 0;
        adj_network = ones(ngene,ngene);
    else
        for ct=2:nhop
            adj_network = adj_network*adj_network;
            adj_network = min(adj_network,1);
            adj_network = max(adj_network,adj_network');
        end
    end
    degree = sum(adj_network);
    
    fprintf('cancer_type=%s finished compute adjnetwork\n',cancer_type);
    
    if strcmp(cancer_type,'HNSC')
        % for HPV, we only consider HPV negative patient
        HPV_status = read_TCGA_HPV('../data/tumor_stage/hpv_combo.txt',p2i);
        selected_pat = find(HPV_status==0);
    else
        npat = size(pat_info,1);
        selected_pat = 1:npat;
    end
    [g2hgnc,i2hgnc,hgnc2i,hgnc2g] = map_gene_name(g2i,i2g,gene_name_mapping_file);
    p2g = pat_info(selected_pat,:);
    [npat,ngene] = size(p2g);
    p2g = full(p2g);
    
    s = sum(pat_info);
    high_freq_gene = find(s>5);
    	low_freq_gene = setdiff(1:ngene,high_freq_gene);
	adj_network(:,low_freq_gene)=0;
	adj_network(low_freq_gene,:)=0;
%     if length(high_freq_gene)>1100
%         high_freq_gene = find(s>10);
%     elseif length(high_freq_gene)<100
%         high_freq_gene = find(s>1);
%     end
    
    
    fprintf('cancer_type=%s ngene=%d\n',cancer_type,length(high_freq_gene));
    
    [pv_left_mat,odds_mat] = fast_local_ME_test(high_freq_gene,p2g,adj_network,1,cancer_type);
    pv_left_mat = full(pv_left_mat);
    
    [vi,vj] = find(pv_left_mat<0.5);
    nedge = length(vi);
    fout_net = fopen(['../inter_result/ME_network/',cancer_type,'_ME_',num2str(nhop),'.txt'],'w');
    for i=1:nedge
        fprintf(fout_net,'%d\t%s\t%d\t%s\t%f\t%f\n',vi(i),i2hgnc(vi(i)),vj(i),i2hgnc(vj(i)),pv_left_mat(vi(i),vj(i)),odds_mat(vi(i),vj(i)));
    end
    fclose(fout_net);
end