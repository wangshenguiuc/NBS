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
nclst = 6;
diff_mat_type = 3;
random_network_type_l = {'Random_network_preserve_degree','Random_network_preserve_degree_distribution'};
for ri=2:2
    random_network_type = random_network_type_l{ri};
    %Random network(preserve degree of each node)
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
    our_chisq = zeros(nperm,max_nclst);
    for iter=1:nperm
        fprintf('%s:%d\n',cancer_type,iter);
        random_net = get_shuffle_network(adj_network,random_network_type);
        
        random_pv_left_mat = max(pv_left_mat,1 - random_net);
        
        [o,v] = sort(random_pv_left_mat(:),'ascend');
        pv_thres = o(topk);
        all_pv_net = random_pv_left_mat <= min(0.99999,pv_thres);
        pat_diff = run_diffusion(all_pv_net,rst,pat_info);
        nan_pat = isnan(sum(pat_diff,2));
        pat_diff(nan_pat,:) = 0;
        
        [our_clust_pv,indClust,our_chisq(iter,:),our_cluster_max,our_cluster_min]=...
            our_NBS_cluster(pat_info,{pat_diff},pat_survival,rst,...
            selected_pat,i2p,degree,cancer_type,max_nclst,2,low_d_type,input_mat_type,clst_type,dist_type,dim,random_network_type);
    end
 dlmwrite(['../inter_result/single_gene_analysis/',cancer_type,random_network_type,'0902.txt'],our_chisq,'delimiter','\t');

    our_chisq_mean = mean(our_chisq);
    random_clust_pv = zeros(1,max_nclst);
    
    NBS_prefix = [cancer_type,'_',num2str(max_nclst),'_',num2str(rst),'_',num2str(diff_mat_type),'_',num2str(input_mat_type), ...
        '_',num2str(clst_type), ...
        '_',num2str(dim(end)),'_',num2str(dist_type),random_network_type,'.txt'];
    cluster_pv_folder = ['../output/NBS_comparision/',cancer_type,'/'];
    if exist(cluster_pv_folder,'dir')~=7
        mkdir(cluster_pv_folder)
    end
    
   
    fout = fopen([cluster_pv_folder,NBS_prefix],'w');
    for nc = 2:max_nclst
        random_clust_pv = 1 - chi2cdf(our_chisq_mean(nc),nc-1);
%         fprintf(flog,'nc=%d\tpv=%f\n',nc,random_clust_pv);
%         for p=1:nperm
%             fprintf(flog,'pv=%f\n',our_chisq(p,nc-1));
%         end
        fprintf(fout,'%s\t%d\t%d\t%f\t',cancer_type,dim(end),nc,rst);
        output_vec = zeros(1,16);
        output_vec(12) = random_clust_pv;
        output_vec(16) = our_chisq_mean(nc);
        for r = 1:16
            fprintf(fout,'%e\t',output_vec(r));
        end
        fprintf(fout,'\n');
        
    end
    fclose(fout);
    
end
