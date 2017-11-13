addpath('matlab_submine')
addpath('evaluation')
addpath('graph_search')
addpath('find_subnet')
addpath('permute_matrix')
addpath('util')
addpath('NBS')
addpath('Hotnet')
gene_name_mapping_file = '../data/util/gene_name_mapping.txt';
use_NBS = false;
write_mut2exp = false;
test_our_method = false;
perform_subtype_analysis = true;

depth = 5;

analysis_ME = false;
use_rwr_network = false;
compare_with_hotnet2 = false;
%(3,50,nc=2,4)
%0,3,1,2,2,
method_name = 'ME_Network_based';
% method_name = 'Mutation_profile';
rst = 0.5;
network_file = '../data/network/InBio-Map_Entrez.sif';
network_type = 'eqweight_dnnet_InBio_Map';
clst_type = 2;
dist_type = 1;
low_d_type = 1;
dim = 50;
input_mat_type = 4;
diff_mat_type = 3;
nclst = 6;
max_nroot_gene = 10;

dataset = 'TCGA';
if strcmp(dataset,'TCGA')
    all_survival_dir = '../data/survival/raw_survival/';
    filter_no_survival_pat = true;
    min_mutation = 10;
elseif  strcmp(dataset,'ICGC')
    all_survival_dir = '../data/survival/ICGC/';
    filter_no_survival_pat = false;
    min_mutation = 10;
else
    error('wrong dataset');
end
use_filter_diff_mat = false;
network_type = 'eqweight_dnnet_InBio_Map';
use_existing_DCA_sim_matrix = false;
max_nclst = 6;
for cv = {'KIPAN','KIRC','LIHC','LUAD','LUSC','OV','PAAD','SARC','UCEC','BLCA','HNSC'}
    cancer_type = char(cv);
    nclst_l = 2:max_nclst;
    pat_label = [];
    for nclst = nclst_l
        mutation_file = [all_survival_dir,cancer_type,'_mut.txt'];
        survival_file = [all_survival_dir,cancer_type,'_surv.txt_clean'];
        [all_train_p2m,~,net,all_i2p,i2g,all_p2i,g2i] = read_data_no_survival(mutation_file,network_file,survival_file,min_mutation,false,false,false,filter_no_survival_pat);
        network = diag(sum(net,2).^-0.5)*net*diag(sum(net,2).^-0.5);
        pat_diff = run_diffusion(network,rst_prob,all_train_p2m);
        if strcmp(cancer_type,'HNSC')
            HPV_status = read_TCGA_HPV('../data/tumor_stage/hpv_combo.txt',all_p2i);
            selected_pat = find(HPV_status==0);
            %%we only selected HPV negative patients
        else
            npat = size(all_train_p2m,1);
            selected_pat = 1:npat;
        end
        p2m = all_train_p2m(selected_pat,:);
      indClust_prefix = [cancer_type,'_',num2str(nclst),'_',num2str(rst),'_',num2str(diff_mat_type),'_',num2str(input_mat_type), ...
            '_',num2str(clst_type), ...
            '_',num2str(dim(end)),'_',num2str(dist_type),method_name,'.txt'];
        subtype_folder =  ['../output/subtype/',cancer_type,'/'];
        label_file = [subtype_folder, indClust_prefix];
        [~,~,train_sur,~,~,~,p2i,i2p] = ...
            read_data(mutation_file,label_file,network_file);
        
        [g2hgnc,i2hgnc,hgnc2i,hgnc2g] = map_gene_name(g2i,i2g,gene_name_mapping_file);
        [npat,ngene] = size(p2m);
        degree = full(sum(net));
        [cancer_gene,non_cancer_gene] = read_cancer_driver_gene(cancer_type,'../data/subgraph/cancer_gene.txt',hgnc2i);
        pat_label = [pat_label,train_sur(:,3)];
    end
    subtype_size = zeros(max_nclst,max_nclst);
    subtype_map = zeros(max_nclst,max_nclst);
    subtype_map_pv = ones(max_nclst,max_nclst);
    for i=2:max_nclst
        pat_label_up = pat_label(:,i-1);
        for j=1:i
            subtype_size(i,j) = length(find(pat_label_up==j));
        end
    end
    if ~exist(['../output/clinical/',method_name,'/'],'dir')
        mkdir(['../output/clinical/',method_name,'/']);
    end
    survival_file = ['../output/clinical/',method_name,'/',cancer_type,'_',num2str(var_remove),'_',num2str(rst_prob),'_',num2str(diff_mat_type),'_',num2str(norm_type),'_',num2str(input_mat_type),'_',num2str(clst_type),'_',num2str(dim),'_',num2str(dist_type),'.map_size'];
    dlmwrite(survival_file,subtype_size,'delimiter','\t');
    
    for i=3:max_nclst
        subtype_map_pv = zeros(max_nclst,max_nclst);
        pat_label_up = pat_label(:,i-2);
        pat_label_down = pat_label(:,i-1);
        for j=1:i-1
            sub_pat_j = find(pat_label_up == j);
            for k=1:i
                sub_pat_k = find(pat_label_down == k);
                
%                 tbl = crosstab(sub_pat_j,sub_pat_k);
%                 [h,pv,stats] = fishertest(tbl,'Tail','right');
                subtype_map_pv(j,k) = jaccard_sim( sub_pat_k,sub_pat_j);
            end
            
        end
        survival_file = ['../output/clinical/',method_name,'/',cancer_type,'_',num2str(i),'.map_DAG_index'];
        dlmwrite(survival_file,subtype_map_pv,'delimiter','\t');
        %     survival_file = ['../output/clinical/',cancer_type,'_',num2str(var_remove),'_',num2str(rst_prob),'_',num2str(diff_mat_type),'_',num2str(norm_type),'_',num2str(input_mat_type),'_',num2str(clst_type),'_',num2str(dim),'_',num2str(dist_type),'.map_index'];
        %     dlmwrite(survival_file,subtype_map,'delimiter','\t');
        %     survival_file = ['../output/clinical/',cancer_type,'_',num2str(var_remove),'_',num2str(rst_prob),'_',num2str(diff_mat_type),'_',num2str(norm_type),'_',num2str(input_mat_type),'_',num2str(clst_type),'_',num2str(dim),'_',num2str(dist_type),'.map_pv'];
        %     dlmwrite(survival_file,subtype_map_pv,'delimiter','\t');
    end
end
