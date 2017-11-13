addpath('matlab_submine')
addpath('steiner_tree')
addpath('evaluation')
addpath('graph_search')
addpath('find_subnet')
addpath('permute_matrix')
addpath('simulated_experiment')
addpath('util')
addpath('NBS')
addpath('Hotnet')
gene_name_mapping_file = '../data/util/gene_name_mapping.txt';
use_NBS = false;
write_mut2exp = false;
test_our_method = false;
perform_subtype_analysis = true;


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


analysis_ME = false;
use_rwr_network = false;
use_existing_DCA_sim_matrix = true;
compare_with_hotnet2 = false;
%(3,50,nc=2,4)
%0,3,1,2,2,
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
    if ~exist(['../output/clinical/',method_name,'/'],'dir')
        mkdir(['../output/clinical/',method_name,'/']);
    end
fasso = fopen(['../output/clinical/',method_name,'/','association.txt'],'w');
for cv = {'BLCA','SARC','STAD','UCEC','GBM','HNSC','KIPAN','KIRC','LIHC','LUAD','LUSC','OV','PAAD'}
    cancer_type = char(cv);
    nclst_l = [2,3,4,5,6];
    for nclst = nclst_l
        mutation_file = [all_survival_dir,cancer_type,'_mut.txt'];
        survival_file = [all_survival_dir,cancer_type,'_surv.txt_clean'];
        
        [pat_info,pat_survival,adj_network,i2p,i2g,p2i,g2i] = ...
            MBS_read_data(mutation_file,network_file,survival_file,min_mutation,false,filter_no_survival_pat,false);
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
        indClust_prefix = [cancer_type,'_',num2str(nclst),'_',num2str(rst),'_',num2str(diff_mat_type),'_',num2str(input_mat_type), ...
            '_',num2str(clst_type), ...
            '_',num2str(dim(end)),'_',num2str(dist_type),method_name,'.txt'];
        subtype_folder =  ['../output/subtype/',cancer_type,'/'];
        label_file = [subtype_folder, indClust_prefix];
        [~,~,train_sur,~,~,~,p2i,i2p] = ...
            read_data(mutation_file,label_file,network_file);
        
        [g2hgnc,i2hgnc,hgnc2i,hgnc2g] = map_gene_name(g2i,i2g,gene_name_mapping_file);
 
        
        pat_label = train_sur(:,3);
        clinical_file = ['../data/cBioPortal_clinical/',lower(char(cv)),'_clinical.txt'];
        T = readtable(clinical_file,'ReadVariableNames',true,'ReadRowNames',true);
        pat_name = T.Properties.RowNames;
        filter = isKey(p2i,pat_name);
        pat_id = cell2mat(values(p2i,pat_name(filter)));
        nclinical = size(T,2);
        clinical_mat = -1*ones(npat,nclinical);
        clinical_mat(pat_id,:) = table2array(T(filter,:));
        check_clincal_variable_overlap(pat_label,clinical_mat,T.Properties.VariableNames,cancer_type,nclst,fasso,true,method_name);
    end
end
fclose(fasso);