addpath('matlab_submine')
addpath('steiner_tree_util')
addpath('steiner_tree_main')
addpath('evaluation')
addpath('graph_search')
addpath('analysis_subnet')
addpath('find_subnet')
addpath('permute_matrix')
addpath('util')
addpath('NBS')
addpath('Hotnet')
max_subnet_size = 20;
max_subnet = 10;
DCA_node_topk = 10;
rank_mute_freq_pv_thres = 0.05;
mute_freq_pv = 0.05;
node_weight = 2.5;
warning off
max_nroot_gene = 100;
max_nprop_gene = 1000;
adjust_cluster_method = 'frequency';

gene_name_mapping_file = '../data/util/gene_name_mapping.txt';
use_NBS = false;
write_mut2exp = false;
test_our_method = false;
perform_subtype_analysis = true;

method_name = 'ME_Network_based';
% method_name = 'Mutation_profile';
rst_prob = 0.5;
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

use_existing_DCA_sim_matrix = true;


% for cv = {'BLCA','KIRC','OV','READ','LUAD','HNSC','UCEC','GBM','SARC'}
%     for cv = {'BLCA'};
cancer_type = char(cv);

if strcmp(cancer_type,'GBM')
    nclst_l = [2,3];
elseif strcmp(cancer_type,'BLCA')
    nclst_l = [3,4];
elseif strcmp(cancer_type,'UVM')
    nclst_l = [2,3,4,5,6];
elseif strcmp(cancer_type,'LUAD')
    nclst_l = [4,5];
elseif strcmp(cancer_type,'OV')
    nclst_l = [2,3];
elseif strcmp(cancer_type,'LIHC')
    nclst_l = [2,3];
elseif strcmp(cancer_type,'HNSC')
    nclst_l = [3,4];
elseif strcmp(cancer_type,'SARC')
    nclst_l = [2];
elseif strcmp(cancer_type,'THCA')
    nclst_l = [2];
else
    nclst_l = 2;
end
for nclst = nclst_l
    if ~exist('../inter_result/subnetwork_analysis/','dir')
        mkdir('../inter_result/subnetwork_analysis/');
    end
    fout = fopen(['../inter_result/subnetwork_analysis/',cancer_type,num2str(nclst),method_name,'.txt'],'w');
    
    mutation_file = [all_survival_dir,cancer_type,'_mut.txt'];
    survival_file = [all_survival_dir,cancer_type,'_surv.txt_clean'];
    [all_train_p2m,pat_survival,net,all_i2p,i2g,all_p2i,g2i] = ...
        MBS_read_data(mutation_file,network_file,survival_file,min_mutation,false,filter_no_survival_pat,false);
    
    network = diag(sum(net,2).^-0.5)*net*diag(sum(net,2).^-0.5);
    
    if strcmp(cancer_type,'HNSC')
        HPV_status = read_TCGA_HPV('../data/tumor_stage/hpv_combo.txt',all_p2i);
        selected_pat = find(HPV_status==0);
        %%we only selected HPV negative patients
    else
        npat = size(all_train_p2m,1);
        selected_pat = 1:npat;
    end
    
    subtype_folder = ['../output/subtype/',cancer_type,'/'];
    NBS_prefix = [cancer_type,'_',num2str(nclst),'_',num2str(rst_prob),'_',num2str(diff_mat_type),...
        '_',num2str(input_mat_type), ...
        '_',num2str(clst_type), ...
        '_',num2str(dim(end)),'_',num2str(dist_type),method_name,'.txt'];
    label_file = [subtype_folder, NBS_prefix];
    [~,~,train_sur,~,~,~,p2i,i2p] = ...
        read_data(mutation_file,label_file,network_file);
    
    [g2hgnc,i2hgnc,hgnc2i,hgnc2g] = map_gene_name(g2i,i2g,gene_name_mapping_file);
    
    degree = full(sum(net));
    pat_label = train_sur(:,3);
    ngene = size(network,1);
    
    [cancer_gene,non_cancer_gene] = read_cancer_driver_gene(cancer_type,'../data/subgraph/cancer_gene.txt',hgnc2i);
    
    %     [large_ME_net,~] = MBS_read_ME_network(cancer_type,rst_prob,all_train_p2m,i2hgnc,inf,2);
    %     fprintf('n background edge=%d\n',length(find(large_ME_net(:)>0)));
    [ME_net,pat_diff] = MBS_read_ME_network(cancer_type,rst_prob,all_train_p2m,i2hgnc,100);
    %     if ~strcmp( method_name,'ME_Network_based')
    %         ME_net = net;
    %         large_ME_net = net;
    %     end
    %     net = large_ME_net;
    %         [our_clust_pv,indClust]=...
    %             our_NBS_cluster(all_train_p2m,{[pat_diff,all_train_p2m(:,rare_gene_set)]},pat_survival,0,rst_prob,...
    %             selected_pat,all_i2p,degree,cancer_type,6,2,low_d_type,input_mat_type,clst_type,dist_type,dim);
    ME_net_degree = sum(ME_net);
    %     isequal(pat_label,indClust(:,nclst));
    pat_label = train_sur(:,3);
    p2m = all_train_p2m(selected_pat,:);
    pat_diff = pat_diff(selected_pat,:);
    [npat,ngene] = size(p2m);
    pat_diff = full(pat_diff);
    diff_prop_pv = ones(nclst,ngene);
    %
    %     diff_prop_edge_pv = ones(nclst,ngene);
    %     for c=1:nclst
    %         cls_label = double(pat_label==c);
    %         for g = unique(rare_gene_set)
    %             gene_label = p2m(:,g);
    %             tbl = crosstab(cls_label,gene_label);
    %             [~,diff_prop_edge_pv(c,g)] = fishertest(tbl,'tail','right');
    %         end
    %     end
    
    diff_mute_sign_gene = cell(1,nclst);
    for c=1:nclst
        for i=1:ngene
            cls1 = find(pat_label==c);
            cls2 = setdiff(1:npat,cls1);
            vec_c = pat_diff(cls1,i);
            vec_non_c = pat_diff(cls2,i);
            [~,diff_prop_pv(c,i)] = ttest2(vec_c,vec_non_c,'Tail','right');
        end
        diff_mute_sign_gene{c} = find(diff_prop_pv(c,:)<0.05);
    end
    
    
    %     cluster_heatmap_file = [cancer_type,'_',num2str(nclst),'.png'];
    %     pat_diff_all = full(run_diffusion(network,0.8,p2m));
    %     plot_exp_signature(diff_mute_sign_gene,pat_diff,pat_label,cluster_heatmap_file,cancer_type,2);
    %
    mutation_count = zeros(nclst,ngene);
    relevant_mutation_count = zeros(nclst,ngene);
    c_pat = cell(1,nclst);
    for c=1:nclst
        c_pat{c} = find(pat_label==c);
        mutation_count(c,:)= sum(p2m(c_pat{c},:));
        relevant_mutation_count(c,:)= sum(p2m(c_pat{c},:))./length(c_pat{c})-sum(p2m)/npat;
    end
    
    for c=1:nclst
        gene_ranking_folder = ['../output/gene_ranking/',cancer_type,'/'];
        if exist(gene_ranking_folder,'dir')~=7
            mkdir(gene_ranking_folder)
        end
        node_gene_rank = diff_prop_pv(c,:);
        %             steiner_tree_write_gene_score(1:ngene,[node_gene_rank;diff_prop_pv(c,:);degree;ME_net_degree;mutation_count(c,:);diff_prop_edge_pv(c,:)],i2hgnc,[gene_ranking_folder,'/',cancer_type,'_',num2str(nclst),'_',num2str(c),'_',num2str(rst_prob),'_',num2str(input_mat_type),'_',num2str(dim),'_',num2str(dist_type)]);
    end
    ME_net_gene = find(sum(ME_net)>0);
    max_nroot_gene = 20;
    for c=1:nclst
        diff_sign_gene = find(diff_prop_pv(c,:)<0.05);
        prop_gene = intersect(ME_net_gene,diff_sign_gene);
        [o,v] = sort(relevant_mutation_count(c,prop_gene),'descend');
        
        root_gene_l = prop_gene(v(1:min(max_nroot_gene,length(v))));
        ct = 1;
        expand_gene_set = [];
        aggregate_path_clst = {};
        for ct=1:min(max_nroot_gene,length(root_gene_l))
            
            file_name = [cancer_type,'test'];
            subnet = root_gene_l(ct);
            if isempty(find(prop_gene==subnet, 1))
                continue
            end
            root = root_gene_l(ct);
            iter = 1;
            while iter<100
                nsubnet = length(subnet);
                ngh = intersect(prop_gene,find(sum(ME_net(subnet,:),1)>0));
                subnet = [subnet,ngh];
                subnet = unique(subnet);
                if length(subnet) == nsubnet
                    break
                end
                iter = iter + 1;
            end
            prop_gene = setdiff(prop_gene,subnet);
            
            sign_gene_pv = 1;
            sign_gene_pv_l = ones(1,ngene);
            for g=subnet
                gene_pv = plot_survival_CNA(train_sur,[cancer_type,'test'],p2m(:,g)+1,0);
                gene_pv = gene_pv(4);
                if gene_pv<sign_gene_pv && sum(p2m(:,g))>10
                    sign_gene_pv = gene_pv;
                    sign_gene_id = g;
                end
                sign_gene_pv_l(g) = gene_pv;
            end
            
            surv_pv = calculate_subnet_pvalue( subnet,ME_net,p2m,rst_prob,train_sur,file_name,pat_diff,cancer_type,i2hgnc,subnet);
            
            fprintf('c:%d,npath:%d,it:%d,%s,ngene:%d,%f\t%f\t%f\t%f\n',c,ct,iter,i2hgnc(root),length(subnet),surv_pv);
            fprintf(fout,'c:%d,npath:%d,it:%d,%s,ngene:%d,%f\t%f\t%f\t%f\n',c,ct,i2hgnc(root),iter,length(subnet),surv_pv);
            if min(surv_pv(1),surv_pv(2))<min(0.05,sign_gene_pv)
                
                for g=subnet
                    if isKey(i2hgnc,g)
                        fprintf('%s(%f %d)\t',i2hgnc(g),sign_gene_pv_l(g),sum(p2m(:,g)));
                        fprintf(fout,'%s(%f %d)\t',i2hgnc(g),sign_gene_pv_l(g),sum(p2m(:,g)));
                    end
                end
                fprintf('\n');
                fprintf(fout,'\n');
                
            end
            aggregate_path_clst = [aggregate_path_clst,subnet];
        end
        unq_aggregate_path_clst = {};
        for ii=1:length(aggregate_path_clst)
            overlap = false;
            for jj=ii+1:length(aggregate_path_clst)
                if isequal(aggregate_path_clst{ii},aggregate_path_clst{jj})
                    overlap = true;
                    break
                end
            end
            if ~overlap
                unq_aggregate_path_clst = [unq_aggregate_path_clst,aggregate_path_clst{ii}];
            end
        end
        aggregate_path_clst = unq_aggregate_path_clst;
        subnet_folder = ['../output/new_subnet/',cancer_type,'/'];
        if exist(subnet_folder,'dir')~=7
            mkdir(subnet_folder)
        end
        output_file = [subnet_folder,num2str(nclst),'_',num2str(c),'_expanded'];
        pathway_benchmark_file ={ '../data/subgraph/nci_pathway_hgnc.txt'};
        [ our_prec,our_recall,our_f1,our_subnet_pv] =  evaluate_pathway_enrichment([mutation_count(c,:);diff_prop_pv(c,:)],degree,aggregate_path_clst,hgnc2i,output_file,i2hgnc,500,0.05,cancer_gene,non_cancer_gene,pathway_benchmark_file);
    end
    fclose(fout);
end
% end
