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
max_subnet = 20;
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
max_nroot_gene = 30;

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
    nclst_l = [2,3,4];
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
    %         if  ~exist('DCA_based_net','var') || size(DCA_based_net,1)~=ngene
    %             if ~use_existing_DCA_sim_matrix || ~exist(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio_top',num2str(DCA_node_topk),'.mat'],'file')~=0
    %                 gene_sim = dlmread(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio']);
    %                 [o,v] = sort(gene_sim,2,'descend');
    %                 u = (1:ngene);
    %                 st = repmat(u,1,DCA_node_topk);
    %                 ed = v(:,1:DCA_node_topk);
    %                 wt = o(:,1:DCA_node_topk);
    %                 DCA_based_net = sparse(st(:),ed(:),wt(:),ngene,ngene);
    %                 DCA_based_net = DCA_based_net.*(net>0);
    %                 save(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio_top',num2str(DCA_node_topk)],'DCA_based_net');
    %             else
    %                 load(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio_top',num2str(DCA_node_topk)]);
    %             end
    %         end
    
    [cancer_gene,non_cancer_gene] = read_cancer_driver_gene(cancer_type,'../data/subgraph/cancer_gene.txt',hgnc2i);
    
    [large_ME_net,~] = MBS_read_ME_network(cancer_type,rst_prob,all_train_p2m,i2hgnc,inf,2);
    fprintf('n background edge=%d\n',length(find(large_ME_net(:)>0)));
    [ME_net,pat_diff] = MBS_read_ME_network(cancer_type,rst_prob,all_train_p2m,i2hgnc,100);
    if ~strcmp( method_name,'ME_Network_based')
        ME_net = net;
        large_ME_net = net;
    end
    net = large_ME_net;
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
        %         steiner_tree_write_gene_score(1:ngene,[node_gene_rank;diff_prop_pv(c,:);degree;ME_net_degree;mutation_count(c,:);diff_prop_edge_pv(c,:)],i2hgnc,[gene_ranking_folder,'/',cancer_type,'_',num2str(nclst),'_',num2str(c),'_',num2str(rst_prob),'_',num2str(input_mat_type),'_',num2str(dim),'_',num2str(dist_type)]);
    end
    ME_net_gene = find(sum(ME_net)>0);
    for c=1:nclst
        diff_sign_gene = find(diff_prop_pv(c,:)<0.1);
        prop_gene = diff_sign_gene;
        
        root_gene = intersect(ME_net_gene,diff_sign_gene);
		[o,v] = sort(mutation_count(c,root_gene),'descend');
        root_gene_l = root_gene(v(1:min(max_nroot_gene,length(v))));
        ct = 1;
        expand_gene_set = [];
        aggregate_path_clst = {};
        for ct=1:min(max_nroot_gene,length(root_gene_l))
            root_gene = root_gene_l(ct);
            file_name = ['real_',cancer_type,'_',num2str(nclst),'_',num2str(c),'_',num2str(ct)];
            
            [our_path_clst,tree_net] = write_msg_steiner_input(hgnc2i,ME_net,ME_net,net,diff_prop_pv(c,:),degree,root_gene, ...
                prop_gene,i2hgnc,file_name,0.5,cancer_type,4);
            
            if isempty(our_path_clst)
                continue
            end
            
            subnet = adjust_cluster_size(our_path_clst{1},net,net,root_gene, ...
                mutation_count(c,:),20,adjust_cluster_method);
            subnet = reshape(subnet,1,length(subnet));
           
            subnet_kmplot_folder = ['../output/subnet_kmplot/',cancer_type,'_',num2str(nclst),'_',num2str(c),'/'];
            if exist(subnet_kmplot_folder,'dir')~=7
                mkdir(subnet_kmplot_folder)
            end
%             subnet_cluster_method = 1;
%              if subnet_cluster_method==1
%                  sub_pat_diff = run_diffusion(ME_net(subnet,subnet),rst_prob,all_train_p2m(:,subnet));
%                  nan_pat = isnan(sum(sub_pat_diff,2));
%                  sub_pat_diff(nan_pat,:) = 0;
%                  subnet_pat_label = kmeans(squareform(pdist(sub_pat_diff)),2,'replicates',100,'distance','cosine');
%              end      
%               ICGC_write_to_file = false;
%            check_ICGC
%             subnet_pat_label = kmeans(squareform(pdist([pat_diff(:,subnet)])),2,'replicates',100,'distance','cosine');
 %           file_name = [subnet_kmplot_folder,num2str(nclst),'_',num2str(c),'_',num2str(max_subnet_size),'_',char(i2hgnc(root_gene))];
%            subnet_surv_pv = plot_survival_CNA(train_sur,file_name,subnet_pat_label,0);
%            subnet_surv_pv = subnet_surv_pv(4);
            
%            subnet_heat_folder = ['../output/subnet_heatmap/',cancer_type,'/'];
%            if exist(subnet_heat_folder,'dir')~=7
 %               mkdir(subnet_heat_folder)
 %           end
           expand_gene_set = [expand_gene_set,subnet];
%             prop_gene = setdiff(prop_gene,expand_gene_set);
%             fprintf(fout,'ct:%d,root gene:%s subnet_surv_pv:%f,ngene=%d\n',ct,char(values(i2hgnc,{root_gene})),...
%                 subnet_surv_pv,length(subnet));
%             fprintf('ct:%d,root gene:%s subnet_surv_pv:%f,ngene=%d\n',ct,char(values(i2hgnc,{root_gene})),...
%                 subnet_surv_pv,length(subnet));
%             if subnet_surv_pv < 0.05
%                 for gg = subnet
%                     gene_pv =  plot_survival_CNA(train_sur,[cancer_type,'test'],all_train_p2m(:,gg),0);
%                     gene_pv = gene_pv(4);
%                     fprintf(fout,'%s\t',i2hgnc(gg));
%                     fprintf('%s\t%f\t%d\n',i2hgnc(gg),gene_pv,full(sum(all_train_p2m(:,gg))));
%                 end
%                 fprintf(fout,'\n');
%                 fprintf('\n');
%             end
%             file_name = [subnet_heat_folder,num2str(nclst),'_',num2str(c),'_',num2str(max_subnet_size),'_',char(i2hgnc(root_gene))];
%             
            aggregate_path_clst = [aggregate_path_clst,subnet];
            if length(aggregate_path_clst)>max_subnet
                break;
            end
        end
        
        subnet_folder = ['../output/new_subnet/',cancer_type,'/'];
        if exist(subnet_folder,'dir')~=7
            mkdir(subnet_folder)
        end
        output_file = [subnet_folder,num2str(nclst),'_',num2str(c),method_name,'0831'];
        pathway_benchmark_file ={ '../data/subgraph/nci_pathway_hgnc.txt'};
        [ our_prec,our_recall,our_f1,our_subnet_pv] =  evaluate_pathway_enrichment([mutation_count(c,:);diff_prop_pv(c,:)],degree,aggregate_path_clst,hgnc2i,output_file,i2hgnc,500,0.05,cancer_gene,non_cancer_gene,pathway_benchmark_file);
        
    end
    fclose(fout);
end
% end
