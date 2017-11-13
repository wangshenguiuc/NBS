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
max_subnet = 100;
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
use_extend = false;
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

only_1_hop = true;
fout = fopen('../inter_result/subnetwork_analysis/subnet_score_0928.txt','w');
for cv ={'BLCA'}%,'HNSC','LUAD','LIHC','OV','GBM','SARC','THCA'}
    %cv = 'LUAD';
    %nclst_l = 5;
    % cv = 'LIHC';
    nclst_l=4:4;
    cancer_type = char(cv);
    
    %         if strcmp(cancer_type,'GBM')
    %             nclst_l = [2,3];
    %         elseif strcmp(cancer_type,'BLCA')
    %             nclst_l = [2,3,4];
    %         elseif strcmp(cancer_type,'UVM')
    %             nclst_l = [2,3,4,5,6];
    %         elseif strcmp(cancer_type,'LUAD')
    %             nclst_l = [2,3,4,5,6];
    %         elseif strcmp(cancer_type,'OV')
    %             nclst_l = [2];
    %         elseif strcmp(cancer_type,'LIHC')
    %             nclst_l = [2,3];
    %         elseif strcmp(cancer_type,'HNSC')
    %             nclst_l = [3];
    %         elseif strcmp(cancer_type,'SARC')
    %             nclst_l = [2,3,4];
    %         elseif strcmp(cancer_type,'THCA')
    %             nclst_l = [2,3,4];
    %         else
    %             nclst_l = 2;
    %         end
    % nclst_l = 2:6;
    for nclst = nclst_l
        if ~exist('../inter_result/subnetwork_analysis/','dir')
            mkdir('../inter_result/subnetwork_analysis/');
        end
        
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
        %             if  ~exist('DCA_based_net','var') || size(DCA_based_net,1)~=ngene
        %                 if ~use_existing_DCA_sim_matrix || ~exist(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio_top',num2str(DCA_node_topk),'.mat'],'file')~=0
        %                     gene_sim = dlmread(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio']);
        %                     [o,v] = sort(gene_sim,2,'descend');
        %                     u = (1:ngene);
        %                     st = repmat(u,1,DCA_node_topk);
        %                     ed = v(:,1:DCA_node_topk);
        %                     wt = o(:,1:DCA_node_topk);
        %                     DCA_based_net = sparse(st(:),ed(:),wt(:),ngene,ngene);
        %                     DCA_based_net = DCA_based_net.*(net>0);
        %                     save(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio_top',num2str(DCA_node_topk)],'DCA_based_net');
        %                 else
        %                     load(['../inter_result/DCA_sim_matrix/',cancer_type,'_Inbio_top',num2str(DCA_node_topk)]);
        %                 end
        %             end
        
        [cancer_gene,non_cancer_gene] = read_cancer_driver_gene(cancer_type,'../data/subgraph/cancer_gene.txt',hgnc2i);
        
        %     [large_ME_net,~] = MBS_read_ME_network(cancer_type,rst_prob,all_train_p2m,i2hgnc,inf,2);
        %     fprintf('n background edge=%d\n',length(find(large_ME_net(:)>0)));
        [ME_net,pat_diff] = MBS_read_ME_network(cancer_type,rst_prob,all_train_p2m,i2hgnc,100);
        
        if only_1_hop
            ME_net = ME_net.*net;
            ME_net = double(ME_net>0);
        end
        norm_pat_diff = log(pat_diff+1/(ngene*ngene))-log(1/(ngene*ngene));
        sdiff = repmat(sum(norm_pat_diff,2),1,size(norm_pat_diff,2));
        norm_pat_diff = norm_pat_diff./sdiff;
        %     pat_diff = norm_pat_diff;
        %     degree_norm_ME_net =  diag(sum(ME_net,2).^-0.5)*ME_net*diag(sum(ME_net,2).^-0.5);
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
                %              diff_prop_pv(c,i) = ranksum(vec_c,vec_non_c,'Tail','right');
            end
            diff_mute_sign_gene{c} = find(diff_prop_pv(c,:)<0.05);
        end
        %     diff_prop_pv = diff_prop_pv*ngene;
        
        %     cluster_heatmap_file = [cancer_type,'_',num2str(nclst),'.png'];
        %     pat_diff_all = full(run_diffusion(network,0.8,p2m));
        %     plot_exp_signature(diff_mute_sign_gene,pat_diff,pat_label,cluster_heatmap_file,cancer_type,2);
        %
        mutation_count = zeros(nclst,ngene);
        relevant_mutation_count = ones(nclst,ngene);
        fisher_test_pv = ones(nclst,ngene);
        c_pat = cell(1,nclst);
        c_not_pat = cell(1,nclst);
        for c=1:nclst
            c_pat{c} = find(pat_label==c);
            c_not_pat{c} = find(pat_label~=c);
            mutation_count(c,:)= sum(p2m(c_pat{c},:));
            relevant_mutation_count(c,:)= sum(p2m(c_pat{c},:))./length(c_pat{c})-sum(p2m)/npat;
            for g=1:ngene
                px = p2m(:,g);
                if sum(px)==0
                    continue
                end
                gx = pat_label==c;
                tbl = crosstab(px,gx);
                [h,pv,stats] = fishertest(tbl,'Tail','right');
                fisher_test_pv(c,g) = pv;
            end
        end
        
        %     for c=1:nclst
        %         gene_ranking_folder = ['../output/gene_ranking/',cancer_type,'/'];
        %         if exist(gene_ranking_folder,'dir')~=7
        %             mkdir(gene_ranking_folder)
        %         end
        %         node_gene_rank = diff_prop_pv(c,:);
        %         %         steiner_tree_write_gene_score(1:ngene,[node_gene_rank;diff_prop_pv(c,:);degree;ME_net_degree;mutation_count(c,:);diff_prop_edge_pv(c,:)],i2hgnc,[gene_ranking_folder,'/',cancer_type,'_',num2str(nclst),'_',num2str(c),'_',num2str(rst_prob),'_',num2str(input_mat_type),'_',num2str(dim),'_',num2str(dist_type)]);
        %     end
        
        %     pv_cutoff_l = [0.1];
        if nclst<=2
            pv_cutoff_l = 0.01;
            mutation_cutoff = 3;
        elseif nclst==3
            pv_cutoff_l = 0.05;
            mutation_cutoff = 3;
        else
            pv_cutoff_l = 0.1;
            mutation_cutoff = 1;
        end
        %         mutation_cutoff = 0.02*npat;
        %     pv_cutoff_l = 0.05;
        ME_net_gene = find(sum(ME_net)>0);
        
        for c=1:nclst
            subnet_score = [];
            aggregate_path_clst = {};
           
            pv_cutoff = pv_cutoff_l;
            diff_sign_gene = find(diff_prop_pv(c,:)<pv_cutoff);
            ME_net_gene = find(sum(ME_net)>0);
            prop_gene = intersect(ME_net_gene,diff_sign_gene);
            prop_gene = intersect(prop_gene,find(mutation_count(c,:)>mutation_cutoff));
            
            
            %                 hop1_ME_net = sparse(zeros(ngene,ngene));
            %                 for g1 = prop_gene
            %                     for g2 = prop_gene
            %                         if ME_net(g1,g2)==1
            %                             if net(g1,g2)>0
            %                                 hop1_ME_net(g1,g2)=1;
            %                                 continue
            %                             end
            %                             if sum(mutation_count(:,g1))<mutation_cutoff || sum(mutation_count(:,g2))<mutation_cutoff
            %                                 continue
            %                             end
            %                             hop1_ME_net(g1,g2) = 1;
            %                             g1_ngh = find(net(g1,:)>0);
            %                             g2_ngh = find(net(g2,:)>0);
            %                             shared_ngh = intersect(g1_ngh,g2_ngh);
            % %                             shared_ngh = intersect(shared_ngh,find(diff_prop_pv(c,:)<pv_cutoff_l));
            %                             if isempty(shared_ngh)
            %                                 continue
            %                             end
            %                             [o,v] = sort(mutation(c,shared_ngh),'ascend');
            %                             hop1_ME_net(g1,shared_ngh(v(1))) = 1;
            %                             hop1_ME_net(g2,shared_ngh(v(1))) = 1;
            %                             diff_sign_gene = [diff_sign_gene,shared_ngh(v(1))];
            %                         end
            %                     end
            %                 end
            %                 diff_sign_gene = unique(diff_sign_gene);
            %                 hop1_ME_net = max(hop1_ME_net,hop1_ME_net');
            hop1_ME_net = ME_net;
            %             ME_net_gene = find(sum(hop1_ME_net)>0);
            %             prop_gene = intersect(ME_net_gene,diff_sign_gene);
            
            
            connect_mat = double(hop1_ME_net(prop_gene,prop_gene));
            iter = 1;
            while true
                last_connect_mat = connect_mat;
                connect_mat = last_connect_mat + last_connect_mat * double(hop1_ME_net(prop_gene,prop_gene));
                connect_mat = double(connect_mat>0);
                if isequal(connect_mat,last_connect_mat)
                    break;
                end
                iter = iter + 1;
            end
            ME_net_extend = sparse(zeros(ngene,ngene));
            ME_net_extend(prop_gene,prop_gene) = connect_mat;
            nME_Edge = [];
            subnet_pv = [];
            mutation_pv = [];
            subnet_gene_set = [];
            while ~isempty(prop_gene)
                root = prop_gene(1);
                ngh = find(ME_net_extend(root,:)>0);
                subnet = [root,ngh];
                subnet = intersect(prop_gene,subnet);
                prop_gene = setdiff(prop_gene,subnet);
                nME_Edge = [nME_Edge,sum(sum(hop1_ME_net(subnet,subnet)))];
                chi_vals = -2.*log(min(1,diff_prop_pv(c,subnet)));
                group_pval = 1 - chi2cdf(sum(chi_vals),2*length(diff_prop_pv(c,subnet)));
                subnet_pv = [subnet_pv,group_pval];
                chi_vals = -2.*log(min(1,fisher_test_pv(c,subnet)));
                group_pval = 1 - chi2cdf(sum(chi_vals),2*length(fisher_test_pv(c,subnet)));
                mutation_pv = [mutation_pv,group_pval];
                aggregate_path_clst = [aggregate_path_clst,subnet];
                subnet_gene_set = [subnet_gene_set,subnet];
            end
            nME_Edge = full(nME_Edge);
            chi_vals = full(chi_vals);
%             if use_extend
 single_gene = setdiff(find(diff_prop_pv(c,:)<pv_cutoff_l),unique(subnet_gene_set));
 [o,v] = sort(mutation_count(c,single_gene),'descend');

                for gi = 1:5-length(aggregate_path_clst)
                    g = single_gene(v(gi));
                    nME_Edge = [nME_Edge,0];
                    aggregate_path_clst = [aggregate_path_clst,g];
                    subnet_pv = [subnet_pv,diff_prop_pv(c,g)];
                    mutation_pv = [mutation_pv,fisher_test_pv(c,g)];
                end
%             end
            [subnet_o,subnet_v]=sort(nME_Edge,'descend');
            nME_Edge = nME_Edge(subnet_v);
            subnet_pv = subnet_pv(subnet_v);
            aggregate_path_clst = aggregate_path_clst(subnet_v);
            subnet_folder = ['../output/new_subnet/',cancer_type,'/'];
            if exist(subnet_folder,'dir')~=7
                mkdir(subnet_folder)
            end
            output_file = [subnet_folder,num2str(nclst),'_',num2str(c),method_name,'0904_robust_pv'];
            if only_1_hop
                output_file = [output_file,'_only_1_hop'];
            end
            if use_extend
                output_file = [output_file,'.extend'];
            end
            for ni = 1:length(aggregate_path_clst)
                subnet = aggregate_path_clst{ni};
                if length(subnet)==1
                    continue
                end
                fprintf(fout,'%s,nclst:%d,c:%d,p:%d,number of ME edge:%d,propogate pvalue:%e,mutatin pv:%e\n',cancer_type,nclst,c,ni,nME_Edge(ni),subnet_pv(ni),mutation_pv(ni));
            end
            
            pathway_benchmark_file ={ '../data/subgraph/nci_pathway_hgnc.txt'};
            [ our_prec,our_recall,our_f1,our_subnet_pv] =  evaluate_pathway_enrichment([mutation_count(c,:);ones(1,ngene)*npat;mutation_count(c,:)/npat;sum(mutation_count); sum(mutation_count)/npat;diff_prop_pv(c,:)],degree,aggregate_path_clst,hgnc2i,output_file,i2hgnc,500,0.05,cancer_gene,non_cancer_gene,pathway_benchmark_file);
            
        end
    end
end
fclose(fout);