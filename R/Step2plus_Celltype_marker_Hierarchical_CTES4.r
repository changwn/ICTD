#tg_R1_lists
#data_CORS_cancer_s
#data.matrix=data_logp1n2_s
#tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]]
#cell_type_enrich_cut=0.5
#tg_R1_cut=6
#ext_cn=5
#CT_balance=1



Step2plus_Celltype_marker_Hierarchical_CTES4<-function(tg_R1_lists,data_CORS_cancer,data.matrix,tg_R1_list_stat,cell_type_enrich_cut=0.45,tg_R1_cut=6,ext_cn=0,CT_balance=1)
{
        tg_selected_R4<-tg_R1_lists
        tg_R1_lists_st0<-tg_R1_lists
        for(i in 1:length(tg_R1_lists))
        {
          tg_R1_lists_st0[[i]]<-tg_R1_lists_st0[[i]][1:min(length(tg_R1_lists_st0[[i]]),tg_R1_cut)]
        }
        stat_selected_R4_RR<-compute_IM_stat(tg_list_c=tg_R1_lists_st0)
        stat_selected_R4_RR_max<-apply(stat_selected_R4_RR,1,max)
        names(tg_R1_lists_st0)<-names(stat_selected_R4_RR_max)

        tg_selected_R4<-tg_R1_lists_st0

        tg_genes_all<-c()
        for(i in 1:length(tg_selected_R4))
        {
                tg_genes_all<-c(tg_genes_all,tg_selected_R4[[i]])
        }
        tg_genes_all<-unique(tg_genes_all)
        data.matrix0 <- data.matrix
        data.matrix0_s<-data.matrix0[tg_genes_all,]
        data_23_s<-data_CORS_cancer[tg_genes_all,]

        BCV_stat_c<-BCV_ttest3(data_23_s,rounds=100,slice0=4,maxrank0=30,msep_cut=0.0001)
        dim_tt<-sum(BCV_stat_c[[1]]<0.01)
        print("Total Cell Dim")
        print(dim_tt)

        NN_ES_table_R4<-Explanation_BASE(data_23_s,tg_selected_R4)
        tg_ES_scores<-diag(NN_ES_table_R4)

        root_leaf<-compute_CompRowspace_NN_selflist(tg_data=data_23_s,tg_data_ng=data.matrix0_s,tg_list=tg_selected_R4,ROUNDS=3)
        names(root_leaf)
        root_leaf[["Leaf_CT"]]
        root_leaf[["Root_CT"]]
        root_leaf[["Other_leat_CT"]]
	  tg_possible_base_ids<-c()
	  tg_selected_R4_DFB<-list()
        if(length(root_leaf[["Leaf_CT"]])>0)
	  {
       	tg_possible_base_ids<-sort(root_leaf[["Leaf_CT"]])
        	tg_selected_R4_DFB<-select_R_base(tg_selected_R4,tg_possible_base_ids)
		STAT_c<-list()
        	for(i in 1:length(tg_possible_base_ids))
        	{
                STAT_c[[i]]<-tg_R1_list_stat[[tg_possible_base_ids[i]]]
        	}
        	R1_selectedCM_step2_results_HCTES<-Step2plus_Celltype_marker_inference_Hclust_CTES(tg_selected_R4_DFB,data_23_s,tg_R1_cut=tg_R1_cut,tg_R1_list_stat=STAT_c,cell_type_enrich_cut=0.5,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.9,IM_reso_level=IM_reso_level,extra_ctn=0)
        	CTES0<-c(R1_selectedCM_step2_results_HCTES[[1]],R1_selectedCM_step2_results_HCTES[[2]])
	  }
	  else
	  {
		CTES0<-list()
	  }
	  #c9_Hclust_CTES<-Compute_Rbase_SVD(data.matrix,CTES0)
        #d9<-cor(t(c9_Hclust_CTES),t(tProp))
        #apply(d9,2,max)
        #c9n_Hclust_CTES<-Compute_Rbase_SVD(data.matrix,tg_selected_R4_DFB)
        #d9n<-cor(t(c9n_Hclust_CTES),t(tProp))
        #apply(d9n,2,max)

        #############################################
        dim_diff<-dim_tt-length(CTES0)+ext_cn
        print(paste("old_version, dim_diff:", dim_diff, sep=""))
        ######dim_diff<-(dim_tt-length(CTES0))*((dim_tt-length(CTES0))>0)+ext_cn
        dim_diff_new <-  (dim_tt-length(CTES0))*((dim_tt-length(CTES0))>0)+ext_cn
        print(paste("new_version, dim_diff:", dim_diff_new, sep=""))

	  if(length(CTES0)>0)
	  {
        #par(mfcol=c(2,6))
        	ROUNDS<-3
        	tg_list<-tg_selected_R4
        	tg_selected_R4_SS<-CTES0
        	Rbase_selected_R4<-Compute_Rbase_SVD(data_23_s,tg_selected_R4_SS)
        	Base_c<-Rbase_selected_R4
        	aaa<-matrix(0,length(tg_list),ROUNDS)
        	bbb<-matrix(0,length(tg_list),ROUNDS)
        	ddd<-matrix(0,length(tg_list),ROUNDS)
        	for(j in 1:length(tg_list))
        	{
                tg_data_c<-data_23_s[tg_list[[j]],]
                ccc<-c()
                ccc<-compute_CompRowspace_NN_uni(tg_data_c,module_RB=Base_c,ROUNDS=3)
               # print(ccc)
                if(length(ccc)>0)
                {
                        tg_id1<-ccc[1,]
                        aaa[j,1:length(tg_id1)]<-tg_id1
                        bbb[j,1:length(tg_id1)]<-ccc[2,]
                        ddd[j,1:length(tg_id1)]<-ccc[3,]
                }
        	}
        	tg_ids_c<-intersect(order(-bbb[,3]),which(bbb[,3]>0.1))#################
        	tg_ids_c<-intersect(tg_ids_c,which(stat_selected_R4_RR_max>cell_type_enrich_cut))#############
		#hist(bbb[tg_ids_c,3])
        	aaa<-names(tg_list)[tg_ids_c]
        	stat_c<-table(aaa)
        	aaa0<-rep(0,length(aaa))
        	names(aaa0)<-aaa
        	aaa1<-unique(aaa)
        	sss_c<-rep(1,length(aaa1))
        	names(sss_c)<-aaa1
        	sss_b<-1/(table(names(tg_selected_R4_SS))+1)
        	sss_c[intersect(names(sss_c),names(sss_b))]<-sss_b[intersect(names(sss_c),names(sss_b))]
        	tg_scores_c<-bbb[tg_ids_c,3]
        	for(i in 1:length(tg_ids_c))
        	{
                tg_scores_c[i]<-tg_scores_c[i]+sss_c[names(tg_list)[tg_ids_c[i]]]
        	}
       	tg_ids_cs<-tg_ids_c[which(tg_scores_c==max(tg_scores_c))][1]
        }
	  else
	  {
		tg_list<-tg_selected_R4
		tg_selected_R4_SS<-CTES0
		tg_ids_c<-intersect(1:length(tg_list),which(stat_selected_R4_RR_max>cell_type_enrich_cut))
		#hist(bbb[tg_ids_c,3])
        	aaa<-names(tg_list)[tg_ids_c]
        	stat_c<-table(aaa)
        	aaa0<-rep(0,length(aaa))
        	names(aaa0)<-aaa
        	aaa1<-unique(aaa)
        	sss_c<-rep(1,length(aaa1))
        	names(sss_c)<-aaa1
        	sss_b<-1/(table(names(tg_selected_R4_SS))+1)
        	sss_c[intersect(names(sss_c),names(sss_b))]<-sss_b[intersect(names(sss_c),names(sss_b))]
        	tg_scores_c<--tg_ES_scores
        	for(i in 1:length(tg_ids_c))
        	{
                tg_scores_c[i]<-tg_scores_c[i]+sss_c[names(tg_list)[tg_ids_c[i]]]
        	}
       	tg_ids_cs<-tg_ids_c[which(tg_scores_c==max(tg_scores_c))][1]
        	bbb<-matrix(1,length(tg_list),ROUNDS)
        }	
        #hist(tg_scores_c)
        NN<-0
        node_add<-c()
        while((NN<dim_diff)&(max(bbb[,3]>0.1))&(length(tg_ids_c)>0))
        {
                NN<-NN+1
                ii<-length(tg_selected_R4_SS)+1
                tg_selected_R4_SS[[ii]]<-tg_list[[tg_ids_cs]]
                names(tg_selected_R4_SS)[ii]<-names(tg_list)[tg_ids_cs]
                Rbase_selected_R4<-Compute_Rbase_SVD(data_23_s,tg_selected_R4_SS)
                node_add<-c(node_add,tg_ids_cs)
                Base_c<-Rbase_selected_R4
                aaa<-matrix(0,length(tg_list),ROUNDS)
                bbb<-matrix(0,length(tg_list),ROUNDS)
                ddd<-matrix(0,length(tg_list),ROUNDS)
        	for(j in 1:length(tg_list))
        	{
                tg_data_c<-data_23_s[tg_list[[j]],]
                ccc<-c()
                ccc<-compute_CompRowspace_NN_uni(tg_data_c,module_RB=Base_c,ROUNDS=3)
               # print(ccc)
                if(length(ccc)>0)
                {
                        tg_id1<-ccc[1,]
                        aaa[j,1:length(tg_id1)]<-tg_id1
                        bbb[j,1:length(tg_id1)]<-ccc[2,]
                        ddd[j,1:length(tg_id1)]<-ccc[3,]
                }
                 }
                tg_ids_c<-intersect(order(-bbb[,3]),which(bbb[,3]>0.1))#################
                tg_ids_c<-intersect(tg_ids_c,which(stat_selected_R4_RR_max>0.25))#############

                #hist(bbb[tg_ids_c,3])
                aaa<-names(tg_list)[tg_ids_c]
                stat_c<-table(aaa)
                aaa0<-rep(0,length(aaa))
                names(aaa0)<-aaa
                aaa1<-unique(aaa)
                sss_c<-rep(1,length(aaa1))
                names(sss_c)<-aaa1
                sss_b<-1/(table(names(tg_selected_R4_SS))+CT_balance)
                sss_c[intersect(names(sss_c),names(sss_b))]<-sss_b[intersect(names(sss_c),names(sss_b))]
                tg_scores_c<-bbb[tg_ids_c,3]
                for(i in 1:length(tg_ids_c))
                {
                        tg_scores_c[i]<-tg_scores_c[i]+sss_c[names(tg_list)[tg_ids_c[i]]]
                }
                tg_ids_cs<-tg_ids_c[which(tg_scores_c==max(tg_scores_c))][1]
                #hist(tg_scores_c)
        }
        #c12<-Compute_Rbase_SVD(data.matrix,tg_selected_R4_SS)
        #d12<-cor(t(c12),t(tProp))
        #apply(d12,2,max)
        ccc<-list(tg_selected_R4_SS,root_leaf,R1_selectedCM_step2_results_HCTES,node_add)
        names(ccc)<-c("Selected_Bases","First_step_leaf_nodes","Second_step_filtered_leaf_nodes","Further_added_other_leaves")
        return(ccc)
}