
ICTD <- function(data_bulk)
{
  data.matrix = data_bulk

  ###########################################################
  if(length(colnames(data.matrix)) == 0) {
  	warning("input data do NOT have colnames")
  	colnames(data.matrix) <- paste( "Setsample", 1:ncol(data.matrix), sep="")
  }
  data.matrix <- rm_zero_row(data.matrix)

  ###########  RNAseq log(X+1)
  if(max(data.matrix) > 20){
  	d.matrix<-log(data.matrix + 1)
  	d.matrix<-as.matrix(d.matrix)
  }else{
  	d.matrix <- as.matrix(data.matrix)
  }

  data0<-d.matrix

  ########### prepare input data for cancer module identification
  data2<-data0

  ###########take the highest expressed probe for the genes with duplicated probes
  data21<-data2[order(-apply(data2,1,mean)),]
  data22<-data21[unique(rownames(data21)),]
  data22_code<-data22[intersect(rownames(data22),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
  data23<-data22[intersect(rownames(data22),TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,4]),]
  data23<-normalize_data2(data23)		#BUG source: get NaN value after "normalize_data2" because some zero rows
  ##############

  #data_ccc<-data_CORS_cancer
  data_CORS_cancer <- data23
  data_ccc <- data23
  list_c1<-MRHCA_IM_compute_MR(data_CORS_cancer=data_ccc,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)	#MR_1_rank_test_function_v1.1
  MR_IM_result_new_c<-MRHCA_IM_compute_full_pub_new(data_CORS_cancer=data_ccc,list_c=list_c1,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
  tg_key="nonono"
  list_new_c2<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key,cell_type_enrich_cut=0.4,cor_cut0=0.7,num_cut=7,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
  R1_filter_step1_results_new<-R1_list_filtering_step1_new(list_new_c2,data_CORS_cancer=data_ccc,max_cut=7,cutn0=7,cut10=0.7,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)

  tg_R1_lists<-R1_filter_step1_results_new[[4]]
  #print(length(tg_R1_lists))
  #R1_selectedCM_step2_results_new<-Step2plus_Celltype_marker_inference(tg_R1_lists,data_CORS_cancer,tg_R1_cut=6,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],resolution_level=resolution_level0)
  #R1_selectedCM_step2_results_new <- Step2plus_Celltype_marker_inference_new(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=20)
  #comb <- c(R1_selectedCM_step2_results_new[[1]],R1_selectedCM_step2_results_new[[2]])

  #R1_selectedCM_step2_results_Hclust<-Step2plus_Celltype_marker_inference_Hclust(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=40)
  #R1_selectedCM_step2_results_fixedCT<-Step2plus_Celltype_marker_inference_fixedCT(tg_R1_lists,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]])
  #R1_selectedCM_step2_results_HP<-Step2plus_Celltype_marker_inference_Hclust_plus(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.5,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.85,IM_reso_level=IM_reso_level)
  #R1_selectedCM_step2_results_HCTES<-Step2plus_Celltype_marker_inference_Hclust_CTES(tg_R1_lists,data_CORS_cancer,tg_R1_cut=5,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.6,resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.95,IM_reso_level=IM_reso_level,extra_ctn=3)
  #CTES<-c(R1_selectedCM_step2_results_HCTES[[1]],R1_selectedCM_step2_results_HCTES[[2]])

  #R1_selectedCM_step2_results_Hierarchical_CTES<-Step2plus_Celltype_marker_Hierarchical_CTES(tg_R1_lists=tg_R1_lists,data_CORS_cancer=data_CORS_cancer,data.matrix=data.matrix,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.45,tg_R1_cut=6)
  R1_selectedCM_step2_results_Hierarchical_CTES<- Step2plus_Celltype_marker_Hierarchical_CTES4(tg_R1_lists,data_CORS_cancer,data.matrix,tg_R1_list_stat=R1_filter_step1_results_new[[3]][[2]],cell_type_enrich_cut=0.45,tg_R1_cut=6,ext_cn=10,CT_balance=1)  #72056_ext=6, ctBALANCE=1
  CTES3<-R1_selectedCM_step2_results_Hierarchical_CTES[[1]]


  #c9_Hclust_CTES<-Compute_Rbase_SVD(data.matrix,tg_list)
  #d9<-cor(t(c9_Hclust_CTES),t(tProp))
  #apply(d9,2,max)

  ##############################
  #Leaf node selection
  ##############################
  tg_selected_R4<-CTES3	#!!!!!!!!!!!
  tg_genes_all<-c()
  for(i in 1:length(tg_selected_R4))
  {
  	tg_genes_all<-c(tg_genes_all,tg_selected_R4[[i]])
  }
  tg_genes_all<-unique(tg_genes_all)
  data.matrix0 <- data.matrix
  data.matrix0_s<-data.matrix0[tg_genes_all,]
  data23 = data_CORS_cancer
  data_23_s<-data23[tg_genes_all,]

  root_leaf<-compute_CompRowspace_NN_selflist(tg_data=data_23_s,tg_data_ng=data.matrix0_s,tg_list=CTES3,ROUNDS=3)
  names(root_leaf)
  root_leaf[["Leaf_CT"]]
  root_leaf[["Root_CT"]]
  root_leaf[["Other_leat_CT"]]
  tg_possible_base_ids<-sort(c(root_leaf[["Leaf_CT"]],root_leaf[["Other_leat_CT"]]))
  NMF_selected_R1<-select_R_base(CTES3,tg_possible_base_ids)		#just like tg_selected_R4_RR
  c10_Hclust_CTES<-Compute_Rbase_SVD(data.matrix,NMF_selected_R1)
  Prop <- c10_Hclust_CTES
  colnames(Prop) <- colnames(data.matrix)

  #ICTD
  #qnmf_result<-NMF_method1_test_version3_pkg(tg_list=NMF_selected_R1,data_ng=data.matrix0_s,data_normalized=data_23_s, max_ES_cut=0.3,NMF_RR=10)
  #NMF_indi_all = qnmf_result[[3]][["NMF_indi_all"]]
  #X1 = qnmf_result[[1]][["X1"]]
  #V = qnmf_result[[1]][["V"]]


  #return(list(ictd.proportion=V,ictd.marker.R1=tg_R1_lists))
  return(list(ictd.proportion=Prop,ictd.marker=NMF_selected_R1))

}

