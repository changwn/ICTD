#load
#load("Functions8.1.RData")
#load("TCGA_ensem_annotation.RData")
#load("IM_markers_20190302.RData")
#load("Key_data_TCGA_IM_1rank_marker_identification.RData")
#load("Cancer_module_common_genes.RData")
#source("MR_functions_v1.0.r")
#source("key_functions_02122018.r")
#source("MR_1_rank_test_functions_v1.1.r")

#source("./estimateScore_0.1.r")
#source("./coexpression_WGCNA_0.2.R")
#source("./2018_05_05_find_cancer_module.r")

#library(gplots)
#library(bcv)
#library(igraph)
#library(estimate)
#library(WGCNA)
#library(tictoc)  #to use timer

.onLoad <- function(libname, pkgname) {
	utils::data(IM_markers_20190302, package = pkgname, envir = parent.env(environment()))
	marker_stats1_uni <- ICTD::marker_stats1_uni
	marker_stats20_uni<-marker_stats1_uni
	
	utils::data(TCGA_ensem_annotation, package = pkgname, envir = parent.env(environment()))	#load for pkg
	TCGA_ensem_annotation <- ICTD::TCGA_ensem_annotation	
	
	utils::data(Cancer_module_common_genes, package = pkgname, envir = parent.env(environment()))	#load for pkg
	Cancer_module_common_genes <- ICTD::Cancer_module_common_genes	
	
	
	
	#############
	IM_id_list<-list()
  for(i in 1:ncol(marker_stats20_uni))
  {
  	IM_id_list[[i]]<-i
  }
  colnames(marker_stats20_uni)
  IM_id_list[[c(i+1)]]<-c(3,4)
  IM_id_list[[c(i+2)]]<-c(3,4,12)
  IM_id_list[[c(i+3)]]<-c(2,3,4,12)
  IM_id_list[[c(i+4)]]<-c(8,9)
  IM_id_list[[c(i+5)]]<-c(8,9,11)
  names(IM_id_list)<-c(colnames(marker_stats20_uni),"T","TNK","Lymp","MM","Myeloid")

  marker_stats_p<-marker_stats0[,colnames(marker_stats20_uni)]
  immune_cell_uni_table0<-1/marker_stats_p
  immune_cell_uni_table0[which(marker_stats_p==0)]<-0
  immune_cell_uni_table0_GPL570<-immune_cell_uni_table0

  ###########################
  immune_cell_uni_table0_GS<-marker_stats20_uni

	utils::data(Functions8.1, package = pkgname, envir = parent.env(environment()))
	GPL570_id_symbol0 <- ICTD::GPL570_id_symbol0
  #############################
  ccc<-marker_stats0[intersect(rownames(marker_stats0),GPL570_id_symbol0[,1]),colnames(marker_stats_p)]
  rownames(ccc)<-GPL570_id_symbol0[rownames(ccc),2]
  aaa<-unique(rownames(ccc))
  ddd<-c()
  for(j in 1:length(aaa))
  {
  	if(sum(rownames(ccc)==aaa[j])==1)
  	{
  		bbb<-ccc[aaa[j],]
  		bbb[which(bbb==0)]<-100
  		bbb<-1/bbb*(bbb<=3)
  		ddd<-rbind(ddd,bbb)
  	}
  	else
  	{
  		bbb<-ccc[which(rownames(ccc)==aaa[j]),]
  		bbb[which(bbb==0)]<-100
  		bbb<-1/bbb*(bbb<=3)
  		ddd<-rbind(ddd,apply(bbb,2,mean))
  	}
  }
  rownames(ddd)<-aaa
  marker_stats20_uni_GPL570_GS<-ddd
  immune_cell_uni_table0_GPL570_GS<-marker_stats20_uni_GPL570_GS

  ################################
  #load("TCGA_ensem_annotation.RData")
						
  TCGA_ensem_annotation0<-TCGA_ensem_annotation
  rownames(TCGA_ensem_annotation0)<-TCGA_ensem_annotation0[,2]
  tg_uni_genes<-names(which(table(TCGA_ensem_annotation[,2])==1))
  tg_uni_genes0<-names(which(TCGA_ensem_annotation0[tg_uni_genes,3]=="protein_coding"))
  rownames(TCGA_ensem_annotation)<-TCGA_ensem_annotation[,7]

  ##############
  tg_top_markers<-list()
  for(i in 1:ncol(marker_stats20_uni))
  {
  	tg_top_markers[[i]]<-names(which(marker_stats20_uni[,i]==1))
  }
  names(tg_top_markers)<-colnames(marker_stats20_uni)
  tg_top_markers_GS<-tg_top_markers

  tg_top_markers<-list()
  for(i in 1:ncol(marker_stats_p))
  {
  	tg_top_markers[[i]]<-names(which(marker_stats_p[,i]==1))
  }
  names(tg_top_markers)<-colnames(marker_stats_p)
  tg_top_markers_GPL570<-tg_top_markers

  ccc<-marker_stats0[intersect(rownames(marker_stats0),GPL570_id_symbol0[,1]),]
  rownames(ccc)<-GPL570_id_symbol0[rownames(ccc),2]
  tg_top_markers<-list()
  for(i in 1:ncol(ccc))
  {
  	tg_top_markers[[i]]<-names(which(ccc[,i]==1))
  }
  names(tg_top_markers)<-colnames(ccc)
  tg_top_markers_GPL570_GS<-tg_top_markers


  ############################
  Cancer_module_common_genes_GPL570<-extract_data_symbol(GPL570_id_symbol,Cancer_module_common_genes)[,1]

  colors = c(0:100)/100
  my_palette <- grDevices::colorRampPalette(c("white", "blue"))(n =100)


  assign("IM_id_list", IM_id_list, envir = parent.env(environment()))
  assign("immune_cell_uni_table0_GS", immune_cell_uni_table0_GS, envir = parent.env(environment()))

  
  .onLoad_temp_func3()
}



