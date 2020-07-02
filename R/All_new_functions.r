
MRHCA_IM_compute_full_pub_new<-function(data_CORS_cancer,list_c,IM_id_list,immune_cell_uni_table=marker_stats20_uni,step_size0=20)
{
        MR_M<-list_c[[1]]
        cor_c<-list_c[[2]]
        down_lim<-list_c[[3]]
        print("Compute MR IM genes")
        MR_IM_result_c<-list()
        print(nrow(MR_M))
        step_size<-step_size0
        IM_id_list0<-IM_id_list
        MR_IM_result_c<-c()
      for(ii in 1:nrow(list_c[[1]]))
      {
                tg_gene<-rownames(MR_M)[ii]
                x0<-sqrt(sort(MR_M[tg_gene,]))
                tg_growth_rate<-calculate_growth_rate2(x0,step=step_size)
                tg_ccc1<-which(tg_growth_rate<down_lim)
                aaa<-sqrt(sort(MR_M[tg_gene,]))
                bbb3<-aaa/(1:length(aaa))
                bbb4<-cor_c[tg_gene,names(sort(MR_M[tg_gene,]))]

                ddd<-cbind(1:length(tg_growth_rate),tg_growth_rate,bbb3,bbb4)
                colnames(ddd)[1:4]<-c("MR_order_ID","growth_rate","sorted_MR","Corr")
                tg_ccc3<-intersect(tg_ccc1,which(bbb4>0.7))

                fff<-""
                if(length(tg_ccc3)>1)
                {
                        fff<-as.matrix(ddd[1:max(tg_ccc3),])
                }
                if(ii%%500==1)
                {
                        print(ii)
                }
                MR_IM_result_c[[ii]]<-list(fff,tg_ccc3)
      }
      names(MR_IM_result_c)<-rownames(MR_M)[1:length(MR_IM_result_c)]
      return(MR_IM_result_c)
}




clean_rank1_module_new<-function(data_c,module_info,module_rank,st0=8,RR=50)
{
        N<-0
        nc<-c()
        module_new<-list()
        for(i in 1:length(module_info))
        {
                if(module_rank[i]==1)
                {
                        N<-N+1
                        nc<-c(nc,names(module_info)[i])
                        module_new[[N]]<-module_info[[i]]
                }
                if(module_rank[i]>1)
                {
                        ccc<-module_info[[i]]
                        st<-st0
                        rr<-1
                        while((rr==1)&(st<=length(ccc)))
                        {
                                tg_genes<-c(ccc[1:st])
                                pp<-BCV_ttest2(data_c[tg_genes,],rounds=RR,maxrank0=5)
                                rr<-sum(pp<0.001)
                                st<-st+1
                        }
                        tg_genes<-tg_genes[-length(tg_genes)]
                        if(length(tg_genes)>st0)
                        {
                                N<-N+1
                                nc<-c(nc,names(module_info)[i])
                                module_new[[N]]<-tg_genes
                        }
                }
        }
        names(module_new)<-nc
        return(module_new)
}


R1_list_filtering_step1_new<-function(list_c2,data_CORS_cancer,max_cut=20,cutn0=20,cut10=0.8,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
{
tg_1_rank_markers0<-list_c2[[1]]
tg_m_names0<-list_c2[[2]]
RMSE_on_CORS_cancer_c<-RMSE_row(data_CORS_cancer)
tg_marker_lists<-list()
for(i in 1:length(list_c2[[1]]))
{
        ccc<-c(names(list_c2[[1]])[i],rownames(list_c2[[1]][[i]]))
        tg_marker_lists[[i]]<-ccc
}
names(tg_marker_lists)<-names(list_c2[[1]])

pp_all<-c()
for(i in 1:length(tg_marker_lists))
{
        pp<-sum(BCV_ttest2(data_CORS_cancer[tg_marker_lists[[i]],],maxrank0=20)<0.001, na.rm=T)
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f1<-clean_rank1_module_new(data_CORS_cancer,tg_marker_lists,pp_all,st0=6)
pp_R1_marker_list_f1.5<-cut_modules(pp_R1_marker_list_f1,cutn=cutn0)

stat_p1<-list()
for(i in 1:length(pp_R1_marker_list_f1.5))
{
	tg_gene_c<-pp_R1_marker_list_f1.5[[i]][-1]
	stat_p1[[i]]<-tg_1_rank_markers0[[names(pp_R1_marker_list_f1.5)[i]]][tg_gene_c,]
}

names(tg_marker_lists)<-names(pp_R1_marker_list_f1.5)
tg_genes_all<-c()
for(i in 1:length(pp_R1_marker_list_f1.5))
{
        tg_genes_all<-c(tg_genes_all,pp_R1_marker_list_f1.5[[i]])
}
tg_genes_all<-unique(sort(tg_genes_all))

print("filter 1 and stat done!")
#ccc<-MAP_GPL570_genes2(R1_markers_f1)
#names(ccc)<-names(Top_cell_proportion)
#table(names(Top_cell_proportion))
R1_markers_f1<-pp_R1_marker_list_f1.5

cut1<-cut10
ccc<-compute_min_jaccard(R1_markers_f1)
ccc0<-ccc>cut1

stat_cc<-c(1:nrow(ccc0))
names(stat_cc)<-1:nrow(ccc0)
for(i in 1:nrow(ccc0))
{
        for(j in 1:ncol(ccc0))
        {
                if((i<j)&(ccc0[i,j]>0))
                {
                        nn<-max(i,j)
                        stat_cc[which(stat_cc==i)]<-nn
                        stat_cc[which(stat_cc==j)]<-nn
                }
        }
}
table(stat_cc)
tg_ccc<-unique(stat_cc)
R1_marker_list_f2<-list()
N<-0
for(i in 1:length(tg_ccc))
{
        N<-N+1
        tg_ids<-as.numeric(names(stat_cc)[which(stat_cc==tg_ccc[i])])
        ccc<-c()
        for(j in 1:length(tg_ids))
        {
                ccc<-c(ccc,R1_markers_f1[[tg_ids[j]]])
        }
        ccc<-unique(ccc)
        R1_marker_list_f2[[N]]<-ccc
}

R1_marker_list_f2.5_stat<-rank_based_module_sorting(data_CORS_cancer,R1_marker_list_f2,IM_id_list,immune_cell_uni_table=immune_cell_uni_table)
R1_marker_list_f2.5<-R1_marker_list_f2.5_stat[[1]]

pp_all<-c()
for(i in 1:length(R1_marker_list_f2.5))
{
        pp<-sum(BCV_ttest(data_CORS_cancer[R1_marker_list_f2.5[[i]],],maxrank0=20)<0.001, na.rm=T)
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-clean_rank1_module(data_CORS_cancer,R1_marker_list_f2.5,pp_all,st0=6)
pp_R1_marker_list_f3.5<-cut_modules(pp_R1_marker_list_f3,cutn=cutn0)
R1_marker_list_f3.5_stat<-rank_based_module_sorting(data_CORS_cancer,pp_R1_marker_list_f3.5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table)

print("filter 2 done!")
ccc<-c()
nn<-c()
for(i in 1:length(pp_R1_marker_list_f3.5))
{
	ccc0<-c()
	for(j in 1:length(IM_id_list))
	{
      if(length(IM_id_list[[j]])>1)
      {
            cc0<-apply(immune_cell_uni_table[pp_R1_marker_list_f3.5[[i]],IM_id_list[[j]]],1,sum)/sum((1/(1:length(IM_id_list[[j]]))))
      }
      else
      {
           cc0<-immune_cell_uni_table[pp_R1_marker_list_f3.5[[i]],IM_id_list[[j]]]
      }
		ccc0<-cbind(ccc0,cc0)
	}
	colnames(ccc0)<-names(IM_id_list)
      ddd<-apply(ccc0,2,mean)
      ccc<-rbind(ccc,ddd)
      nn<-c(nn,colnames(ccc0)[which(ddd==max(ddd))[1]])
}
rownames(ccc)<-nn
cell_enrich_stat<-ccc
names(pp_R1_marker_list_f3.5)<-nn
rrr<-rep(1,length(pp_R1_marker_list_f3.5))
Filter_1_result_list<-list(pp_R1_marker_list_f1,R1_marker_list_f2,R1_marker_list_f3.5_stat,pp_R1_marker_list_f3.5,rrr,cell_enrich_stat)
names(Filter_1_result_list)<-c("R1_marker_list_f1","R1_marker_list_f2","R1_marker_list_f3.5_stat","R1_marker_list_f3.5","R1_marker_list_rank","R1_marker_list_f3.5_cell_enrich_stat")
return(Filter_1_result_list)
}



rank_based_module_sorting<-function(data_c,tg_list,IM_id_list,immune_cell_uni_table=marker_stats20_uni)
{
	tg_list_new<-list()
	tg_list_stat<-list()
	for(i in 1:length(tg_list))
	{
		rr<-svd(data_c[tg_list[[i]],])$v[,1]
		if(mean(cor(t(data_c[tg_list[[i]],]),rr))<0)
		{
			rr<--rr
		}
		tg_genes_cc<-tg_list[[i]][order(-cor(t(data_c[tg_list[[i]],]),rr))]

		ccc<-immune_cell_uni_table[tg_genes_cc,]^(1/2)
		ccc0<-c()
		for(j in 1:length(IM_id_list))
		{
                       			 if(length(IM_id_list[[j]])>1)
                       			 {
                               		 xx<-sum((1/(1:length(IM_id_list[[j]])))^(1/2))
                               		 cc0<-(cumsum(apply(ccc[,IM_id_list[[j]]],1,sum))/xx/c(1:nrow(ccc)))
                       			 }
                        		else
                        		{
                                		cc0<-(cumsum(ccc[,IM_id_list[[j]]])/c(1:nrow(ccc)))
                        		}
                       			 ccc0<-cbind(ccc0,cc0)
		}
		colnames(ccc0)<-names(IM_id_list)
		tg_list_new[[i]]<-tg_genes_cc
		tg_list_stat[[i]]<-ccc0
	}
	return(list(tg_list_new,tg_list_stat))
}

Process_MR_IM_result_new<-function(MR_IM_result_c=MR_IM_result_c,tg_key_c=tg_key_c,cor_cut0=0.7,cell_type_enrich_cut=0.5,num_cut=5,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
{
tg_1_rank_markers<-list()
tg_m_names<-c()
N<-0
print("Select Marker!")
for(i in 1:length(MR_IM_result_c))
{
        if(length(MR_IM_result_c[[i]][[2]])>=num_cut)
        {
                ss<-scoring_MR_order(MR_IM_result_c[[i]][[2]])
                if(ss>=num_cut)
                {
                        ccc<-MR_IM_result_c[[i]][[1]][1:ss,]
                        if(sum(ccc[,4]>cor_cut0)>=num_cut)
                        {
                                tg_ccc<-names(which(ccc[,4]>cor_cut0))
                                if(length(tg_ccc)>=num_cut)
                                {
						    ccc1<-ccc[order(-ccc[,4]),]
						    tg_gene<-rownames(ccc1)
						   ccc<-immune_cell_uni_table[tg_gene,]^(1/2)
                					ccc0<-c()
               				 for(j in 1:length(IM_id_list))
                				{
                       			 if(length(IM_id_list[[j]])>1)
                       			 {
                               		 xx<-sum((1/(1:length(IM_id_list[[j]])))^(1/2))
                               		 cc0<-(cumsum(apply(ccc[,IM_id_list[[j]]],1,sum))/xx/c(1:nrow(ccc)))
                       			 }
                        		else
                        		{
                                		cc0<-(cumsum(ccc[,IM_id_list[[j]]])/c(1:nrow(ccc)))
                        		}
                       			 ccc0<-cbind(ccc0,cc0)
                				}
               				 colnames(ccc0)<-names(IM_id_list)
						ccc3<-cbind(ccc1,ccc0)
                                        ddd<-apply(ccc0,2,mean)
                                        if(max(ddd)>cell_type_enrich_cut)
						    {
								tg_c_ids<-names(which(ddd==max(ddd)))[1]
								tg_c_id0<-IM_id_list[[tg_c_ids]]
								eee<-ccc[,tg_c_id0]
								if(length(tg_c_id0)>1)
								{
									eee<-apply(ccc[,tg_c_id0],1,mean)
								}

                                        	if(sum(eee>0.5)>=num_cut2)
                                        	{
                                                N<-N+1
                                                tg_1_rank_markers[[N]]<-ccc3[which(ccc3[,4]>cor_cut0),]
                                                tg_m_names<-c(tg_m_names,names(MR_IM_result_c)[i])
                                      	  }
                                        }
                                }
                        }
                }
        }
}
print("Select Marker Done!")
#tg_RF2<-paste(tg_key_c,"_1rankmarker_cell_type_consistency.pdf",sep="")
##library(gplots)
#colors = c(0:100)/100
#my_palette <- grDevices::colorRampPalette(c("white","white", "blue"))(n =100)
#pdf(tg_RF2)
#for(i in 1:length(tg_1_rank_markers))
#{
#       aaa<-tg_1_rank_markers[[i]][,-c(1:4)]
#       heatmap.2(aaa,Rowv=F,Colv =F,scale="none",main=tg_m_names[i],
#       col=my_palette,breaks=colors,density.info="none",dendrogram="both",
#       trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)
#}
#dev.off()
names(tg_1_rank_markers)<-tg_m_names
list_cc<-list(tg_1_rank_markers,tg_m_names)
return(list_cc)
}




Process_MR_IM_result_GPL570_new<-function(MR_IM_result_c=MR_IM_result_c,tg_key_c=tg_key_c,cor_cut0=0.7,cell_type_enrich_cut=0.5,num_cut=5,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
{
tg_1_rank_markers<-list()
tg_m_names<-c()
N<-0
print("Select Marker!")
for(i in 1:length(MR_IM_result_c))
{
        if(length(MR_IM_result_c[[i]][[2]])>=num_cut)
        {
                ss<-scoring_MR_order(MR_IM_result_c[[i]][[2]])
                if(ss>=num_cut)
                {
                        ccc<-MR_IM_result_c[[i]][[1]][1:ss,]
                        if(sum(ccc[,4]>cor_cut0)>=num_cut)
                        {
                                tg_ccc<-names(which(ccc[,4]>cor_cut0))
			              tg_ccc2<-unique(GPL570_id_symbol0[intersect(GPL570_id_symbol[,1],c(names(MR_IM_result_c)[i],tg_ccc)),2])
                                if(length(tg_ccc2)>=num_cut)
                                {
						    ccc1<-ccc[order(-ccc[,4]),]
						    tg_gene<-rownames(ccc1)
						   ccc<-immune_cell_uni_table[tg_gene,]^(1/2)
                					ccc0<-c()
               				 for(j in 1:length(IM_id_list))
                				{
                       			 if(length(IM_id_list[[j]])>1)
                       			 {
                               		 xx<-sum((1/(1:length(IM_id_list[[j]])))^(1/2))
                               		 cc0<-(cumsum(apply(ccc[,IM_id_list[[j]]],1,sum))/xx/c(1:nrow(ccc)))
                       			 }
                        		else
                        		{
                                		cc0<-(cumsum(ccc[,IM_id_list[[j]]])/c(1:nrow(ccc)))
                        		}
                       			 ccc0<-cbind(ccc0,cc0)
                				}
               				 colnames(ccc0)<-names(IM_id_list)
						ccc3<-cbind(ccc1,ccc0)
                                        ddd<-apply(ccc0,2,mean)
                                        if(max(ddd)>cell_type_enrich_cut)
						    {
								tg_c_ids<-names(which(ddd==max(ddd)))[1]
								tg_c_id0<-IM_id_list[[tg_c_ids]]
								eee<-ccc[,tg_c_id0]
								if(length(tg_c_id0)>1)
								{
									eee<-apply(ccc[,tg_c_id0],1,mean)
								}
                                        }
                                        if(sum(eee>0.5)>=num_cut2)
                                        {
                                                N<-N+1
                                                tg_1_rank_markers[[N]]<-ccc3[which(ccc3[,4]>cor_cut0),]
                                                tg_m_names<-c(tg_m_names,names(MR_IM_result_c)[i])
                                        }
                                }
                        }
                }
        }
}
print("Select Marker Done!")
#tg_RF2<-paste(tg_key_c,"_1rankmarker_cell_type_consistency.pdf",sep="")
##library(gplots)
#colors = c(0:100)/100
#my_palette <- grDevices::colorRampPalette(c("white","white", "blue"))(n =100)
#pdf(tg_RF2)
#for(i in 1:length(tg_1_rank_markers))
#{
#       aaa<-tg_1_rank_markers[[i]][,-c(1:4)]
#       heatmap.2(aaa,Rowv=F,Colv =F,scale="none",main=tg_m_names[i],
#       col=my_palette,breaks=colors,density.info="none",dendrogram="both",
#       trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)
#}
#dev.off()
names(tg_1_rank_markers)<-tg_m_names
list_cc<-list(tg_1_rank_markers,tg_m_names)
return(list_cc)
}


