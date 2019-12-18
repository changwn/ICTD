
RMSE_row<-function(data0)
{
        return(apply(data0^2,1,mean))
}

normalize_data2<-function(data0)
{
        data1<-data0
        for(i in 1:nrow(data1))
        {
                data1[i,]<-(data1[i,]-mean(data1[i,]))/sd(data1[i,])
        }
        return(data1)
}


Compute_base_ES_score<-function(ES_c,tg_selected_c)
{
	ccc<-c()
	for(i in 1:length(tg_selected_c))
	{	
		ccc<-c(ccc,mean(ES_c[tg_selected_c[[i]]]))
	}
	names(ccc)<-names(tg_selected_c)
	return(ccc)
}

compute_CompRowspace_NN_selflist<-function(tg_data,tg_data_ng,tg_list,ROUNDS=3)
{
	Rbase_selected_R4_og<-Compute_Rbase_SVD(tg_data_ng,tg_list)
	Rbase_selected_R4<-Compute_Rbase_SVD(tg_data,tg_list)
	aaa<-matrix(0,length(tg_list),ROUNDS)
	bbb<-matrix(0,length(tg_list),ROUNDS)
	ddd<-matrix(0,length(tg_list),ROUNDS)
      for(j in 1:length(tg_list))
      {
		tg_data_c<-tg_data[tg_list[[j]],]
        	Base_c<-Rbase_selected_R4[-j,]
		#print(tg_data_c[,1:5])
		#print(j)
		#print(Base_c[,1:5])
		#print(j)
		tg_id_kk<-setdiff(1:length(tg_list),j)
		ccc<-c()
		ccc<-compute_CompRowspace_NN_uni(tg_data_c,module_RB=Base_c,ROUNDS=3)
		#print(ccc)
		if(length(ccc)>0)
		{
			tg_id1<-tg_id_kk[ccc[1,]]
			aaa[j,1:length(tg_id1)]<-tg_id1
			bbb[j,1:length(tg_id1)]<-ccc[2,]
			ddd[j,1:length(tg_id1)]<-ccc[3,]
		}
  	}
	tg_ids<-which(apply(bbb,1,min)<0.2)
	c_cell_depend<-list()
	N<-0
	nn<-c()
	for(j in 1:length(tg_ids))
	{
		tg_c<-aaa[tg_ids[j],which(ddd[tg_ids[j],]>0.1)]
		if(length(tg_c)>0)
		{
			data_ccc<-t(Rbase_selected_R4_og[c(tg_ids[j],tg_c),])
			colnames(data_ccc)<-c("Y",tg_c)
			rownames(data_ccc)<-1:nrow(data_ccc)
			data_ccc<-as.data.frame(data_ccc)
			lm_cc<-lm(Y~.+0,data=data_ccc)
			ccc<-summary(lm_cc)$coefficients
			cc1<-which((ccc[,4]<0.01)&ccc[,1]>0)
			cc2<-which((ccc[,4]<0.01)&ccc[,1]<0)
			dd1<-c()
			dd2<-c()
			if(length(cc1)>0)
			{
				dd1<-tg_c[cc1]
			}
			if(length(cc2)>0)
			{
				dd2<-tg_c[cc2]
			}
			dd3<-setdiff(tg_c,c(dd1,dd2))
			c_l<-c(tg_ids[j],dd2)
			c_r<-c(dd1)
			cc<-list(c_l,c_r)
			names(cc)<-c("right","left")
			N<-N+1
			c_cell_depend[[N]]<-cc
			nn<-c(nn,tg_ids[j])
		}
	}
	names(c_cell_depend)<-nn
	baba<-c()
	didi<-c()
	for(i in 1:length(c_cell_depend))
	{
		cc<-c_cell_depend[[i]]
		if((length(cc[[1]])==1)&(length(cc[[2]])>1))
		{
			baba<-c(baba,cc[[1]])
			didi<-c(didi,cc[[2]])
		}
	}
	baba<-unique(sort(baba))
	didi<-setdiff(unique(sort(didi)),baba)
	daye<-setdiff(1:length(tg_list),unique(sort(c(didi,baba))))
	ccc<-list(didi,baba,daye,c_cell_depend)
	names(ccc)<-c("Leaf_CT","Root_CT","Other_leat_CT","Cell_dependency")
      return(ccc)
}

compute_CompRowspace_NN_uni<-function(tg_data_c,module_RB=Base_c,ROUNDS=3)
{
          tg_list_c<-list()
          tg_list_c[[1]]<-rownames(tg_data_c)
          dd<-Compute_Rbase_SVD(tg_data_c,tg_list_c)
        dd<-dd/sd(dd)
          tg_list_selected<-c()
          ee_all<-c()
        ee_all<-c(ee_all,RMSE_row(dd))
        for(ii in 1:ROUNDS)
        {
        ccc<-cor(t(dd),t(module_RB))
                tg_id_c<-which(ccc==max(ccc))[1]
                if(ccc[tg_id_c]>0)
                {
                  tg_list_selected<-c(tg_list_selected,tg_id_c)
                ttt_ccc<-module_RB[tg_id_c,]%*%t(module_RB[tg_id_c,])/sum((module_RB[tg_id_c,])^2)
                dd<-dd-dd%*%ttt_ccc
                  ee_all<-c(ee_all,RMSE_row(dd))
                }
        }
        ee_all<-ee_all/ee_all[1]
          ee_all2<-ee_all[-1]
          if(length(ee_all2)>0)
          {
                for(i in 2:length(ee_all))
                {
                        ee_all2[i-1]<-ee_all[i-1]-ee_all[i]
                }
          }
          ee_all<-ee_all[-1]
          ccc<-c()
          if(length(ee_all)>0)
          {
             ccc<-rbind(tg_list_selected,ee_all,ee_all2)
 colnames(ccc)<-1:ncol(ccc)
          }  
         return(ccc)
}

compute_CompRowspace_NN<-function(data_cc=tg_data,module_RB=module_selected[[4]],ROUNDS=3)
{
        data_CORS_cancer<-data_cc
        for(ii in 1:ROUNDS)
        {
        ccc<-cor(t(data_CORS_cancer),t(module_RB))
        if(nrow(module_RB)>1)
        {
                ddd<-apply(ccc,1,order)
                eee<-ddd[nrow(ddd),]
                eee<-eee*(apply(ccc,1,max)>0)
                tg_ids<-sort(setdiff(unique(eee),0))
                if(length(tg_ids)>0)
                {
                        for(i in 1:length(tg_ids))
                        {
                                ttt_ccc<-module_RB[tg_ids[i],]%*%t(module_RB[tg_ids[i],])/sum((module_RB[tg_ids[i],])^2)
                                data_CORS_cancer[names(which(eee==tg_ids[i])),]<-data_CORS_cancer[names(which(eee==tg_ids[i])),]-data_CORS_cancer[names(which(eee==tg_ids[i])),]%*%ttt_ccc
                        }
                }
        }
        else
        {
                ddd<-apply(ccc,1,order)
                eee<-ddd
                eee<-eee*(as.vector(ccc)>0)
                tg_ids<-sort(setdiff(unique(eee),0))
                if(length(tg_ids)>0)
                {
                        for(i in 1:length(tg_ids))
                        {
                                ttt_ccc<-module_RB[tg_ids[i],]%*%t(module_RB[tg_ids[i],])/sum((module_RB[tg_ids[i],])^2)
                                data_CORS_cancer[names(which(eee==tg_ids[i])),]<-data_CORS_cancer[names(which(eee==tg_ids[i])),]-data_CORS_cancer[names(which(eee==tg_ids[i])),]%*%ttt_ccc
                        }
                }
        }
        }
        return(data_CORS_cancer)
}



compute_IM_stat<-function(tg_list_c,immune_cell_uni_table=immune_cell_uni_table0_GS,IM_id_list0=ICTD::IM_id_list)
{
	ccc<-c()
	nn<-c()
	for(i in 1:length(tg_list_c))
	{       
        ccc0<-c()
        for(j in 1:length(ICTD::IM_id_list))
        {       
      if(length(ICTD::IM_id_list[[j]])>1)
      {
            cc0<-apply(immune_cell_uni_table[tg_list_c[[i]],IM_id_list0[[j]]],1,sum)/sum((1/(1:length(IM_id_list0[[j]]))))
      }
      else
      {
           cc0<-immune_cell_uni_table[tg_list_c[[i]],IM_id_list0[[j]]]
      }
                ccc0<-cbind(ccc0,cc0)
        }
        colnames(ccc0)<-names(ICTD::IM_id_list)
      ddd<-apply(ccc0,2,mean)
      ccc<-rbind(ccc,ddd)
      nn<-c(nn,colnames(ccc0)[which(ddd==max(ddd))[1]])
	}
	rownames(ccc)<-nn
	return(ccc)
}


identify_max_base<-function(tg_R1_c,data.matrix,tProp)
{
	cc<-Compute_Rbase_SVD(data.matrix,tg_R1_c)
	tg_R1_c_list<-list()
	dd<-cor(t(cc),t(tProp))
	for(i in 1:ncol(dd))
	{
		tg_id<-which(dd[,i]==max(dd[,i]))[1]
		tg_R1_c_list[[i]]<-tg_R1_c[[tg_id]]
	}
	names(tg_R1_c_list)<-colnames(dd)
	return(tg_R1_c_list)
}

identify_max_base2<-function(tg_R1_c,data.matrix,tProp)
{
	cc<-Compute_Rbase_SVD(data.matrix,tg_R1_c)
	tg_R1_c_list<-list()
	dd<-cor(t(cc),t(tProp))
	ff<-c()
	for(i in 1:ncol(dd))
	{
		tg_id<-which(dd[,i]==max(dd[,i]))[1]
		tg_R1_c_list[[i]]<-tg_R1_c[[tg_id]]
		ff<-c(ff,tg_id)
	}
	print(ff)
	names(tg_R1_c_list)<-colnames(dd)
	return(tg_R1_c_list)
}

select_R_base<-function(tg_R1_c,tg_id)
{
	tg_R1_cc<-list()
	for(i in 1:length(tg_id))
	{
		tg_R1_cc[[i]]<-tg_R1_c[[tg_id[i]]]
	}
	names(tg_R1_cc)<-names(tg_R1_c)[tg_id]
	return(tg_R1_cc)
}