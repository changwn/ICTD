
NMF_method1<-function(tg_list,data_ng,data_normalized,max_ES_cut=0.3)
{
	tg_selected_R4_RR<-tg_list
	data.matrix0_s<-data_ng
	data_23_s<-data_normalized
	Rbase_selected_R4_RR<-Compute_Rbase_SVD(data.matrix0,tg_selected_R4_RR)
	stat_selected_R4_RR<-compute_IM_stat(tg_list_c=tg_selected_R4_RR)
	stat_selected_R4_RR_max<-apply(stat_selected_R4_RR,1,max)
	tg_I_id<-which(stat_selected_R4_RR_max==max(stat_selected_R4_RR_max))[1]
	tg_selected_id<-tg_I_id
	X_I<-data.matrix0[tg_selected_R4_RR[[tg_I_id]],]
	N<-0
	NMF_list_c<-list()
	N<-N+1
	NMF_list_c[[N]]<-tg_selected_R4_RR[[tg_I_id]]
	names(NMF_list_c)<-"T"
	NMF_self_c1<-Building_NMF_input_no_cancer(tg_data=X_I,tg_list=NMF_list_c,tg_list_add=list())
	V_I<-t(as.matrix(Rbase_selected_R4_RR[tg_I_id,]))
	V_I_n<-normalize_data2(V_I)
	V_c<-V_I_n
	ccc<-compute_CompRowspace_NN(data_23_s,module_RB=V_c,ROUNDS=min(nrow(V_I),3))
	ES_I<-RMSE_row(ccc)/RMSE_row(data_23_s)
	ES_base_c<-Compute_base_ES_score(ES_I,tg_selected_R4_RR)
	ES_max<-max(ES_base_c)
	while(ES_max>max_ES_cut)
	{
		N<-N+1
		tg_id<-which(abs(ES_base_c-max(ES_base_c))<0.01)
		tg_id<-setdiff(tg_id,tg_selected_id)
		tg_id_add<-tg_id[which(stat_selected_R4_RR_max[tg_id]==max(stat_selected_R4_RR_max[tg_id]))[1]]
		tg_selected_id<-c(tg_selected_id,tg_id_add)
		NMF_list_c[[N]]<-tg_selected_R4_RR[[tg_id_add]]
		names(NMF_list_c)<-names(tg_selected_R4_RR)[tg_selected_id]
		NMF_self_c1<-Building_NMF_input_no_cancer(tg_data=data.matrix0_s,tg_list=NMF_list_c,tg_list_add=list())
		qnmf_result_c <- run_NMF(NMF_self_c1,RR0=0.8, maxIter=2000)
		X1 = qnmf_result_c[["X1"]]
		U = qnmf_result_c[["U"]]
		V = qnmf_result_c[["V"]]
		ttt1 = qnmf_result_c[["ttt1"]]
		########check correlation between: V and true Proportion
		ictd.ccc= t(cor(t(tProp),V) ) 
		o6_predict_ture = apply(ictd.ccc, 1, max)
		o6_predict_ture2 = apply(ictd.ccc, 2, max)
		print(N)
		print(o6_predict_ture2)
		print(apply(cor(t(tProp),t(Rbase_selected_R4_RR[tg_selected_id,])),1,max))
		V_c<-normalize_data2(t(V))
		ccc<-compute_CompRowspace_NN(data_23_s,module_RB=V_c,ROUNDS=min(nrow(V_I),3))
		ES_I<-RMSE_row(ccc)/RMSE_row(data_23_s)
		ES_base_c<-Compute_base_ES_score(ES_I,tg_selected_R4_RR)
		ES_max<-max(ES_base_c)
	}
	return(qnmf_result_c)
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

#######################
run_NMF<-function(NMF_input1=NMF_input1,RR0=0.8, maxIter=20000, tProp){

###########Preprocess the data
P_fix=NMF_input1$S_indi
NMF_indi_all = NMF_input1$NMF_indi_all

data_t=NMF_input1$NMF_data
addP=NMF_input1$NMF_P_pre
NMF_indi_all0<-NMF_indi_all[which(apply(NMF_indi_all,1,sum)<=2),]
NMF_indi_all01<-NMF_indi_all0[,which(apply(NMF_indi_all0,2,sum)>0 | P_fix==0)] #?????0#######if 5 is okay, originally it was 10
NMF_indi_all01<-NMF_indi_all01[which(apply(NMF_indi_all01,1,sum)>0),]
NMF_indi_all=NMF_indi_all01[rownames(NMF_indi_all01)%in%rownames(data_t),]
X1=data_t[match(rownames(NMF_indi_all),rownames(data_t)),]###take only those rows that correspond to rows in NMF_indi_all
K=ncol(NMF_indi_all)
indiS=1-NMF_indi_all
indiS[,which(P_fix==0)]=0*indiS[,which(P_fix==0)]
###########
        
###########Parameter settings
RR=RR0##penalty parameter for constraints on NMF_indi_all
indiS_method="nonprdescent" ##the updating scheme for the structural constraints
iter=maxIter
alpha=beta=gamma=roh=0
#gamma=0.1
#roh=0.1 ##roh has to be equal to 0 in order to implement addP
UM=VM=NULL
qq=1
epslog=8
nPerm=2
initial_U=initial_V=NULL
mscale=0###########1
cnormalize=1
###########

###########Run the constrained qNMF
ttt1=qnmf_indisS_all_revise_addP_RR(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,RR,qq,iter,epslog,mscale,cnormalize,addP, P_fix, tProp)
#ttt1=qnmf_indisS_all_revise_addP_old(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale,cnormalize,addP, P_fix)
###########
X1=ttt1$X1
U=ttt1$U
V=ttt1$V
extraOut=ttt1$extraOut
#sum((X1-(ttt1$U)%*%t(ttt1$V))^2)
#diag(cor(ttt1$U,NMF_indi_all, method="spearman"))
###########Run cross validation of the constrained qNMF
#devs1=cv_qnmf(X1,indiS_method,NMF_indi_all,alpha,beta,gamma,roh,theta,qq,iter,epslog,nPerm,mscale,cnormalize)
#############################################

return(list(X1=X1, U=U, V=V, extraOut=extraOut, ttt1=ttt1))


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
