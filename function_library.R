

print("Loading Libraries")
library(ade4)
library(adegraphics)
library(compositions)
library(dendextend)
library(dplyr)
library(dunn.test)
library(effsize)
library(ggplot2)
library(ggrepel)
library(gplots)
library(igraph)
library(MASS)
library(metafor)
#library(metap)
library(pcaPP)
library(pROC)
library(psych)
library(randomForestExplainer)
library(RColorBrewer)
library(robumeta)
library(sfsmisc)
library(vegan)
library(ccrepe)


## Stage 1##
cat("Stage 1")
cat("Defining Functions")

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

compare_features <- function(data,Group)
{
	comparison_matrix <- as.data.frame(matrix(NA,ncol(data),2))
	rownames(comparison_matrix) <- colnames(data)
	colnames(comparison_matrix) <- c("estimate","pvalue")
	feature_list <- colnames(data)
	for(i in 1:length(feature_list))
	{
		feature_name <- feature_list[i]
		if(grepl("detect_",feature_name))
		{	
			#Do a Fishers' test
			mat <- matrix(NA,2,2)
			mat[1,1] <- length(which(data[Group,feature_name]==1)) + 1
			mat[1,2] <- length(which(data[Group,feature_name]==0)) + 1
			mat[2,1] <- length(which(data[!(Group),feature_name]==1)) + 1
			mat[2,2] <- length(which(data[!(Group),feature_name]==0)) + 1
			temp_fisher <- fisher.test(mat)
			comparison_matrix[feature_name,1] <- log((as.numeric(temp_fisher$estimate)+0.00001),10)
			print(mat)
			print(temp_fisher$estimate)
			print(log((as.numeric(temp_fisher$estimate)+0.00001),10))
			comparison_matrix[feature_name,2] <- temp_fisher$p.value
		}
		else
		{
			comparison_matrix[feature_name,1] <- effsize::cohen.d(data[Group,feature_name],data[!(Group),feature_name])$estimate
			comparison_matrix[feature_name,2] <- as.numeric(wilcox.test(data[,feature_name]~Group)$p.value)
		}
	}
	return(comparison_matrix)
}

lm_mat1 <- function(data)
{
	est_matrix <- matrix(0,ncol(data),ncol(data))
	rownames(est_matrix) <- colnames(data)
	colnames(est_matrix) <- rownames(data)
	pval_matrix <- matrix(1,ncol(data),ncol(data))
	rownames(pval_matrix) <- colnames(data)
	colnames(pval_matrix) <- rownames(data)
	for(i in 1:ncol(data))
	{
		feature1 <- colnames(data)[i]
		for(j in 1:ncol(data))
		{
			feature2 <- colnames(data)[j]
			taxa1 <- sub("detect_","",sub("abund_","",feature1))
			taxa2 <- sub("detect_","",sub("abund_","",feature2))
			if(taxa1 != taxa2)
			{
				if(grepl("detect_",feature1))
				{
					temp_lm <- summary(glm(as.formula(paste0(feature1,"~",feature2)),data=data,family="binomial"))
					est_matrix[i,j] <- temp_lm$coefficients[2,3]
					pval_matrix[i,j] <- temp_lm$coefficients[2,4]
				}
				else
				{
					temp_lm <- summary(glm(as.formula(paste0(feature1,"~",feature2)),data=data))
					est_matrix[i,j] <- temp_lm$coefficients[2,3]
					pval_matrix[i,j] <- temp_lm$coefficients[2,4]
				}
			}
		}
	}
	return_list <- list("est"=est_matrix,"p.value"=pval_matrix)
	return(return_list)
}

compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
	colnames(temp_meta) <- c("dataset","ri","ni")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		dat1 <- data[data[,grouping_variable]==group,var1]
		dat2 <- data[data[,grouping_variable]==group,var2]
		temp_meta[i,2] <- cor.fk(dat1,dat2)
		temp_meta[i,3] <- length(dat1)
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	print(temp_meta)
	temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
	temp_meta <- temp_meta[!is.na(temp_meta$ri),]
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
	detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
	rownames(detection_matrix) <- var1_list
	colnames(detection_matrix) <- grouping_list
	for(i in 1:length(var1_list))
	{
		var1 <- var1_list[i]
		for(j in 1:length(grouping_list))
		{
			group <- grouping_list[j]
			detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]>0))/length(data[data[,grouping_variable]==group,var1])
		}
	}
	return(detection_matrix)
}

compute_hedges_g <- function(data,var1_list,var2,grouping_variable,grouping_list)
{
	hedges_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
	colnames(hedges_matrix) <- grouping_list
	rownames(hedges_matrix) <- var1_list
	p_value_matrix <- data.frame(matrix(1,length(var1_list),length(grouping_list)))
	colnames(p_value_matrix) <- grouping_list
	rownames(p_value_matrix) <- var1_list
	for(i in 1:length(var1_list))
	{
		var1 <- var1_list[i]
		for(j in 1:length(grouping_list))
		{
			group <- grouping_list[j]
			data_group_Case <- data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1]
			data_group_Control <- data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1]
			hedges_matrix[i,j] <- as.numeric(effsize::cohen.d(data_group_Case,data_group_Control,hedges.correction=TRUE)$estimate)
			p_value_matrix[i,j] <- as.numeric(wilcox.test(data_group_Case,data_group_Control)$p.value)
		}
	}
	return_list = list("hedges"=hedges_matrix,"p_value"=p_value_matrix)
	return(return_list)
}

compute_meta_effsize_net <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),7))
	colnames(temp_meta) <- c("dataset","m1i","m2i","sd1i","sd2i","n1i","n2i")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(paste0(i,",",group))
		data_group_Case <- data[(data[,var2]!="G1")&(data[,grouping_variable]==group),var1]
		data_group_Control <- data[(data[,var2]=="G1")&(data[,grouping_variable]==group),var1]
		temp_meta[i,2] <- mean(data_group_Case)
		temp_meta[i,3] <- mean(data_group_Control)
		temp_meta[i,4] <- sd(data_group_Case)
		temp_meta[i,5] <- sd(data_group_Control)
		temp_meta[i,6] <- length(data_group_Case)
		temp_meta[i,7] <- length(data_group_Control)
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		print(temp_meta)
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="SMD",m1i=m1i,m2i=m2i,sd1i=sd1i,sd2i=sd2i,n1i=n1i,n2i=n2i,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	#res$slabs <- rownames(temp_meta)
	return(res)
}

compute_meta_lm <- function(data,var1,var2,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		vec1 <- data[data[,grouping_variable]==group,var1]
		vec2 <- data[data[,grouping_variable]==group,var2]
		#print(paste0(group,",",length(vec1[vec1>0])))
		if((length(vec1[vec1>0]) > 0)&&(length(vec2[vec2>0]) > 0))
		{
			#print(data[data[,grouping_variable]==group,c(var1,var2)])
			f <- as.formula(paste0(var1,"~",var2))
			temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
			summary_temp_rlm <- summary(temp_rlm)
			temp_meta[i,2] <- summary_temp_rlm$coefficients[2,3]
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
			temp_meta[i,6] <- sign(temp_meta[i,2])
			#levels <- unique(data[data[,grouping_variable]==group,var2])
			#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
			#temp <- ifelse(is.nan(temp),0,temp)
			#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
			#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
			#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
			#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
			#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		}
		else
		{
			temp_meta[i,2] <- 0
			temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
			temp_meta[i,4] <- 1
			temp_meta[i,5] <- 1
			temp_meta[i,6] <- 1
		}
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	print(temp_meta)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

compute_meta_lm_single_adjust <- function(data,var1,var2,var3,grouping_variable,grouping_list)
{
	temp_meta <- data.frame(matrix(0,length(grouping_list),6))
	colnames(temp_meta) <- c("dataset","ti","ni","mi","pi","di")
	for(i in 1:length(grouping_list))
	{
		group <- grouping_list[i]
		temp_meta[i,1] <- group
		print(group)
		f <- as.formula(paste0(var1,"~",var3,"+",var2))
		temp_rlm <- rlm(f,data=data[data[,grouping_variable]==group,])
		summary_temp_rlm <- summary(temp_rlm)
		temp_meta[i,2] <- summary_temp_rlm$coefficients[3,3]
		temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
		temp_meta[i,4] <- 1
		temp_meta[i,5] <- f.robftest(temp_rlm,var=var2)$p.value
		temp_meta[i,6] <- sign(temp_meta[i,2])
		#levels <- unique(data[data[,grouping_variable]==group,var2])
		#temp <- effsize::cohen.d(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1],data[(data[,var2]=="Control")&(data[,grouping_variable]==group),var1])$estimate
		#temp <- ifelse(is.nan(temp),0,temp)
		#temp <- ifelse(abs(temp)>1,0.99*sign(temp),temp)
		#print(data[(data[,var2]!="Control")&(data[,grouping_variable]==group),var1])
		#temp_meta[i,2] <- ifelse(is.nan(temp),0,temp)
		#temp_meta[i,2] <- #cor(data[data[,grouping_variable]==group,var1],data[data[,grouping_variable]==group,var2])
		#temp_meta[i,3] <- nrow(data[data[,grouping_variable]==group,])
	}
	temp_meta <- mutate(temp_meta,study_id=grouping_list)
	rownames(temp_meta) <- grouping_list
	#temp_meta <- temp_meta %>% select(study_id, ri:ni)
	temp_meta <- escalc(measure="ZPCOR",mi=mi,ni=ni,ti=ti,data=temp_meta)
	res <- rma(yi, vi, data=temp_meta)
	res$ids <- rownames(temp_meta)
	res$slabs <- rownames(temp_meta)
	return_list <- list("df_studies"=temp_meta,"model"=res)
	return(return_list)
}

compute_meta_lm_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
	return_out <- as.data.frame(matrix(NA,length(feature_list),9))
	rownames(return_out) <- feature_list
	colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir")
	for(i in 1:length(feature_list))
	{
		species_name <- feature_list[i]
		temp_res <- compute_meta_lm(data,species_name,metadata_var,grouping_var,grouping_list)
		return_out[i,"beta"] <- temp_res$model$beta
		return_out[i,"pval"] <- temp_res$model$pval
		return_out[i,"ci.ub"] <- temp_res$model$ci.ub
		return_out[i,"ci.lb"] <- temp_res$model$ci.lb
		return_out[i,"tau2"] <- temp_res$model$tau2
		return_out[i,"QE"] <- temp_res$model$QE
		return_out[i,"QEp"] <- temp_res$model$QEp
	}
	return_out$qval <- p.adjust(return_out$pval,method="fdr")
	return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.05,2*sign(return_out$beta),sign(return_out$beta)))
	return(return_out)
}

batch_rlm_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	df_est <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_est) <- grouping_list
	colnames(df_est) <- variable_group
	df_p_val <- as.data.frame(matrix(1,length(grouping_list),length(variable_group)))
	rownames(df_p_val) <- grouping_list
	colnames(df_p_val) <- variable_group
	for(i in 1:length(grouping_list))
	{
		study_name <- grouping_list[i]
		#print(study_name)
		study_samples <- rownames(metadata[metadata$study_name == study_name,])
		for(j in 1:length(variable_group))
		{
			species_name <- variable_group[j]
			species_val <- data[study_samples,species_name]
			#print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				temp_rlm <- rlm(data[study_samples,species_name]~metadata[study_samples,metadata_feature])
				summary_temp_rlm <- summary(temp_rlm)
				df_est[i,j] <- summary_temp_rlm$coefficients[2,3]
				df_p_val[i,j] <- f.robftest(temp_rlm)$p.value
			}
		}
		
	}
	df_q_val <- apply(df_p_val,2,p.adjust)
	l_fisher <- p.adjust(apply(df_q_val,2,function(x)(sumlog(x)$p)),method="fdr")
	print("l_fisher generated")
	df_dir <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
	rownames(df_dir) <- grouping_list
	colnames(df_dir) <- variable_group
	for(i in 1:length(grouping_list))
	{
		for(j in 1:length(variable_group))
		{
			df_dir[i,j] <- ifelse(df_q_val[i,j]<=0.10,3*sign(df_est[i,j]),ifelse(df_p_val[i,j]<=0.05,2*sign(df_est[i,j]),1*sign(df_est[i,j])))
		}
	}
	return_list <- list("est"=df_est,"p.value"=df_p_val,"q.value"=df_q_val,"fisher"=l_fisher,"dir"=df_dir)
	return(return_list)
}

batch_rem_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
	common_rows <- intersect(rownames(data),rownames(metadata))
	df_est <- as.data.frame(matrix(0,length(variable_group),5))
	rownames(df_est) <- variable_group
	colnames(df_est) <- c("est","pval","qval","dir",metadata_feature)
	for(j in 1:length(variable_group))
	{
			species_name <- variable_group[j]
			print(species_name)
			species_val <- data[common_rows,species_name]
			metadata_val <- metadata[common_rows,metadata_feature]
			metadata_grouping <- metadata[common_rows,grouping_feature]
			print(length(species_val[species_val>0]))
			if(length(species_val[species_val>0])>0)
			{
				df_temp <- data.frame(sp=species_val,meta_val=metadata_val,meta_grp=metadata_grouping,row.names=common_rows)
				temp_res <- compute_meta_lm(df_temp,"sp","meta_val","meta_grp",grouping_list)
				print(temp_res$model)
				df_est[j,1] <- as.numeric(temp_res$model$beta)
				df_est[j,2] <- temp_res$model$pval
			}
		
	}
	df_est[,3] <- p.adjust(df_est[,2],method="fdr")
	df_est[,4] <- ifelse(df_est[,3]<=0.10,3*sign(df_est[,1]),ifelse(df_est[,2]<=0.05,2*sign(df_est[,1]),1*sign(df_est[,1])))
	df_est[,5] <- metadata_feature
	return(df_est)
}

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

range_scale=function(x)
{
	y <- (x-min(x))/(max(x)-min(x));
	return(y);
}

rem_network1 <- function(data,species_list,group_name,study_list,q_threshold,c_threshold)
{
        species_data <- data[,species_list]
        species_data$group <- data[,group_name]
        est_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
		rownames(est_matrix) <- species_list
		colnames(est_matrix) <- species_list
		pval_matrix <- as.data.frame(matrix(1,length(species_list),length(species_list)))
		rownames(pval_matrix) <- species_list
		colnames(pval_matrix) <- species_list
		consistency_matrix <- as.data.frame(matrix(1,length(species_list),length(species_list)))
		rownames(consistency_matrix) <- species_list
		colnames(consistency_matrix) <- species_list
        for(i in 1:length(species_list))
        {
                species1 <- species_list[i]
                #print(species1)
                for(j in 1:length(species_list))
                {
                        species2 <- species_list[j]
                        print(paste0(species1,",",species2))
                        if(species1 != species2)
                        {
								tryCatch(               
								expr = {                     
										temp_rem <- compute_meta_corr(data,species1,species2,group_name,study_list)
										#print(temp_rem)
										est_matrix[species1,species2] <- as.numeric(temp_rem$model$beta)
										pval_matrix[species1,species2] <- temp_rem$model$pval
										consistency_matrix[species1,species2] <- length(which(sign(temp_rem$df_studies$ri)==sign(as.numeric(temp_rem$model$beta))))/ nrow(temp_rem$df_studies)
									},
									error = function(e){ 
										print(e)
										print("Error observed. Moving to next")
									},
									finally = {            
										print("finally Executed")
									}
								)
                               
                        }
                }
        }
        qval_matrix <- apply(pval_matrix,2,function(x)(p.adjust(x,method="bonferroni")))
		dir_matrix <- as.data.frame(matrix(0,length(species_list),length(species_list)))
		rownames(dir_matrix) <- species_list
		colnames(dir_matrix) <- species_list
		for(i in 1:length(species_list))
        {
          for(j in 1:length(species_list))
          {
			dir_matrix[i,j] <- ifelse((qval_matrix[i,j]<=q_threshold)&(consistency_matrix[i,j]>=c_threshold),sign(est_matrix[i,j]),0)
		  }
		}
		return_list <- list("est"=est_matrix,"pval"=pval_matrix,"qval"=qval_matrix,"dir"=dir_matrix)
		return(return_list)
}

rem_network2 <- function(data,species_list,feature_list,group_name,study_list)
{
        species_data <- data[,species_list]
        species_data$group <- data[,group_name]
        est_matrix <- as.data.frame(matrix(0,length(species_list),length(feature_list)))
		rownames(est_matrix) <- species_list
		colnames(est_matrix) <- feature_list
		pval_matrix <- as.data.frame(matrix(1,length(species_list),length(feature_list)))
		rownames(pval_matrix) <- species_list
		colnames(pval_matrix) <- feature_list
        for(i in 1:length(species_list))
        {
                species1 <- species_list[i]
                #print(species1)
                for(j in 1:length(feature_list))
                {
                        species2 <- feature_list[j]
                        print(paste0(species1,",",species2))
                        if(species1 != species2)
                        {
								tryCatch(               
								expr = {                     
										temp_rem <- compute_meta_lm(data,species1,species2,group_name,study_list)
										#print(temp_rem)
										est_matrix[species1,species2] <- as.numeric(temp_rem$model$beta)
										pval_matrix[species1,species2] <- temp_rem$model$pval
		
									},
									error = function(e){         
										print("Error observed. Moving to next")
										print("e")
									},
									finally = {            
										print("finally Executed")
									}
								)
                               
                        }
                }
        }
        qval_matrix <- apply(pval_matrix,2,function(x)(p.adjust(x,method="fdr")))
		dir_matrix <- as.data.frame(matrix(0,length(species_list),length(feature_list)))
		rownames(dir_matrix) <- species_list
		colnames(dir_matrix) <- feature_list
		for(i in 1:length(species_list))
        {
          for(j in 1:length(feature_list))
          {
			dir_matrix[i,j] <- ifelse(qval_matrix[i,j]<=0.15,2*sign(est_matrix[i,j]),ifelse(pval_matrix[i,j]<=0.05,sign(est_matrix[i,j]),0))
		  }
		}
		return_list <- list("est"=est_matrix,"pval"=pval_matrix,"qval"=qval_matrix,"dir"=dir_matrix)
		return(return_list)
}

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
        detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
        rownames(detection_matrix) <- var1_list
        colnames(detection_matrix) <- grouping_list
        for(i in 1:length(var1_list))
        {
                var1 <- var1_list[i]
                for(j in 1:length(grouping_list))
                {
                        group <- grouping_list[j]
                        detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]!=0))/length(data[data[,grouping_variable]==group,var1])
                }
        }
        return(detection_matrix)
}

wilcox_batch = function(x,y)
{
	p_array <- NULL;
	type_array <- NULL;
	mean1_array <- NULL;
	mean2_array <- NULL;
	x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
	y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
	z <- intersect(rownames(x),rownames(y));
	for(i in 1:length(z))
	{
		p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
		type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
		mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
		mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
		i <- i + 1;
	}
	out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array,method="fdr"),mean1_array,mean2_array));
	rownames(out) <- z;
	out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
	return(t(out));
}

library(randomForest)

iterative_rf = function(data,window,disease,control,iter)
{
	set.seed(100);
	featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
	AUCArray <- NULL;
	SensitivityArray <- NULL;
	SpecificityArray <- NULL;
	#threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
	for(i in 1:iter)
	{
		trainDisease <- sample(intersect(window,disease),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		if(as.integer(length(intersect(window,disease))/2) < length(intersect(window,control)))
		{
			trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=FALSE);
		}
		else
		{
			trainControl <- sample(intersect(window,control),as.integer(length(intersect(window,disease))/2),replace=TRUE);
		}
		tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
		TrainDiseaseTags <- NULL;
		TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
		TrainControlTags <- NULL;
		TrainControlTags[1:length(trainDisease)] <- "Control";
		TrainTags <- c(TrainDiseaseTags,TrainControlTags);
		rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
		featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
		testDisease <- setdiff(intersect(window,disease),trainDisease);
		testControl <- setdiff(intersect(window,control),trainControl);
		tempTest <- rbind(data[testDisease,],data[testControl,]);
		TestDiseaseTags <- NULL;
		TestDiseaseTags[1:length(testDisease)] <- "Diseased";
		TestControlTags <- NULL;
		TestControlTags[1:length(testControl)] <- "Control";
		TestTags <- c(TestDiseaseTags,TestControlTags);
		rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE)
		AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
		SensitivityArray[i] <- length(which(predict(rfTempComp,data[testDisease,])=="Diseased"))/length(predict(rfTempComp,data[testDisease,]));
		SpecificityArray[i] <- length(which(predict(rfTempComp,data[testControl,])=="Control"))/length(predict(rfTempComp,data[testControl,]));
		print(i);
		i <- i + 1;
	}
	colnames(featureProfile) <- colnames(tempTrain);
	returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray);
	return(returnList);
		
}

batch_dunns = function(data,factor)
{
	species_list <- colnames(data);
	length_factor <- length(unique(factor));
	pair_factor <- length_factor*(length_factor-1)/2;
	comparisons <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	comparisons_adjusted <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	z_changes <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	temp_comp <- NULL;
	k_test <- NULL;
	for(i in 1:length(species_list))
	{
		print(i)
		temp_dunn <- dunn.test(data[,species_list[i]],factor,method="bh");
		comparisons[i,] <- temp_dunn$P;
		comparisons_adjusted[i,] <- temp_dunn$P.adjusted;
		z_changes[i,] <- temp_dunn$Z;	
		k_test[i] <- kruskal.test(data[,species_list[i]]~factor)$p.value
		temp_comp <- temp_dunn$comparisons;
		i <- i + 1;
	}
	names(k_test) <- species_list
	rownames(comparisons) <- species_list;
	colnames(comparisons) <- temp_comp;
	rownames(z_changes) <- species_list;
	colnames(z_changes) <- temp_comp;
	#comparisons_adjusted <- apply(comparisons,2,function(x)(p.adjust(x,method="fdr")));
	final_trends <- as.data.frame(matrix(NA,length(species_list),pair_factor));
	for(i in 1:length(species_list))
	{
		for(j in 1:length(temp_comp))
		{
			final_trends[i,j] <- ifelse(comparisons_adjusted[i,j] <= 0.10,2*sign(z_changes[i,j]),sign(z_changes[i,j]));
			j <- j + 1;
		}
		i <- i + 1;
	}
	rownames(final_trends) <- species_list;
	colnames(final_trends) <- temp_comp;
	return_list <- list("kruskal"=k_test,"PValue"=comparisons,"CorrectedP"=comparisons_adjusted,"Z"=z_changes,"final_trends"=final_trends);
	return(return_list);
}
