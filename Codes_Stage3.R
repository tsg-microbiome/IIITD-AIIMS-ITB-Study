print("Performing Marker Identification using Random Forest and Principal Coordinate Analysis")

common_rows <- intersect(rownames(df_clr_Microbiome_Species_AIIMS_ITB),rownames(df_clr_Fungome_Species_AIIMS_ITB))

AIIMS_ITB_All_Taxa_Abundance_Combined <- as.data.frame(cbind(df_clr_Microbiome_Genus_AIIMS_ITB[common_rows,],df_clr_Microbiome_Species_AIIMS_ITB[common_rows,],df_clr_Fungome_Genus_AIIMS_ITB[common_rows,],df_clr_Fungome_Species_AIIMS_ITB[common_rows,]))

AIIMS_ITB_All_Taxa_Detection_Combined <- as.data.frame(cbind(apply(AIIMS_ITB_All_Taxa_Abundance_Combined,2,function(x)(ifelse(x>0,1,0)))))

FungomeRows <- intersect(colnames(AIIMS_ITB_All_Taxa_Abundance_Combined),c(colnames(df_clr_Fungome_Genus_AIIMS_ITB),colnames(df_clr_Fungome_Species_AIIMS_ITB)))

MicrobiomeRows <- intersect(colnames(AIIMS_ITB_All_Taxa_Abundance_Combined),c(colnames(df_clr_Microbiome_Genus_AIIMS_ITB),colnames(df_clr_Microbiome_Species_AIIMS_ITB)))

colnames(AIIMS_ITB_All_Taxa_Detection_Combined) <- paste0("detect_",colnames(AIIMS_ITB_All_Taxa_Detection_Combined))

colnames(AIIMS_ITB_All_Taxa_Abundance_Combined) <- paste0("abund_",colnames(AIIMS_ITB_All_Taxa_Abundance_Combined))

AIIMS_ITB_All_Taxa_Properties_Combined <- as.data.frame(cbind(AIIMS_ITB_All_Taxa_Abundance_Combined,AIIMS_ITB_All_Taxa_Detection_Combined))

AIIMS_ITB_All_Taxa_Properties_Combined$Group <- NA
AIIMS_ITB_All_Taxa_Properties_Combined[intersect(common_rows,AIIMS_ITB_Controls_Fungome),"Group"] <- "Control"
AIIMS_ITB_All_Taxa_Properties_Combined[intersect(common_rows,AIIMS_ITB_CD_Cases_Fungome),"Group"] <- "CD"
AIIMS_ITB_All_Taxa_Properties_Combined[intersect(common_rows,AIIMS_ITB_Cases_Fungome),"Group"] <- "ITB"
#AIIMS_ITB_All_Taxa_Properties_Combined$Group <- as.factor(AIIMS_ITB_All_Taxa_Properties_Combined$Group)

CD_ITB_All_Taxa_Properties_Combined <- AIIMS_ITB_All_Taxa_Properties_Combined[AIIMS_ITB_All_Taxa_Properties_Combined$Group %in% c("CD","ITB"),]
CD_ITB_All_Taxa_Properties_Combined$Group <- ifelse(CD_ITB_All_Taxa_Properties_Combined$Group == "ITB",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(700);
	numb<-numb_features[i];
	print(numb)
	RF_CD_ITB_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,CD_ITB_All_Taxa_Properties_Combined);SelectedFeatures <- names(tail(RF_CD_ITB_Full$importance[order(RF_CD_ITB_Full$importance),],numb));print(numb);
	RF_CD_ITB_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,CD_ITB_All_Taxa_Properties_Combined[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_CD_ITB_Selected$y,as.numeric(RF_CD_ITB_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_CD_ITB_Selected$y,as.numeric(RF_CD_ITB_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}

CD_ITB_roc_vec <- roc_vec
set.seed(700);CD_ITB_Selected_Features <- names(tail(RF_CD_ITB_Full$importance[order(RF_CD_ITB_Full$importance),],numb_features[which(CD_ITB_roc_vec==max(CD_ITB_roc_vec))]));RF_CD_ITB_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,CD_ITB_All_Taxa_Properties_Combined[,c(CD_ITB_Selected_Features,"Group")],localImp = TRUE);
CD_ITB_Feature_Importance <- as.data.frame(measure_importance(RF_CD_ITB_Selected))
plot_multi_way_importance(CD_ITB_Feature_Importance, no_of_labels = 30,y_measure="times_a_root",x_measure="accuracy_decrease")

Controls_CD_All_Taxa_Properties_Combined <- AIIMS_ITB_All_Taxa_Properties_Combined[AIIMS_ITB_All_Taxa_Properties_Combined$Group %in% c("Control","CD"),]
Controls_CD_All_Taxa_Properties_Combined$Group <- ifelse(Controls_CD_All_Taxa_Properties_Combined$Group == "CD",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(700);
	numb<-numb_features[i];
	print(numb)
	RF_Controls_CD_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,Controls_CD_All_Taxa_Properties_Combined)
	SelectedFeatures <- SelectedFeatures <- names(tail(RF_Controls_CD_Full$importance[order(RF_Controls_CD_Full$importance),],numb));
	RF_Controls_CD_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Controls_CD_All_Taxa_Properties_Combined[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_Controls_CD_Selected$y,as.numeric(RF_Controls_CD_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_Controls_CD_Selected$y,as.numeric(RF_Controls_CD_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}
Controls_CD_roc_vec <- roc_vec
set.seed(700);Controls_CD_Selected_Features <- names(tail(RF_Controls_CD_Full$importance[order(RF_Controls_CD_Full$importance),],numb_features[which(Controls_CD_roc_vec==max(Controls_CD_roc_vec))][1]));RF_Controls_CD_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Controls_CD_All_Taxa_Properties_Combined[,c(Controls_CD_Selected_Features,"Group")],localImp = TRUE);
Controls_CD_Feature_Importance <- as.data.frame(measure_importance(RF_Controls_CD_Selected))
plot_multi_way_importance(Controls_CD_Feature_Importance, no_of_labels = 40,y_measure="times_a_root",x_measure="accuracy_decrease")

Controls_ITB_All_Taxa_Properties_Combined <- AIIMS_ITB_All_Taxa_Properties_Combined[AIIMS_ITB_All_Taxa_Properties_Combined$Group %in% c("Control","ITB"),]
Controls_ITB_All_Taxa_Properties_Combined$Group <- ifelse(Controls_ITB_All_Taxa_Properties_Combined$Group == "ITB",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(700);
	numb<-numb_features[i];
	print(numb)
	RF_Controls_ITB_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,Controls_ITB_All_Taxa_Properties_Combined)
	SelectedFeatures <- names(tail(RF_Controls_ITB_Full$importance[order(RF_Controls_ITB_Full$importance),],numb));
	RF_Controls_ITB_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Controls_ITB_All_Taxa_Properties_Combined[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_Controls_ITB_Selected$y,as.numeric(RF_Controls_ITB_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_Controls_ITB_Selected$y,as.numeric(RF_Controls_ITB_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}
Controls_ITB_roc_vec <- roc_vec
set.seed(700);Controls_ITB_Selected_Features <- names(tail(RF_Controls_ITB_Full$importance[order(RF_Controls_ITB_Full$importance),],numb_features[which(Controls_ITB_roc_vec==max(Controls_ITB_roc_vec))][1]));RF_Controls_ITB_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Controls_ITB_All_Taxa_Properties_Combined[,c(Controls_ITB_Selected_Features,"Group")],localImp = TRUE);
Controls_ITB_Feature_Importance <- as.data.frame(measure_importance(RF_Controls_ITB_Selected))
plot_multi_way_importance(Controls_ITB_Feature_Importance, no_of_labels = 10,y_measure="times_a_root",x_measure="accuracy_decrease")

ITB_Others_All_Taxa_Properties_Combined <- AIIMS_ITB_All_Taxa_Properties_Combined[AIIMS_ITB_All_Taxa_Properties_Combined$Group %in% c("Control","ITB","CD"),]
ITB_Others_All_Taxa_Properties_Combined$Group <- ifelse(ITB_Others_All_Taxa_Properties_Combined$Group == "ITB",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(700);
	numb<-numb_features[i];
	print(numb)
	RF_ITB_Others_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,ITB_Others_All_Taxa_Properties_Combined)
	SelectedFeatures <- names(tail(RF_ITB_Others_Full$importance[order(RF_ITB_Others_Full$importance),],numb));
	RF_ITB_Others_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,ITB_Others_All_Taxa_Properties_Combined[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_ITB_Others_Selected$y,as.numeric(RF_ITB_Others_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_ITB_Others_Selected$y,as.numeric(RF_ITB_Others_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}
ITB_Others_roc_vec <- roc_vec
set.seed(700);ITB_Others_Selected_Features <- names(tail(RF_ITB_Others_Full$importance[order(RF_ITB_Others_Full$importance),],numb_features[which(ITB_Others_roc_vec==max(ITB_Others_roc_vec))][1]));RF_ITB_Others_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,ITB_Others_All_Taxa_Properties_Combined[,c(ITB_Others_Selected_Features,"Group")],localImp = TRUE);
ITB_Others_Feature_Importance <- as.data.frame(measure_importance(RF_ITB_Others_Selected))
plot_multi_way_importance(ITB_Others_Feature_Importance, no_of_labels = 10,y_measure="times_a_root",x_measure="accuracy_decrease")

df_MostPredictiveFeatures <- data.frame(Controls_ITB = Controls_ITB_roc_vec,CD_ITB = CD_ITB_roc_vec,Controls_CD = Controls_CD_roc_vec,row.names=numb_features)

CD_ITB_Iterative_top <- iterative_rf(CD_ITB_All_Taxa_Properties_Combined[,CD_ITB_Selected_Features],rownames(CD_ITB_All_Taxa_Properties_Combined),rownames(CD_ITB_All_Taxa_Properties_Combined[CD_ITB_All_Taxa_Properties_Combined$Group,]),rownames(CD_ITB_All_Taxa_Properties_Combined[!CD_ITB_All_Taxa_Properties_Combined$Group,]),50)

CD_ITB_Iterative_non_top <- iterative_rf(CD_ITB_All_Taxa_Properties_Combined[,setdiff(colnames(CD_ITB_All_Taxa_Properties_Combined),CD_ITB_Selected_Features)],rownames(CD_ITB_All_Taxa_Properties_Combined),rownames(CD_ITB_All_Taxa_Properties_Combined[CD_ITB_All_Taxa_Properties_Combined$Group,]),rownames(CD_ITB_All_Taxa_Properties_Combined[!CD_ITB_All_Taxa_Properties_Combined$Group,]),50)

Controls_CD_Iterative_top <- iterative_rf(Controls_CD_All_Taxa_Properties_Combined[,Controls_CD_Selected_Features],rownames(Controls_CD_All_Taxa_Properties_Combined),rownames(Controls_CD_All_Taxa_Properties_Combined[Controls_CD_All_Taxa_Properties_Combined$Group,]),rownames(Controls_CD_All_Taxa_Properties_Combined[!Controls_CD_All_Taxa_Properties_Combined$Group,]),50)

Controls_CD_Iterative_non_top <- iterative_rf(Controls_CD_All_Taxa_Properties_Combined[,setdiff(colnames(Controls_CD_All_Taxa_Properties_Combined),Controls_CD_Selected_Features)],rownames(Controls_CD_All_Taxa_Properties_Combined),rownames(Controls_CD_All_Taxa_Properties_Combined[Controls_CD_All_Taxa_Properties_Combined$Group,]),rownames(Controls_CD_All_Taxa_Properties_Combined[!Controls_CD_All_Taxa_Properties_Combined$Group,]),50)

Controls_ITB_Iterative_top <- iterative_rf(Controls_ITB_All_Taxa_Properties_Combined[,Controls_ITB_Selected_Features],rownames(Controls_ITB_All_Taxa_Properties_Combined),rownames(Controls_ITB_All_Taxa_Properties_Combined[Controls_ITB_All_Taxa_Properties_Combined$Group,]),rownames(Controls_ITB_All_Taxa_Properties_Combined[!Controls_ITB_All_Taxa_Properties_Combined$Group,]),50)

Controls_ITB_Iterative_non_top <- iterative_rf(Controls_ITB_All_Taxa_Properties_Combined[,setdiff(colnames(Controls_ITB_All_Taxa_Properties_Combined),Controls_ITB_Selected_Features)],rownames(Controls_ITB_All_Taxa_Properties_Combined),rownames(Controls_ITB_All_Taxa_Properties_Combined[Controls_ITB_All_Taxa_Properties_Combined$Group,]),rownames(Controls_ITB_All_Taxa_Properties_Combined[!Controls_ITB_All_Taxa_Properties_Combined$Group,]),50)

boxplot(100*Controls_ITB_Iterative_top$AUC,100*Controls_ITB_Iterative_non_top$AUC,outline=FALSE,col=c("dodgerblue","gold1"),ylim=c(40,100),cex.axis=2)
boxplot(100*Controls_CD_Iterative_top$AUC,100*Controls_CD_Iterative_non_top$AUC,outline=FALSE,col=c("dodgerblue","gold1"),ylim=c(40,100),cex.axis=2)
boxplot(100*CD_ITB_Iterative_top$AUC,100*CD_ITB_Iterative_non_top$AUC,outline=FALSE,col=c("dodgerblue","gold1"),ylim=c(40,100),cex.axis=2)


AllTopMarkers <- unique(c(CD_ITB_Selected_Features,Controls_CD_Selected_Features,Controls_ITB_Selected_Features))
CompareFeatures_Controls_ITB <- compare_features(Controls_ITB_All_Taxa_Properties_Combined[,AllTopMarkers],Controls_ITB_All_Taxa_Properties_Combined$Group)
CompareFeatures_Controls_CD <- compare_features(Controls_CD_All_Taxa_Properties_Combined[,AllTopMarkers],Controls_CD_All_Taxa_Properties_Combined$Group)
CompareFeatures_CD_ITB <- compare_features(CD_ITB_All_Taxa_Properties_Combined[,AllTopMarkers],CD_ITB_All_Taxa_Properties_Combined$Group)

df_FeaturePatterns <- as.data.frame(matrix(0,length(unique(c(CD_ITB_Selected_Features,Controls_CD_Selected_Features,Controls_ITB_Selected_Features))),3))
rownames(df_FeaturePatterns) <- unique(c(CD_ITB_Selected_Features,Controls_CD_Selected_Features,Controls_ITB_Selected_Features))
colnames(df_FeaturePatterns) <- c("Controls_v_ITB","Controls_v_CD","CD_v_ITB")
for(i in 1:nrow(df_FeaturePatterns))
{
	feature_name <- rownames(df_FeaturePatterns)[i]
	df_FeaturePatterns[feature_name,1] <- ifelse(p.adjust(CompareFeatures_Controls_ITB[feature_name,2],method="fdr")<=0.20,2*sign(CompareFeatures_Controls_ITB[feature_name,1]),sign(CompareFeatures_Controls_ITB[feature_name,1]))
	df_FeaturePatterns[feature_name,2] <- ifelse(p.adjust(CompareFeatures_Controls_CD[feature_name,2],method="fdr")<=0.20,2*sign(CompareFeatures_Controls_CD[feature_name,1]),sign(CompareFeatures_Controls_CD[feature_name,1]))
	df_FeaturePatterns[feature_name,3] <- ifelse(p.adjust(CompareFeatures_CD_ITB[feature_name,2],method="fdr")<=0.20,2*sign(CompareFeatures_CD_ITB[feature_name,1]),sign(CompareFeatures_CD_ITB[feature_name,1]))
	
}

df_FeaturePresence <- as.data.frame(matrix(0,length(unique(c(CD_ITB_Selected_Features,Controls_CD_Selected_Features,Controls_ITB_Selected_Features))),3))
rownames(df_FeaturePresence) <- unique(c(CD_ITB_Selected_Features,Controls_CD_Selected_Features,Controls_ITB_Selected_Features))
colnames(df_FeaturePresence) <- c("Controls_v_ITB","Controls_v_CD","CD_v_ITB")
for(i in 1:nrow(df_FeaturePresence))
{
	feature_name <- rownames(df_FeaturePatterns)[i]
	df_FeaturePresence[feature_name,1] <- ifelse(feature_name %in% Controls_ITB_Selected_Features,1,0)
	df_FeaturePresence[feature_name,2] <- ifelse(feature_name %in% Controls_CD_Selected_Features,1,0)
	df_FeaturePresence[feature_name,3] <- ifelse(feature_name %in% CD_ITB_Selected_Features,1,0)
	
}

mat <- t(df_FeaturePatterns)
hmp_FeaturePatterns <- heatmap.2(mat,density="none",trace="none",Rowv=FALSE,col=c("lightblue4","lightblue1","lightpink1","lightpink4"),lwid=c(1,5),margins=c(15,10),srtCol=90,sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat),rowsep=0:nrow(mat))

mat <- t(df_FeaturePresence)
mat <- mat[rev(colnames(hmp_FeaturePatterns$carpet)),rownames(hmp_FeaturePatterns$carpet)]
heatmap.2(mat[rev(colnames(hmp_FeaturePatterns$carpet)),rownames(hmp_FeaturePatterns$carpet)],density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","lightgreen"),lwid=c(1,5),margins=c(15,10),srtCol=90,sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat),rowsep=0:nrow(mat))

## Microbiome Species
SelectSpecies_AIIMS_ITB <- names(which(colSums(apply(df_clr_Microbiome_Species_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))/nrow(df_clr_Microbiome_Species_AIIMS_ITB)>=0.10))


dfSelectPCo_Microbiome_species = data.frame(PCo1_Kendall=pco_AIIMS_ITB_Kendall_Microbiome_Species$li$A1,
                                            PCo1_Aitchison=pco_AIIMS_ITB_Aitchison_Microbiome_Species$li$A1,
                                            PCo1_Jaccard=pco_AIIMS_ITB_Jaccard_Microbiome_Species$li$A1,
                                            row.names=rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li))


corr_Microbiome_Species_PCo = corr.test(df_clr_Microbiome_Species_AIIMS_ITB[rownames(dfSelectPCo_Microbiome_species),],dfSelectPCo_Microbiome_species,use="pairwise.complete",method="kendall")
SigAssoc_Microbiome_Species <- intersect(names(which(apply(corr_Microbiome_Species_PCo$p,1,function(x)(length(x[x<=0.05])))>=2)),SelectSpecies_AIIMS_ITB)

hmpMicrobiomeSignature_species <- heatmap.2(t(corr_Microbiome_Species_PCo$r[SigAssoc_Microbiome_Species,]),density="none",trace="none")

#heatmap.2(t(corr_Microbiome_Species_PCo$r[SigAssoc_Microbiome_Species,]),density="none",trace="none",margins=c(15,10),lwid=c(0.1,5),cexRow=0.8,Rowv = FALSE, cellnote=apply(corr_Microbiome_Species_PCo$p[SigAssoc_Microbiome_Species,],1,function(x)(ifelse(x<=0.05,"*",""))), col=rev(brewer.pal(4,"PiYG")),notecol="black")

## Microbiome Genus
SelectGenus_AIIMS_ITB <- names(which(colSums(apply(df_clr_Microbiome_Genus_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))/nrow(df_clr_Microbiome_Genus_AIIMS_ITB)>=0.10))
dfSelectPCo_Microbiome_genus = data.frame(PCo1_Kendall=pco_AIIMS_ITB_Kendall_Microbiome_Genus$li$A1,
                                            PCo1_Aitchison=pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li$A1,
                                            PCo1_Jaccard=pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li$A1,
                                            row.names=rownames(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li))



corr_Microbiome_Genus_PCo = corr.test(df_clr_Microbiome_Genus_AIIMS_ITB[rownames(dfSelectPCo_Microbiome_genus),],dfSelectPCo_Microbiome_genus,use="pairwise.complete",method="kendall")
SigAssoc_Microbiome_Genus <- intersect(names(which(apply(corr_Microbiome_Genus_PCo$p,1,function(x)(length(x[x<=0.05])))>=2)),SelectGenus_AIIMS_ITB)

hmpMicrobiomeSignature_genus <- heatmap.2(t(corr_Microbiome_Genus_PCo$r[SigAssoc_Microbiome_Genus,]),density="none",trace="none")

#hmpMicrobiomeSignature_genus <- heatmap.2(t(corr_Microbiome_Genus_PCo$r[SigAssoc_Microbiome_Genus,]),density="none",trace="none", margins=c(15,10),lwid=c(0.1,5),cexRow=0.8,Rowv = FALSE, cellnote=apply(corr_Microbiome_Genus_PCo$p[SigAssoc_Microbiome_Genus,],1,function(x)(ifelse(x<=0.05,"*",""))), col=rev(brewer.pal(4,"PiYG")),notecol="black")

## Fungome Species
SelectFungomeSpecies_AIIMS_ITB <- names(which(colSums(apply(df_clr_Fungome_Species_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))/nrow(df_clr_Fungome_Species_AIIMS_ITB)>=0.10))

dfSelectPCo_Fungome_species = data.frame(PCo2_Kendall=pco_AIIMS_ITB_Kendall_Fungome_Species$li$A2,PCo2_Jaccard=pco_AIIMS_ITB_Jaccard_Fungome_Species$li$A2,PCo3_Jaccard=pco_AIIMS_ITB_Jaccard_Fungome_Species$li$A3,row.names=rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li))

corr_Fungome_Species_PCo = corr.test(df_clr_Fungome_Species_AIIMS_ITB[rownames(dfSelectPCo_Fungome_species),],dfSelectPCo_Fungome_species,use="pairwise.complete",method="kendall")
SigAssoc_Fungome_Species <- intersect(names(which(apply(corr_Fungome_Species_PCo$p,1,function(x)(length(x[x<=0.05])))>=1)),SelectFungomeSpecies_AIIMS_ITB)

hmpFungomeSignature_species <- heatmap.2(t(corr_Fungome_Species_PCo$r[SigAssoc_Fungome_Species,]),density="none",trace="none")

#heatmap.2(t(corr_Fungome_Species_PCo$r[SigAssoc_Fungome_Species,]),density="none",trace="none", margins=c(15,10),lwid=c(0.1,5),Rowv=FALSE,cexRow=0.8, cellnote=apply(corr_Fungome_Species_PCo$p[SigAssoc_Fungome_Species,],1,function(x)(ifelse(x<=0.05,"*",""))), col=rev(brewer.pal(4,"PiYG")),notecol="black")

## Fungome Genus
SelectFungomeGenus_AIIMS_ITB <- names(which(colSums(apply(df_clr_Fungome_Genus_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))/nrow(df_clr_Fungome_Genus_AIIMS_ITB)>=0.10))

dfSelectPCo_Fungome_genus = data.frame(PCo1_Kendall=pco_AIIMS_ITB_Kendall_Fungome_Genus$li$A1,PCo1_Jaccard=pco_AIIMS_ITB_Jaccard_Fungome_Genus$li$A1,row.names=rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li))

corr_Fungome_Genus_PCo = corr.test(df_clr_Fungome_Genus_AIIMS_ITB[rownames(dfSelectPCo_Fungome_genus),],dfSelectPCo_Fungome_genus,use="pairwise.complete",method="kendall")

SigAssoc_Fungome_Genus <- intersect(names(which(apply(corr_Fungome_Genus_PCo$p,1,function(x)(length(x[x<=0.05])))>=1)),SelectFungomeGenus_AIIMS_ITB)

hmpFungomeSignature_genus <- heatmap.2(t(corr_Fungome_Genus_PCo$r[SigAssoc_Fungome_Genus,]),density="none",trace="none")

#heatmap.2(t(corr_Fungome_Genus_PCo$r[SigAssoc_Fungome_Genus,]),density="none",trace="none", margins=c(15,10),lwid=c(0.1,5),Rowv=FALSE,cexRow=0.8, cellnote=apply(corr_Fungome_Genus_PCo$p[SigAssoc_Fungome_Genus,],1,function(x)(ifelse(x<=0.05,"*",""))), col=rev(brewer.pal(4,"PiYG")),notecol="black")

OnlyPCoAssociatedMarkers <- setdiff(paste0("abund_",unique(c(rownames(hmpMicrobiomeSignature_species$carpet),rownames(hmpMicrobiomeSignature_genus$carpet),rownames(hmpFungomeSignature_species$carpet),rownames(hmpFungomeSignature_genus$carpet)))),AllTopMarkers)

ComparePCoFeatures_Controls_ITB <- compare_features(Controls_ITB_All_Taxa_Properties_Combined[,OnlyPCoAssociatedMarkers],Controls_ITB_All_Taxa_Properties_Combined$Group)
ComparePCoFeatures_Controls_CD <- compare_features(Controls_CD_All_Taxa_Properties_Combined[,OnlyPCoAssociatedMarkers],Controls_CD_All_Taxa_Properties_Combined$Group)
ComparePCoFeatures_CD_ITB <- compare_features(CD_ITB_All_Taxa_Properties_Combined[,OnlyPCoAssociatedMarkers],CD_ITB_All_Taxa_Properties_Combined$Group)

df_PCoFeaturePatterns <- as.data.frame(matrix(0,length(OnlyPCoAssociatedMarkers),3))
rownames(df_PCoFeaturePatterns) <- OnlyPCoAssociatedMarkers
colnames(df_PCoFeaturePatterns) <- c("Controls_v_ITB","Controls_v_CD","CD_v_ITB")
for(i in 1:nrow(df_PCoFeaturePatterns))
{
	feature_name <- rownames(df_PCoFeaturePatterns)[i]
	df_PCoFeaturePatterns[feature_name,1] <- ifelse(p.adjust(ComparePCoFeatures_Controls_ITB[feature_name,2],method="fdr")<=0.20,2*sign(ComparePCoFeatures_Controls_ITB[feature_name,1]),sign(ComparePCoFeatures_Controls_ITB[feature_name,1]))
	df_PCoFeaturePatterns[feature_name,2] <- ifelse(p.adjust(ComparePCoFeatures_Controls_CD[feature_name,2],method="fdr")<=0.20,2*sign(ComparePCoFeatures_Controls_CD[feature_name,1]),sign(ComparePCoFeatures_Controls_CD[feature_name,1]))
	df_PCoFeaturePatterns[feature_name,3] <- ifelse(p.adjust(ComparePCoFeatures_CD_ITB[feature_name,2],method="fdr")<=0.20,2*sign(ComparePCoFeatures_CD_ITB[feature_name,1]),sign(ComparePCoFeatures_CD_ITB[feature_name,1]))
	
}
df_PCoFeaturePatterns <- df_PCoFeaturePatterns[!is.na(rowSums(df_PCoFeaturePatterns)),]
mat <- t(df_PCoFeaturePatterns)
hmp_PCoFeaturePatterns <- heatmap.2(mat,density="none",trace="none",Rowv=FALSE,col=c("lightblue4","lightblue1","lightpink1","lightpink4"),lwid=c(1,5),margins=c(15,10),srtCol=90,sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat),rowsep=0:nrow(mat))

common_rows <- intersect(intersect(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li),rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li)),intersect(rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li),rownames(AIIMS_ITB_All_Taxa_Properties_Combined)))

df_SelectPCo_All_Taxa_Combined = data.frame(PCo1_Kendall_Microbiome_Species = pco_AIIMS_ITB_Kendall_Microbiome_Species$li[common_rows,"A1"], PCo1_Aitchison_Microbiome_Species = pco_AIIMS_ITB_Aitchison_Microbiome_Species$li[common_rows,"A1"], PCo1_Jaccard_Microbiome_Species = pco_AIIMS_ITB_Jaccard_Microbiome_Species$li[common_rows,"A1"], PCo2_Kendall_Fungome_Species =pco_AIIMS_ITB_Kendall_Fungome_Species$li[common_rows,"A2"], PCo2_Jaccard_Fungome_Species =pco_AIIMS_ITB_Jaccard_Fungome_Species$li[common_rows,"A2"], PCo3_Jaccard_Fungome_Species =pco_AIIMS_ITB_Jaccard_Fungome_Species$li[common_rows,"A3"], PCo1_Kendall_Microbiome_Genus = pco_AIIMS_ITB_Kendall_Microbiome_Genus$li[common_rows,"A1"], PCo1_Aitchison_Microbiome_Genus = pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li[common_rows,"A1"],PCo1_Jaccard_Microbiome_Genus =pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li[common_rows,"A1"], PCo1_Kendall_Fungome_Genus = pco_AIIMS_ITB_Kendall_Fungome_Genus$li[common_rows,"A1"], PCo1_Jaccard_Fungome_Genus = pco_AIIMS_ITB_Jaccard_Fungome_Genus$li[common_rows,"A1"], row.names=common_rows)

Corr_PCo_With_Markers <- corr.test(AIIMS_ITB_All_Taxa_Properties_Combined[common_rows,OnlyPCoAssociatedMarkers],df_SelectPCo_All_Taxa_Combined,method="kendall",adjust="fdr",use="pairwise.complete")

mat2 <- t(Corr_PCo_With_Markers$r[rownames(hmp_PCoFeaturePatterns$carpet),])
mat3 <- t(Corr_PCo_With_Markers$p.adj[rownames(hmp_PCoFeaturePatterns$carpet),])

#heatmap.2(mat2,density="none",trace="none",Colv=FALSE,margins=c(15,10),lwid=c(1,5),col=brewer.pal(8,"RdYlBu"),Rowv=FALSE,sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat2),rowsep=0:nrow(mat2),cellnote=apply(mat3,2,function(x)(ifelse(x<=0.2,"*",""))),notecol="darkorchid4")

save.image(paste0(directory_name,"Stage3_Workflow_Output.RData"))
