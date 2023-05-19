library(ade4)
library(vegan)
library(randomForest)
library(pcaPP)
library(pROC)
library(adegraphics)
library(randomForestExplainer)
NJC19_Matrix <- read.table("C:\\Projects\\CurrentProjects\\CardioGut\\NJC19_Matrix.txt",sep="\t",row.names=1,header=TRUE)
colnames(NJC19_Matrix) <- sub("Production.","Production",colnames(NJC19_Matrix))
colnames(NJC19_Matrix) <- sub("Consumption.","Consumption",colnames(NJC19_Matrix))
NJC19_Matrix <- NJC19_Matrix[,grep("Consumption..Production",colnames(NJC19_Matrix),value=TRUE,invert=TRUE)]
NJC19_Matrix <- NJC19_Matrix[,grep("Consumption.Production",colnames(NJC19_Matrix),value=TRUE,invert=TRUE)]

Microbiome_InferredMetabolite_AIIMS_ITB <- as.matrix(df_relab_Microbiome_Species_AIIMS_ITB[,intersect(colnames(df_relab_Microbiome_Species_AIIMS_ITB),rownames(NJC19_Matrix))]) %*% as.matrix(NJC19_Matrix[intersect(colnames(df_relab_Microbiome_Species_AIIMS_ITB),rownames(NJC19_Matrix)),])
Microbiome_InferredMetabolite_AIIMS_ITB <- Microbiome_InferredMetabolite_AIIMS_ITB[,colSums(Microbiome_InferredMetabolite_AIIMS_ITB)>0]


Kendall_MicrobiomeInferredMetabolite <- 1-cor.fk(t(Microbiome_InferredMetabolite_AIIMS_ITB))/2
pco_Kendall_MicrobiomeInferredMetabolite <- dudi.pco(as.dist(Kendall_MicrobiomeInferredMetabolite),scannf=FALSE,nf=10)
pco_Kendall_MicrobiomeInferredMetabolite$li <- as.data.frame(pco_Kendall_MicrobiomeInferredMetabolite$li)
pco_Kendall_MicrobiomeInferredMetabolite$li$Group <- NA
pco_Kendall_MicrobiomeInferredMetabolite$li[intersect(rownames(pco_Kendall_MicrobiomeInferredMetabolite$li),AIIMS_ITB_Controls),"Group"] <- "Control"
pco_Kendall_MicrobiomeInferredMetabolite$li[intersect(rownames(pco_Kendall_MicrobiomeInferredMetabolite$li),AIIMS_ITB_CD_Cases),"Group"] <- "CD"
pco_Kendall_MicrobiomeInferredMetabolite$li[intersect(rownames(pco_Kendall_MicrobiomeInferredMetabolite$li),AIIMS_ITB_Cases),"Group"] <- "ITB"

AIIMS_ITB_Inferred_Metabolite <- as.data.frame(Microbiome_InferredMetabolite_AIIMS_ITB)
AIIMS_ITB_Inferred_Metabolite$Group <- NA
AIIMS_ITB_Inferred_Metabolite[intersect(rownames(AIIMS_ITB_Inferred_Metabolite),AIIMS_ITB_Controls),"Group"] <- "Control"
AIIMS_ITB_Inferred_Metabolite[intersect(rownames(AIIMS_ITB_Inferred_Metabolite),AIIMS_ITB_CD_Cases),"Group"] <- "CD"
AIIMS_ITB_Inferred_Metabolite[intersect(rownames(AIIMS_ITB_Inferred_Metabolite),AIIMS_ITB_Cases),"Group"] <- "ITB"
CD_ITB_Inferred_Metabolite <- AIIMS_ITB_Inferred_Metabolite[AIIMS_ITB_Inferred_Metabolite$Group %in% c("CD","ITB"),]
CD_ITB_Inferred_Metabolite$Group <- ifelse(CD_ITB_Inferred_Metabolite$Group == "ITB",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(400);
	numb<-numb_features[i];
	print(numb)
	RF_CD_ITB_Metabolite_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,CD_ITB_Inferred_Metabolite);SelectedFeatures <- names(tail(RF_CD_ITB_Metabolite_Full$importance[order(RF_CD_ITB_Metabolite_Full$importance),],numb));print(numb);
	RF_CD_ITB_Metabolite_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,CD_ITB_Inferred_Metabolite[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_CD_ITB_Metabolite_Selected$y,as.numeric(RF_CD_ITB_Metabolite_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_CD_ITB_Metabolite_Selected$y,as.numeric(RF_CD_ITB_Metabolite_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}
CD_ITB_roc_vec <- roc_vec
set.seed(400);CD_ITB_Selected_Metabolites <- names(tail(RF_CD_ITB_Metabolite_Full$importance[order(RF_CD_ITB_Metabolite_Full$importance),],numb_features[which(CD_ITB_roc_vec==max(CD_ITB_roc_vec))][1]));RF_CD_ITB_Metabolite_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,CD_ITB_Inferred_Metabolite[,c(CD_ITB_Selected_Metabolites,"Group")],localImp = TRUE);
CD_ITB_Metabolites_Feature_Importance <- as.data.frame(measure_importance(RF_CD_ITB_Metabolite_Selected))
plot_multi_way_importance(CD_ITB_Metabolites_Feature_Importance, no_of_labels = length(CD_ITB_Selected_Metabolites),y_measure="times_a_root",x_measure="accuracy_decrease")

Control_ITB_Inferred_Metabolite <- AIIMS_ITB_Inferred_Metabolite[AIIMS_ITB_Inferred_Metabolite$Group %in% c("Control","ITB"),]
Control_ITB_Inferred_Metabolite$Group <- ifelse(Control_ITB_Inferred_Metabolite$Group == "ITB",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(400);
	numb<-numb_features[i];
	print(numb)
	RF_Control_ITB_Metabolite_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,Control_ITB_Inferred_Metabolite);SelectedFeatures <- names(tail(RF_Control_ITB_Metabolite_Full$importance[order(RF_Control_ITB_Metabolite_Full$importance),],numb));print(numb);
	RF_Control_ITB_Metabolite_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Control_ITB_Inferred_Metabolite[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_Control_ITB_Metabolite_Selected$y,as.numeric(RF_Control_ITB_Metabolite_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_Control_ITB_Metabolite_Selected$y,as.numeric(RF_Control_ITB_Metabolite_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}
Control_ITB_roc_vec <- roc_vec
set.seed(400);Control_ITB_Selected_Metabolites <- names(tail(RF_Control_ITB_Metabolite_Full$importance[order(RF_Control_ITB_Metabolite_Full$importance),],numb_features[which(Control_ITB_roc_vec==max(Control_ITB_roc_vec))][2]));RF_Control_ITB_Metabolite_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Control_ITB_Inferred_Metabolite[,c(Control_ITB_Selected_Metabolites,"Group")],localImp = TRUE);
Control_ITB_Metabolites_Feature_Importance <- as.data.frame(measure_importance(RF_Control_ITB_Metabolite_Selected))
plot_multi_way_importance(Control_ITB_Metabolites_Feature_Importance, no_of_labels = length(Control_ITB_Selected_Metabolites),y_measure="times_a_root",x_measure="accuracy_decrease")

Control_CD_Inferred_Metabolite <- AIIMS_ITB_Inferred_Metabolite[AIIMS_ITB_Inferred_Metabolite$Group %in% c("Control","CD"),]
Control_CD_Inferred_Metabolite$Group <- ifelse(Control_CD_Inferred_Metabolite$Group == "CD",TRUE,FALSE)

numb_features <- c(10,20,30,40,50,60,70,80,90,100,150,200,250)
roc_vec <- NULL
for(i in 1:length(numb_features))
{
	set.seed(400);
	numb<-numb_features[i];
	print(numb)
	RF_Control_CD_Metabolite_Full <- randomForest(as.factor(ifelse(Group,1,0))~.,Control_CD_Inferred_Metabolite);SelectedFeatures <- names(tail(RF_Control_CD_Metabolite_Full$importance[order(RF_Control_CD_Metabolite_Full$importance),],numb));print(numb);
	RF_Control_CD_Metabolite_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Control_CD_Inferred_Metabolite[,c(SelectedFeatures,"Group")]);
	temp_roc <- roc(RF_Control_CD_Metabolite_Selected$y,as.numeric(RF_Control_CD_Metabolite_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822");
	print(roc(RF_Control_CD_Metabolite_Selected$y,as.numeric(RF_Control_CD_Metabolite_Selected$predicted),plot=TRUE,legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE, auc.polygon = TRUE, auc.polygon.col = "#377eb822"))
	roc_vec[i] <- as.numeric(temp_roc$auc)
}
Control_CD_roc_vec <- roc_vec
set.seed(400);Control_CD_Selected_Metabolites <- names(tail(RF_Control_CD_Metabolite_Full$importance[order(RF_Control_CD_Metabolite_Full$importance),],numb_features[which(Control_CD_roc_vec==max(Control_CD_roc_vec))][2]));RF_Control_CD_Metabolite_Selected <- randomForest(as.factor(ifelse(Group,1,0))~.,Control_CD_Inferred_Metabolite[,c(Control_CD_Selected_Metabolites,"Group")],localImp = TRUE);
Control_CD_Metabolites_Feature_Importance <- as.data.frame(measure_importance(RF_Control_CD_Metabolite_Selected))
plot_multi_way_importance(Control_CD_Metabolites_Feature_Importance, no_of_labels = length(Control_CD_Selected_Metabolites),y_measure="times_a_root",x_measure="accuracy_decrease")

df_MostPredictiveFeatures_Metabolite <- data.frame(Controls_ITB = Controls_ITB_roc_vec,CD_ITB = CD_ITB_roc_vec,Controls_CD = Controls_CD_roc_vec,row.names=numb_features)

AllTopMetabolites <- unique(c(CD_ITB_Selected_Metabolites,Control_CD_Selected_Metabolites,Control_ITB_Selected_Metabolites))
CompareMetabolites_Control_ITB <- compare_features(Control_ITB_Inferred_Metabolite[,AllTopMetabolites],Control_ITB_Inferred_Metabolite$Group)
#CompareMetabolites_Control_ITB[,2] <- p.adjust(CompareMetabolites_Control_ITB[,2],method="fdr")
CompareMetabolites_Control_CD <- compare_features(Control_CD_Inferred_Metabolite[,AllTopMetabolites],Control_CD_Inferred_Metabolite$Group)
#CompareMetabolites_Control_CD[,2] <- p.adjust(CompareMetabolites_Control_CD[,2],method="fdr")
CompareMetabolites_ITB_CD <- compare_features(CD_ITB_Inferred_Metabolite[,AllTopMetabolites],CD_ITB_Inferred_Metabolite$Group)
#CompareMetabolites_ITB_CD[,2] <- p.adjust(CompareMetabolites_ITB_CD[,2],method="fdr")

df_FeaturePatterns_Metabolites <- as.data.frame(matrix(0,length(unique(c(CD_ITB_Selected_Metabolites,Control_CD_Selected_Metabolites,Control_ITB_Selected_Metabolites))),3))
rownames(df_FeaturePatterns_Metabolites) <- unique(c(CD_ITB_Selected_Metabolites,Control_CD_Selected_Metabolites,Control_ITB_Selected_Metabolites))
colnames(df_FeaturePatterns_Metabolites) <- c("Controls_v_ITB","Controls_v_CD","CD_v_ITB")
for(i in 1:nrow(df_FeaturePatterns_Metabolites))
{
	feature_name <- rownames(df_FeaturePatterns_Metabolites)[i]
	df_FeaturePatterns_Metabolites[feature_name,1] <- ifelse(p.adjust(CompareMetabolites_Control_ITB[feature_name,2],method="fdr")<=0.20,2*sign(CompareMetabolites_Control_ITB[feature_name,1]),ifelse(CompareMetabolites_Control_ITB[feature_name,2] <= 0.05,sign(CompareMetabolites_Control_ITB[feature_name,1]),0))
	df_FeaturePatterns_Metabolites[feature_name,2] <- ifelse(p.adjust(CompareMetabolites_Control_CD[feature_name,2],method="fdr")<=0.20,2*sign(CompareMetabolites_Control_CD[feature_name,1]),ifelse(CompareMetabolites_Control_CD[feature_name,2] <= 0.05,sign(CompareMetabolites_Control_CD[feature_name,1]),0))
	df_FeaturePatterns_Metabolites[feature_name,3] <- ifelse(p.adjust(CompareMetabolites_ITB_CD[feature_name,2],method="fdr")<=0.20,2*sign(CompareMetabolites_ITB_CD[feature_name,1]),ifelse(CompareMetabolites_ITB_CD[feature_name,2] <= 0.05,sign(CompareMetabolites_ITB_CD[feature_name,1]),0))
	
}

BatchDunns_InferredMetabolite <- batch_dunns(AIIMS_ITB_Inferred_Metabolite[,-ncol(AIIMS_ITB_Inferred_Metabolite)],as.factor(AIIMS_ITB_Inferred_Metabolite$Group))
mat <- BatchDunns_InferredMetabolite$final_trends
SelectMetabolites <- names(which(apply(mat,1,function(x)(length(x[x!=0])))>1))
hmpSelectedMetabolites <- heatmap.2(cor(AIIMS_ITB_Inferred_Metabolite[,SelectMetabolites],method="kendall"),density="none",trace="none",col=brewer.pal(8,"PuOr"),key.title="Kendall-Tau")

met_list <- rownames(hmpSelectedMetabolites$carpet)[(124-27):124]
for(i in 1:length(met_list))
{
	MetaboliteName <- met_list[i]
	png(paste0("C:\\Projects\\CurrentProjects\\IBD_THSTI_AIIMS\\Data\\MetaboliteFigures\\",MetaboliteName,".png"),height=480,width=480)
	boxplot(AIIMS_ITB_Inferred_Metabolite[,MetaboliteName]~factor(AIIMS_ITB_Inferred_Metabolite[,"Group"],levels=c("Control","CD","ITB")),ylab=MetaboliteName,xlab="",col=c("Cornflowerblue","Darkgoldenrod2","Firebrick1"),outline=FALSE,cex.axis=2,cex.lab=2)
	dev.off()
}

BileAcidProfiles <- data.frame(GlycoConjugatedConsumption = apply(AIIMS_ITB_Inferred_Metabolite[,GlycoConjugatedConsumption],1,mean),TauroConjugatedConsumption = apply(AIIMS_ITB_Inferred_Metabolite[,TauroConjugatedConsumption],1,mean),CA_CDCA_Production = apply(AIIMS_ITB_Inferred_Metabolite[,c("Cholic.acid_Production","Chenodeoxycholic.acid_Production")],1,mean), DCA_LCA_Production = apply(AIIMS_ITB_Inferred_Metabolite[,c("Lithocholic.acid_Production","Deoxycholic.acid_Production")],1,mean), betaMCA_Production = AIIMS_ITB_Inferred_Metabolite[,c("beta.Muricholic.acid_Production")], Oxo7LCA_Production = AIIMS_ITB_Inferred_Metabolite[,c("X7.Oxolithocholic.acid_Production")], CA_CDCA_Consumption = apply(AIIMS_ITB_Inferred_Metabolite[,c("Cholic.acid_Consumption","Chenodeoxycholic.acid_Consumption")],1,mean), Group = AIIMS_ITB_Inferred_Metabolite[,"Group"])

bile_list <- colnames(BileAcidProfiles)
for(i in 1:length(bile_list))
{
	MetaboliteName <- bile_list[i]
	png(paste0("C:\\Projects\\CurrentProjects\\IBD_THSTI_AIIMS\\Data\\MetaboliteFigures\\SelectMetabolites\\",MetaboliteName,".png"),height=480,width=480)
	boxplot(BileAcidProfiles[,MetaboliteName]~factor(BileAcidProfiles[,"Group"],levels=c("Control","CD","ITB")),ylab=MetaboliteName,xlab="",col=c("Cornflowerblue","Darkgoldenrod2","Firebrick1"),outline=FALSE,cex.axis=2,cex.lab=2)
	dev.off()
}


