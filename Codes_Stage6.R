library(psych)
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

IndianStudyNames <- c("AIIMS_1_IBD","AIIMS_ITB","DhakanDB","GuptaA","logmpie","MicroDiab","PRJDB7616")

df_clr_Microbiome_Genus_IndianStudies <- df_clr_Microbiome_Genus_Combined[df_clr_Microbiome_Genus_Combined$study_name %in% IndianStudyNames,]
df_clr_Microbiome_Species_IndianStudies <- df_clr_Microbiome_Species_Combined[df_clr_Microbiome_Species_Combined$study_name %in% IndianStudyNames,]
PRJDB7616_Japanese_Samples <- grep("^R",rownames(df_clr_Microbiome_Species_IndianStudies),value=TRUE)
df_clr_Microbiome_Species_IndianStudies <- df_clr_Microbiome_Species_IndianStudies[setdiff(rownames(df_clr_Microbiome_Species_IndianStudies),PRJDB7616_Japanese_Samples),]
df_clr_Microbiome_Genus_IndianStudies <- df_clr_Microbiome_Genus_IndianStudies[setdiff(rownames(df_clr_Microbiome_Genus_IndianStudies),PRJDB7616_Japanese_Samples),]

## Loading Metadata and selecting healthy samples
MicroDiab_Metadata <- read.table(paste0(directory_name,"Microdiab_Metadata.txt"),sep="\t",row.names=1,header=TRUE)
MicroDiab_Healthy <- rownames(MicroDiab_Metadata[MicroDiab_Metadata[,3] == "NG",])
GuptaA_Metadata <- read.table(paste0(directory_name,"GuptaA_Metadata.txt"),sep="\t",row.names=1,header=TRUE)
GuptaA_Healthy <- rownames(GuptaA_Metadata[GuptaA_Metadata[,3] == "control",])

DhakanDB_Healthy <- rownames(df_clr_Microbiome_Species_IndianStudies[df_clr_Microbiome_Species_IndianStudies$study_name == "DhakanDB",])
AIIMS_1_IBD_Healthy <- c(AIIMS_1_IBD_Controls,AIIMS_1_IBD_Leh_Controls)
AIIMS_ITB_Healthy <- AIIMS_ITB_Controls
logmpie_Healthy <- rownames(df_clr_Microbiome_Species_IndianStudies[df_clr_Microbiome_Species_IndianStudies$study_name == "logmpie",])
PRJDB7616_Healthy <- setdiff(rownames(df_clr_Microbiome_Species_IndianStudies[df_clr_Microbiome_Species_IndianStudies$study_name == "PRJDB7616",]),PRJDB7616_Japanese_Samples)

All_Indian_Healthy <- c(MicroDiab_Healthy,GuptaA_Healthy,DhakanDB_Healthy,AIIMS_1_IBD_Healthy,AIIMS_ITB_Healthy,logmpie_Healthy,PRJDB7616_Healthy)

df_clr_Microbiome_Species_IndianStudies_Controls <- df_clr_Microbiome_Species_IndianStudies[All_Indian_Healthy,]
df_clr_Microbiome_Genus_IndianStudies_Controls <- df_clr_Microbiome_Genus_IndianStudies[All_Indian_Healthy,]

SpeciesDetectionIndianStudyCohorts <- compute_detection(df_clr_Microbiome_Species_IndianStudies_Controls,colnames(df_clr_Microbiome_Species_IndianStudies)[1:4675],"study_name",IndianStudyNames)
GenusDetectionIndianStudyCohorts <- compute_detection(df_clr_Microbiome_Genus_IndianStudies_Controls,colnames(df_clr_Microbiome_Genus_IndianStudies)[1:1526],"study_name",IndianStudyNames)

select_species <- names(which(apply(SpeciesDetectionIndianStudyCohorts,1,function(x)(length(x[x>=0.2])))>=2))
select_genus <- names(which(apply(GenusDetectionIndianStudyCohorts,1,function(x)(length(x[x>=0.2])))>=2))

df_clr_Microbiome_All_IndianStudies_Controls <- as.data.frame(cbind(df_clr_Microbiome_Species_IndianStudies_Controls[rownames(df_clr_Microbiome_Species_IndianStudies_Controls),select_species],df_clr_Microbiome_Genus_IndianStudies_Controls[rownames(df_clr_Microbiome_Species_IndianStudies_Controls),c(select_genus,"study_name")]))

IndianControlAllTaxaMetwork <- rem_network1(df_clr_Microbiome_All_IndianStudies_Controls,c(select_species,select_genus),"study_name",IndianStudyNames,0.001,0.7)

AllDetectionIndianStudyCohorts <- compute_detection(df_clr_Microbiome_All_IndianStudies_Controls,colnames(df_clr_Microbiome_All_IndianStudies_Controls)[1:365],"study_name",IndianStudyNames)

CoOccurrence_IndianControlAllTaxaMetwork <- graph_from_adjacency_matrix(as.matrix(IndianControlAllTaxaMetwork$dir))

DetectedCore <- names(which(apply(AllDetectionIndianStudyCohorts,1,function(x)(length(x[x>=0.7])))>=4))
DegreeCentralCore <- names(which(rank_scale(degree(graph_from_adjacency_matrix(as.matrix(IndianControlAllTaxaMetwork$dir))))>=0.70))
BetweennessCentralCore <- names(which(rank_scale(betweenness(graph_from_adjacency_matrix(as.matrix(IndianControlSpeciesMetwork$dir))))>=0.70))
DegreeDetectedCore <- names(which(table(c(DetectedCore,DegreeCentralCore))==2))
subgraph_ControlCore_IndianControlAllTaxaMetwork <- subgraph(CoOccurrence_IndianControlAllTaxaMetwork,v=DegreeDetectedCore)
as_edgelist()

df_clr_Microbiome_All_IndianStudies <- as.data.frame(cbind(df_clr_Microbiome_Species_IndianStudies[rownames(df_clr_Microbiome_Species_IndianStudies),select_species],df_clr_Microbiome_Genus_IndianStudies[rownames(df_clr_Microbiome_Species_IndianStudies),c(select_genus,"study_name")]))
IndianAllTaxaMetwork <- rem_network1(df_clr_Microbiome_All_IndianStudies,c(select_species,select_genus),"study_name",IndianStudyNames,0.001,0.7)
AllSamplesDetectionIndianStudyCohorts <- compute_detection(df_clr_Microbiome_All_IndianStudies,colnames(df_clr_Microbiome_All_IndianStudies)[1:365],"study_name",IndianStudyNames)
DetectedCoreAll <- names(which(apply(AllSamplesDetectionIndianStudyCohorts,1,function(x)(length(x[x>=0.7])))>=4))
DegreeCentralCoreAll <- names(which(rank_scale(degree(graph_from_adjacency_matrix(as.matrix(IndianAllTaxaMetwork$dir))))>=0.70))
BetweennessCentralCoreAll <- names(which(rank_scale(betweenness(graph_from_adjacency_matrix(as.matrix(IndianAllTaxaMetwork$dir))))>=0.70))
DegreeDetectedCoreAll<- names(which(table(c(DetectedCoreAll,DegreeCentralCoreAll))==2))
CoOccurrence_IndianAllTaxaMetwork <- graph_from_adjacency_matrix(as.matrix(IndianAllTaxaMetwork$dir))
subgraph_All_IndianAllTaxaMetwork <- subgraph(CoOccurrence_IndianControlAllTaxaMetwork,v=intersect(DegreeDetectedCoreAll,colnames(df_clr_Microbiome_Species_IndianStudies)))
write.table(as_edgelist(subgraph_All_IndianAllTaxaMetwork),paste0(directory_name,"raw_edgelist_Core_IndianAllTaxaMetwork"),sep="\t")




DegreeDetectedCoreProperties <- read.table("C:\\Projects\\CurrentProjects\\IBD_THSTI_AIIMS\\Data\\DegreeDetectedCoreProperties.txt",sep="\t",row.names=1,header=TRUE)
df_relab_Microbiome_Species_AIIMS_ITB <- df_Microbiome_Species_AIIMS_ITB/rowSums(df_Microbiome_Species_AIIMS_ITB)
df_relab_Microbiome_Species_AIIMS_1_IBD <- df_Microbiome_Species_AIIMS_1_IBD/rowSums(df_Microbiome_Species_AIIMS_1_IBD)
CoreScores_Species_AIIMS_ITB <- as.data.frame(apply(df_relab_Microbiome_Species_AIIMS_ITB[,intersect(colnames(df_relab_Microbiome_Species_AIIMS_ITB),rownames(DegreeDetectedCoreProperties))],2,rank_scale) %*% as.matrix(DegreeDetectedCoreProperties[intersect(colnames(df_relab_Microbiome_Species_AIIMS_ITB),rownames(DegreeDetectedCoreProperties)),]))
CoreScores_Species_AIIMS_ITB$Type <- NA
CoreScores_Species_AIIMS_ITB[intersect(rownames(CoreScores_Species_AIIMS_ITB),AIIMS_ITB_Controls),"Type"] <- "Controls"
CoreScores_Species_AIIMS_ITB[intersect(rownames(CoreScores_Species_AIIMS_ITB),AIIMS_ITB_Cases),"Type"] <- "ITB"
CoreScores_Species_AIIMS_ITB[intersect(rownames(CoreScores_Species_AIIMS_ITB),AIIMS_ITB_CD_Cases),"Type"] <- "CD"
CoreScores_Species_AIIMS_1_IBD <- as.data.frame(apply(df_relab_Microbiome_Species_AIIMS_1_IBD[,intersect(colnames(df_relab_Microbiome_Species_AIIMS_1_IBD),rownames(DegreeDetectedCoreProperties))],2,rank_scale) %*% as.matrix(DegreeDetectedCoreProperties[intersect(colnames(df_relab_Microbiome_Species_AIIMS_1_IBD),rownames(DegreeDetectedCoreProperties)),]))
CoreScores_Species_AIIMS_1_IBD$Type <- NA
CoreScores_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_IBD_Controls),"Type"] <- "Controls"
CoreScores_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_IBD_Leh_Controls),"Type"] <- "Controls"
CoreScores_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_IBD_Cases),"Type"] <- "IBD"
CoreScores_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_ASUC_Cases),"Type"] <- "ASUC"
DegreeDetectedCorePropertiesAll <- read.table("C:\\Projects\\CurrentProjects\\IBD_THSTI_AIIMS\\Data\\DegreeDetectedCorePropertiesAll.txt",sep="\t",row.names=1,header=TRUE)
CoreScores_All_Species_AIIMS_ITB <- as.data.frame(apply(df_relab_Microbiome_Species_AIIMS_ITB[,intersect(colnames(df_relab_Microbiome_Species_AIIMS_ITB),rownames(DegreeDetectedCorePropertiesAll))],2,rank_scale) %*% as.matrix(DegreeDetectedCorePropertiesAll[intersect(colnames(df_relab_Microbiome_Species_AIIMS_ITB),rownames(DegreeDetectedCorePropertiesAll)),]))
CoreScores_All_Species_AIIMS_ITB$Type <- NA
CoreScores_All_Species_AIIMS_ITB[intersect(rownames(CoreScores_All_Species_AIIMS_ITB),AIIMS_ITB_Controls),"Type"] <- "Controls"
CoreScores_All_Species_AIIMS_ITB[intersect(rownames(CoreScores_All_Species_AIIMS_ITB),AIIMS_ITB_Cases),"Type"] <- "ITB"
CoreScores_All_Species_AIIMS_ITB[intersect(rownames(CoreScores_All_Species_AIIMS_ITB),AIIMS_ITB_CD_Cases),"Type"] <- "CD"
CoreScores_All_Species_AIIMS_1_IBD <- as.data.frame(apply(df_relab_Microbiome_Species_AIIMS_1_IBD[,intersect(colnames(df_relab_Microbiome_Species_AIIMS_1_IBD),rownames(DegreeDetectedCorePropertiesAll))],2,rank_scale) %*% as.matrix(DegreeDetectedCorePropertiesAll[intersect(colnames(df_relab_Microbiome_Species_AIIMS_1_IBD),rownames(DegreeDetectedCorePropertiesAll)),]))
CoreScores_All_Species_AIIMS_1_IBD$Type <- NA
CoreScores_All_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_IBD_Controls),"Type"] <- "Controls"
CoreScores_All_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_IBD_Leh_Controls),"Type"] <- "Controls"
CoreScores_All_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_IBD_Cases),"Type"] <- "IBD"
CoreScores_All_Species_AIIMS_1_IBD[intersect(rownames(CoreScores_Species_AIIMS_1_IBD),AIIMS_1_ASUC_Cases),"Type"] <- "ASUC"
CoreScores_All_Species_MicroDiab <- as.data.frame(apply(df_relab_Microbiome_Species_MicroDiab[,intersect(colnames(df_relab_Microbiome_Species_MicroDiab),rownames(DegreeDetectedCorePropertiesAll))],2,rank_scale) %*% as.matrix(DegreeDetectedCorePropertiesAll[intersect(colnames(df_relab_Microbiome_Species_MicroDiab),rownames(DegreeDetectedCorePropertiesAll)),]))
CoreScores_All_Species_GuptaA <- as.data.frame(apply(GuptaA_2019_Species[,intersect(colnames(GuptaA_2019_Species),rownames(DegreeDetectedCorePropertiesAll))],2,rank_scale) %*% as.matrix(DegreeDetectedCorePropertiesAll[intersect(colnames(GuptaA_2019_Species),rownames(DegreeDetectedCorePropertiesAll)),]))