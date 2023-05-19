library(ccrepe)
print("Generating CCREPE for Markers")
load(paste0(directory_name,"Stage4_Workflow_Output.RData"))

df_All_Taxa_Properties_for_ccrepe <- AIIMS_ITB_All_Taxa_Properties_Combined[,intersect(colnames(AIIMS_ITB_All_Taxa_Properties_Combined),paste0("abund_",AllUnionMarkerAbundance))]

colnames(df_All_Taxa_Properties_for_ccrepe) <- sub("abund_","",colnames(df_All_Taxa_Properties_for_ccrepe))

ccrepe_AIIMS_ITB_All_Markers <- ccrepe(df_All_Taxa_Properties_for_ccrepe,sim.score=cor,sim.score.args=list(method="kendall",use="pairwise.complete"))

ccrepe_AIIMS_ITB_All_Markers$sim.score <- apply(ccrepe_AIIMS_ITB_All_Markers$sim.score,2,function(x)(ifelse(is.na(x),0,x)))
ccrepe_AIIMS_ITB_All_Markers$p.value <- apply(ccrepe_AIIMS_ITB_All_Markers$p.value,2,function(x)(ifelse(is.na(x),1,x)))
ccrepe_AIIMS_ITB_All_Markers$q.value <- apply(ccrepe_AIIMS_ITB_All_Markers$p.value,2,function(x)(p.adjust(x,method="fdr")))
ccrepe_AIIMS_ITB_All_Markers$dir <- matrix(NA,nrow(ccrepe_AIIMS_ITB_All_Markers$q.value),ncol(ccrepe_AIIMS_ITB_All_Markers$q.value))
rownames(ccrepe_AIIMS_ITB_All_Markers$dir) <- rownames(ccrepe_AIIMS_ITB_All_Markers$q.value)
colnames(ccrepe_AIIMS_ITB_All_Markers$dir) <- colnames(ccrepe_AIIMS_ITB_All_Markers$q.value)
for(i in 1:nrow(ccrepe_AIIMS_ITB_All_Markers$q.value))
{
	for(j in 1:ncol(ccrepe_AIIMS_ITB_All_Markers$q.value))
	{
		ccrepe_AIIMS_ITB_All_Markers$dir[i,j] <- ifelse((ccrepe_AIIMS_ITB_All_Markers$sim.score[i,j] > 0)&(ccrepe_AIIMS_ITB_All_Markers$q.value[i,j] <= 0.10),1,0)
	}
}

co_occurrence_All_AIIMS_ITB <- graph_from_adjacency_matrix(ccrepe_AIIMS_ITB_All_Markers$dir)

save.image(paste0(directory_name,"Stage5_Workflow_Output.RData"))