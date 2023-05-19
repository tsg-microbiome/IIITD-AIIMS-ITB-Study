FungomeSpeciesStudies <- sub("df_clr_Fungome_Species_","",ls(pattern="df_clr_Fungome_Species"))

df1 <- get(paste0("df_clr_Fungome_Species_",FungomeSpeciesStudies[1]))
df2 <- get(paste0("df_clr_Fungome_Species_",FungomeSpeciesStudies[2]))

df_combined <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,-1]
rownames(df_combined) <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,1]
df_combined <- apply(df_combined,1,function(x)(ifelse(is.na(x),0,x)))

cat("Combining Fungome Species")
for(i in 3:length(FungomeSpeciesStudies))
{
	StudyName <- FungomeSpeciesStudies[i]
	df_temp <- get(paste0("df_clr_Fungome_Species_",StudyName))
	temp_combined <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,-1]
	rownames(temp_combined) <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,1]
	temp_combined <- apply(temp_combined,1,function(x)(ifelse(is.na(x),0,x)))
	df_combined <- temp_combined
}
df_clr_Fungome_Species_Combined <- as.data.frame(df_combined)
df_clr_Fungome_Species_Combined$study_name <- NA
for(i in 1:length(FungomeSpeciesStudies))
{
	StudyName <- FungomeSpeciesStudies[i]
	df_temp <- get(paste0("df_clr_Fungome_Species_",StudyName))
	df_clr_Fungome_Species_Combined[rownames(df_temp),"study_name"] <- StudyName
}

rm(df_temp)
rm(df_combined)

FungomeGenusStudies <- sub("df_clr_Fungome_Genus_","",ls(pattern="df_clr_Fungome_Genus"))

df1 <- get(paste0("df_clr_Fungome_Genus_",FungomeSpeciesStudies[1]))
df2 <- get(paste0("df_clr_Fungome_Genus_",FungomeSpeciesStudies[2]))

df_combined <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,-1]
rownames(df_combined) <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,1]
df_combined <- apply(df_combined,1,function(x)(ifelse(is.na(x),0,x)))

cat("Combining Fungome Genus\n")
for(i in 3:length(FungomeGenusStudies))
{
	StudyName <- FungomeGenusStudies[i]
	df_temp <- get(paste0("df_clr_Fungome_Genus_",StudyName))
	temp_combined <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,-1]
	rownames(temp_combined) <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,1]
	temp_combined <- apply(temp_combined,1,function(x)(ifelse(is.na(x),0,x)))
	df_combined <- temp_combined
}
df_clr_Fungome_Genus_Combined <- as.data.frame(df_combined)
df_clr_Fungome_Genus_Combined$study_name <- NA
for(i in 1:length(FungomeGenusStudies))
{
	StudyName <- FungomeGenusStudies[i]
	df_temp <- get(paste0("df_clr_Fungome_Genus_",StudyName))
	df_clr_Fungome_Genus_Combined[rownames(df_temp),"study_name"] <- StudyName
}
rm(df_temp)
rm(df_combined)

MicrobiomeSpeciesStudies <- sub("df_clr_Microbiome_Species_","",ls(pattern="df_clr_Microbiome_Species"))
MicrobiomeSpeciesStudies <- setdiff(MicrobiomeSpeciesStudies,"Chinese_ITB")
df1 <- get(paste0("df_clr_Microbiome_Species_",MicrobiomeSpeciesStudies[1]))
df2 <- get(paste0("df_clr_Microbiome_Species_",MicrobiomeSpeciesStudies[2]))
df_combined <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,-1]
rownames(df_combined) <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,1]
df_combined <- apply(df_combined,1,function(x)(ifelse(is.na(x),0,x)))

cat("Combining Microbiome Species\n")
for(i in 3:length(MicrobiomeSpeciesStudies))
{
	StudyName <- MicrobiomeSpeciesStudies[i]
	df_temp <- get(paste0("df_clr_Microbiome_Species_",StudyName))
	print(paste0("df_clr_Microbiome_Species_",StudyName))
	temp_combined <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,-1]
	rownames(temp_combined) <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,1]
	temp_combined <- apply(temp_combined,1,function(x)(ifelse(is.na(x),0,x)))
	df_combined <- temp_combined
}
df_clr_Microbiome_Species_Combined <- as.data.frame(df_combined)
df_clr_Microbiome_Species_Combined$study_name <- NA
for(i in 1:length(MicrobiomeSpeciesStudies))
{
	StudyName <- MicrobiomeSpeciesStudies[i]
	df_temp <- get(paste0("df_clr_Microbiome_Species_",StudyName))
	df_clr_Microbiome_Species_Combined[rownames(df_temp),"study_name"] <- StudyName
}
rm(df_temp)
rm(df_combined)

MicrobiomeGenusStudies <- sub("df_clr_Microbiome_Genus_","",ls(pattern="df_clr_Microbiome_Genus"))

df1 <- get(paste0("df_clr_Microbiome_Genus_",MicrobiomeGenusStudies[1]))
df2 <- get(paste0("df_clr_Microbiome_Genus_",MicrobiomeGenusStudies[2]))

df_combined <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,-1]
rownames(df_combined) <- merge(t(df1),t(df2),by="row.names",all=TRUE)[,1]
df_combined <- apply(df_combined,1,function(x)(ifelse(is.na(x),0,x)))

cat("Combining Microbiome Genus\n")
for(i in 3:length(MicrobiomeGenusStudies))
{
	StudyName <- MicrobiomeGenusStudies[i]
	df_temp <- get(paste0("df_clr_Microbiome_Genus_",StudyName))
	print(paste0("df_clr_Microbiome_Genus_",StudyName))
	temp_combined <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,-1]
	rownames(temp_combined) <- merge(t(df_combined),t(df_temp),by="row.names",all=TRUE)[,1]
	temp_combined <- apply(temp_combined,1,function(x)(ifelse(is.na(x),0,x)))
	df_combined <- temp_combined
}
df_clr_Microbiome_Genus_Combined <- as.data.frame(df_combined)
df_clr_Microbiome_Genus_Combined$study_name <- NA
for(i in 1:length(MicrobiomeGenusStudies))
{
	StudyName <- MicrobiomeGenusStudies[i]
	df_temp <- get(paste0("df_clr_Microbiome_Genus_",StudyName))
	df_clr_Microbiome_Genus_Combined[rownames(df_temp),"study_name"] <- StudyName
}
rm(df_temp)
rm(df_combined)

cat("Generating MetaNetworks\n")
AllUnionMarkerAbundance <- sub("abund_","",grep("abund_",union(AllTopMarkers,OnlyPCoAssociatedMarkers),value=TRUE))

combined_microbiome_rows <- intersect(rownames(df_clr_Microbiome_Species_Combined),rownames(df_clr_Microbiome_Genus_Combined))
df_Microbiome_MetaNetwork <- as.data.frame(cbind(df_clr_Microbiome_Species_Combined[combined_microbiome_rows,],df_clr_Microbiome_Genus_Combined[combined_microbiome_rows,]))
Microbiome_MetaNetwork <- rem_network1(df_Microbiome_MetaNetwork,intersect(AllUnionMarkerAbundance,colnames(df_Microbiome_MetaNetwork)),"study_name",setdiff(MicrobiomeSpeciesStudies,c("Chinese_ITB","AIIMS_ITB")),0.0001,0.7)
co_occurrence_MetaNetwork_Microbiome <- graph_from_adjacency_matrix(as.matrix(Microbiome_MetaNetwork$dir))

combined_Fungome_rows <- intersect(rownames(df_clr_Fungome_Species_Combined),rownames(df_clr_Fungome_Genus_Combined))
df_Fungome_MetaNetwork <- as.data.frame(cbind(df_clr_Fungome_Species_Combined[combined_Fungome_rows,],df_clr_Fungome_Genus_Combined[combined_Fungome_rows,]))
Fungome_MetaNetwork <- rem_network1(df_Fungome_MetaNetwork,intersect(AllUnionMarkerAbundance,colnames(df_Fungome_MetaNetwork)),"study_name",setdiff(FungomeSpeciesStudies,"AIIMS_ITB"),0.05,0.7)
co_occurrence_MetaNetwork_Fungome <- graph_from_adjacency_matrix(as.matrix(Fungome_MetaNetwork$dir))

write.table(as_edgelist(co_occurrence_MetaNetwork_Fungome),file=paste0(directory_name,"Networks/FungomeMetaNetwork_EdgelistRaw.txt"),sep="\t",row.names=TRUE,col.names=TRUE)
write.table(as_edgelist(co_occurrence_MetaNetwork_Microbiome),file=paste0(directory_name,"Networks/MicrobiomeMetaNetwork_EdgelistRaw.txt"),sep="\t",row.names=TRUE,col.names=TRUE)
## Reprocess the edgelists and then start the Stage 5.

save.image(paste0(directory_name,"Stage4_Workflow_Output.RData"))

