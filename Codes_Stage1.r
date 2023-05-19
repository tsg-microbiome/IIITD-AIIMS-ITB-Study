print("Stage 1: Reading Data \n")
directory_name <- "C:\\Projects\\CurrentProjects\\IBD_THSTI_AIIMS\\Data\\"

source(paste0(directory_name,"function_library.R"))

df_ASV_AIIMS_ITB <- read.table(paste0(directory_name,"ASV-table.tsv"),sep="\t",row.names=1,header=TRUE)
df_Spingo_AIIMS_ITB <- read.table(paste0(directory_name,"ITB_AIIMS_rep_seq.spingo.out.txt"),sep="\t",row.names=1,header=TRUE)
df_Spingo_AIIMS_ITB <- df_Spingo_AIIMS_ITB[,c(4,5)]
df_Spingo_AIIMS_ITB[,1] <- ifelse(df_Spingo_AIIMS_ITB[,1] == "Escherichia/Shigella_coli","Escherichia_coli",df_Spingo_AIIMS_ITB[,1])
df_Spingo_AIIMS_ITB[,1] <- ifelse(df_Spingo_AIIMS_ITB[,2] < 0.70 ,"AMBIGUOUS",df_Spingo_AIIMS_ITB[,1])
df_Spingo_AIIMS_ITB <- df_Spingo_AIIMS_ITB[intersect(rownames(df_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),]
df_Genus_Spingo_AIIMS_ITB <- read.table(paste0(directory_name,"ITB_AIIMS_rep_seq.genus.spingo.txt"),sep="\t",row.names=1,header=TRUE)
df_Genus_Spingo_AIIMS_ITB <- df_Genus_Spingo_AIIMS_ITB[,c(4,5)]
df_Genus_Spingo_AIIMS_ITB[,1] <- ifelse(df_Genus_Spingo_AIIMS_ITB[,1] == "Escherichia/Shigella","Escherichia",df_Genus_Spingo_AIIMS_ITB[,1])

df_ASV_AIIMS_ITB <- df_ASV_AIIMS_ITB[intersect(rownames(df_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),]

df_Microbiome_Species_AIIMS_ITB <- aggregate(df_ASV_AIIMS_ITB,by=list(df_Spingo_AIIMS_ITB[intersect(rownames(df_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),1]),FUN=sum)[,-1]
rownames(df_Microbiome_Species_AIIMS_ITB) <- aggregate(df_ASV_AIIMS_ITB,by=list(df_Spingo_AIIMS_ITB[intersect(rownames(df_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),1]),FUN=sum)[,1]
df_Microbiome_Species_AIIMS_ITB <- as.data.frame(t(df_Microbiome_Species_AIIMS_ITB))
df_Microbiome_Species_AIIMS_ITB <- df_Microbiome_Species_AIIMS_ITB[,setdiff(colnames(df_Microbiome_Species_AIIMS_ITB),"AMBIGUOUS")]
df_clr_Microbiome_Species_AIIMS_ITB <- as.data.frame(as.matrix(clr(df_Microbiome_Species_AIIMS_ITB+0.00001)))
df_clr_Microbiome_Species_AIIMS_ITB <- as.data.frame(t(apply(df_clr_Microbiome_Species_AIIMS_ITB,1,function(x)(x-min(x)))))

df_ASV_AIIMS_ITB <- df_ASV_AIIMS_ITB[intersect(rownames(df_Genus_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),]
df_Microbiome_Genus_AIIMS_ITB <- aggregate(df_ASV_AIIMS_ITB,by=list(df_Genus_Spingo_AIIMS_ITB[intersect(rownames(df_Genus_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),1]),FUN=sum)[,-1]
rownames(df_Microbiome_Genus_AIIMS_ITB) <- aggregate(df_ASV_AIIMS_ITB,by=list(df_Genus_Spingo_AIIMS_ITB[intersect(rownames(df_Genus_Spingo_AIIMS_ITB),rownames(df_ASV_AIIMS_ITB)),1]),FUN=sum)[,1]
df_Microbiome_Genus_AIIMS_ITB <- as.data.frame(t(df_Microbiome_Genus_AIIMS_ITB))
df_Microbiome_Genus_AIIMS_ITB <- df_Microbiome_Genus_AIIMS_ITB[,setdiff(colnames(df_Microbiome_Genus_AIIMS_ITB),"AMBIGUOUS")]
df_clr_Microbiome_Genus_AIIMS_ITB <- as.data.frame(as.matrix(clr(df_Microbiome_Genus_AIIMS_ITB+0.00001)))
df_clr_Microbiome_Genus_AIIMS_ITB <- as.data.frame(t(apply(df_clr_Microbiome_Genus_AIIMS_ITB,1,function(x)(x-min(x)))))

df_Metadata_AIIMS_ITB <- read.table(paste0(directory_name,"metadata-modified-final.tsv"),sep="\t",row.names=1,header=TRUE)
rownames(df_Metadata_AIIMS_ITB) <- sub("sample-","itb_sample_",rownames(df_Metadata_AIIMS_ITB))
AIIMS_ITB_Controls <- rownames(df_Metadata_AIIMS_ITB[df_Metadata_AIIMS_ITB[,1]=="Controls",])
AIIMS_ITB_Cases <- rownames(df_Metadata_AIIMS_ITB[df_Metadata_AIIMS_ITB[,1]=="ITB",])
AIIMS_ITB_CD_Cases <- rownames(df_Metadata_AIIMS_ITB[df_Metadata_AIIMS_ITB[,1]=="CD",])


##Reading PRJDB7616 Bacteriome data
print("PRJDB7616_Bacteriome_data")
df_Microbiome_Species_PRJDB7616 <- read.table(paste0(directory_name,"PRJDB7616_16S.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_PRJDB7616 <- as.data.frame(clr(df_Microbiome_Species_PRJDB7616+0.00001))
df_clr_Microbiome_Species_PRJDB7616 <- as.data.frame(t(apply(df_clr_Microbiome_Species_PRJDB7616,1,function(x)(x-min(x)))))

df_Microbiome_Genus_PRJDB7616 <- read.table(paste0(directory_name,"PRJDB7616_16S.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_PRJDB7616 <- as.data.frame(clr(df_Microbiome_Genus_PRJDB7616+0.00001))
df_clr_Microbiome_Genus_PRJDB7616 <- as.data.frame(t(apply(df_clr_Microbiome_Genus_PRJDB7616,1,function(x)(x-min(x)))))

##Reading Liguori_et_al Bacteriome data
print("Liguori_et_al_Bacteriome_data")
df_Microbiome_Species_Liguori <- read.table(paste0(directory_name,"Liguori_16S_Species_Table.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_Liguori <- as.data.frame(clr(df_Microbiome_Species_Liguori+0.00001))
df_clr_Microbiome_Species_Liguori <- as.data.frame(t(apply(df_clr_Microbiome_Species_Liguori,1,function(x)(x-min(x)))))

df_Microbiome_Genus_Liguori <- read.table(paste0(directory_name,"Liguori_16S_Genus_Table.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_Liguori <- as.data.frame(clr(df_Microbiome_Genus_Liguori+0.00001))
df_clr_Microbiome_Genus_Liguori <- as.data.frame(t(apply(df_clr_Microbiome_Genus_Liguori,1,function(x)(x-min(x)))))

##Reading PRJEB42375 Bacteriome data
print("PRJEB42375_Bacteriome_data")
df_Microbiome_Species_PRJEB42375 <- read.table(paste0(directory_name,"PRJEB42375_16S.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_PRJEB42375 <- as.data.frame(clr(df_Microbiome_Species_PRJEB42375+0.00001))
df_clr_Microbiome_Species_PRJEB42375 <- as.data.frame(t(apply(df_clr_Microbiome_Species_PRJEB42375,1,function(x)(x-min(x)))))

df_Microbiome_Genus_PRJEB42375 <- read.table(paste0(directory_name,"PRJEB42375_16S.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_PRJEB42375 <- as.data.frame(clr(df_Microbiome_Genus_PRJEB42375+0.00001))
df_clr_Microbiome_Genus_PRJEB42375 <- as.data.frame(t(apply(df_clr_Microbiome_Genus_PRJEB42375,1,function(x)(x-min(x)))))


## Reading Batch 1 IBD Data (Kedia Ghosh paper and Das Ghosh paper)##
df_Microbiome_Species_AIIMS_1_IBD <- read.table(paste0(directory_name,"IBD_v1_Species_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_AIIMS_1_IBD <- as.data.frame(as.matrix(clr(df_Microbiome_Species_AIIMS_1_IBD+0.00001)))
df_clr_Microbiome_Species_AIIMS_1_IBD <- as.data.frame(t(apply(df_clr_Microbiome_Species_AIIMS_1_IBD,1,function(x)(x-min(x)))))

df_Microbiome_Genus_AIIMS_1_IBD <- read.table(paste0(directory_name,"IBD_v1_Genus_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_AIIMS_1_IBD <- as.data.frame(as.matrix(clr(df_Microbiome_Genus_AIIMS_1_IBD+0.00001)))
df_clr_Microbiome_Genus_AIIMS_1_IBD <- as.data.frame(t(apply(df_clr_Microbiome_Genus_AIIMS_1_IBD,1,function(x)(x-min(x)))))

df_Metadata_AIIMS_1_IBD <- read.table(paste0(directory_name,"IBD_AIIMS_1_metadata.txt"),sep="\t",row.names=1,header=TRUE)

AIIMS_1_IBD_Controls <- rownames(df_Metadata_AIIMS_1_IBD[(df_Metadata_AIIMS_1_IBD[,1]!="Healthy_Leh")&(df_Metadata_AIIMS_1_IBD[,2]=="Healthy"),])
AIIMS_1_IBD_Leh_Controls <- rownames(df_Metadata_AIIMS_1_IBD[(df_Metadata_AIIMS_1_IBD[,1]=="Healthy_Leh")&(df_Metadata_AIIMS_1_IBD[,2]=="Healthy"),])
AIIMS_1_IBD_Cases <- rownames(df_Metadata_AIIMS_1_IBD[(df_Metadata_AIIMS_1_IBD[,2]=="IBD"),])
AIIMS_1_ASUC_Cases <- rownames(df_Metadata_AIIMS_1_IBD[(df_Metadata_AIIMS_1_IBD[,2]=="ASUC"),])

## Reading ITB Chinese
df_Microbiome_Species_Chinese_ITB <- read.table(paste0(directory_name,"Chinese_ITB_Species_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
rownames(df_Microbiome_Species_Chinese_ITB) <- sub(".spingo.txt","",rownames(df_Microbiome_Species_Chinese_ITB))
df_clr_Microbiome_Species_Chinese_ITB <- as.data.frame(as.matrix(clr(df_Microbiome_Species_Chinese_ITB+0.00001)))
df_clr_Microbiome_Species_Chinese_ITB <- as.data.frame(t(apply(df_clr_Microbiome_Species_Chinese_ITB,1,function(x)(x-min(x)))))

df_Microbiome_Genus_Chinese_ITB <- read.table(paste0(directory_name,"Chinese_ITB_Genus_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
rownames(df_Microbiome_Genus_Chinese_ITB) <- sub(".spingo.txt","",rownames(df_Microbiome_Genus_Chinese_ITB))
df_clr_Microbiome_Genus_Chinese_ITB <- as.data.frame(as.matrix(clr(df_Microbiome_Genus_Chinese_ITB+0.00001)))
df_clr_Microbiome_Genus_Chinese_ITB <- as.data.frame(t(apply(df_clr_Microbiome_Genus_Chinese_ITB,1,function(x)(x-min(x)))))

Chinese_ITB_Controls <- grep("HC_",rownames(df_clr_Microbiome_Species_Chinese_ITB),value=TRUE)
Chinese_ITB_Cases <- grep("ITB_",rownames(df_clr_Microbiome_Species_Chinese_ITB),value=TRUE)

## Reading He et al
df_Microbiome_Species_He <- read.table(paste0(directory_name,"He_Controls_Species_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_He <- as.data.frame(as.matrix(clr(df_Microbiome_Species_He+0.00001)))
df_clr_Microbiome_Species_He <- as.data.frame(t(apply(df_clr_Microbiome_Species_He,1,function(x)(x-min(x)))))

df_Microbiome_Genus_He <- read.table(paste0(directory_name,"He_Controls_Genus_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_He <- as.data.frame(as.matrix(clr(df_Microbiome_Genus_He+0.00001)))
df_clr_Microbiome_Genus_He <- as.data.frame(t(apply(df_clr_Microbiome_Genus_He,1,function(x)(x-min(x)))))

## Reading LogMPie
df_Microbiome_Species_logmpie <- read.table(paste0(directory_name,"LogMPie_Species_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_logmpie <- as.data.frame(as.matrix(clr(df_Microbiome_Species_logmpie+0.00001)))
df_clr_Microbiome_Species_logmpie <- as.data.frame(t(apply(df_clr_Microbiome_Species_logmpie,1,function(x)(x-min(x)))))

df_Microbiome_Genus_logmpie <- read.table(paste0(directory_name,"LogMPie_Genus_Microbiome.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_logmpie <- as.data.frame(as.matrix(clr(df_Microbiome_Genus_logmpie+0.00001)))
df_clr_Microbiome_Genus_logmpie <- as.data.frame(t(apply(df_clr_Microbiome_Genus_logmpie,1,function(x)(x-min(x)))))

## Reading Dhakan DB
print("Reading DhakanDB")
load(paste0(directory_name,"DhakanDB_Microbiome.RData"))
df_clr_Microbiome_Species_DhakanDB <- as.data.frame(as.matrix(clr(DhakanDB_2019_Species+0.00001)))
df_clr_Microbiome_Species_DhakanDB <- as.data.frame(t(apply(df_clr_Microbiome_Species_DhakanDB,1,function(x)(x-min(x)))))
df_clr_Microbiome_Species_DhakanDB <- df_clr_Microbiome_Species_DhakanDB[,setdiff(colnames(df_clr_Microbiome_Species_DhakanDB),"study_name")]

df_clr_Microbiome_Genus_DhakanDB <- as.data.frame(as.matrix(clr(DhakanDB_2019_Genus+0.00001)))
df_clr_Microbiome_Genus_DhakanDB <- as.data.frame(t(apply(df_clr_Microbiome_Genus_DhakanDB,1,function(x)(x-min(x)))))
df_clr_Microbiome_Genus_DhakanDB <- df_clr_Microbiome_Genus_DhakanDB[,setdiff(colnames(df_clr_Microbiome_Genus_DhakanDB),"study_name")]

## Reading GuptaA
print("Reading GuptaA")
load(paste0(directory_name,"GuptaA_Microbiome.RData"))
df_clr_Microbiome_Species_GuptaA <- as.data.frame(as.matrix(clr(GuptaA_2019_Species+0.00001)))
df_clr_Microbiome_Species_GuptaA <- as.data.frame(t(apply(df_clr_Microbiome_Species_GuptaA,1,function(x)(x-min(x)))))
df_clr_Microbiome_Species_GuptaA <- df_clr_Microbiome_Species_GuptaA[,setdiff(colnames(df_clr_Microbiome_Species_GuptaA),"study_name")]

df_clr_Microbiome_Genus_GuptaA <- as.data.frame(as.matrix(clr(GuptaA_2019_Genus+0.00001)))
df_clr_Microbiome_Genus_GuptaA <- as.data.frame(t(apply(df_clr_Microbiome_Genus_GuptaA,1,function(x)(x-min(x)))))
df_clr_Microbiome_Genus_GuptaA <- df_clr_Microbiome_Genus_GuptaA[,setdiff(colnames(df_clr_Microbiome_Genus_GuptaA),"study_name")]

print("Reading MicroDiab India")
df_Microbiome_Species_MicroDiab <- read.table(paste0(directory_name,"MicroDiab_India.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Species_MicroDiab <- as.data.frame(as.matrix(clr(df_Microbiome_Species_MicroDiab+0.00001)))
df_clr_Microbiome_Species_MicroDiab <- as.data.frame(t(apply(df_clr_Microbiome_Species_MicroDiab,1,function(x)(x-min(x)))))

df_Microbiome_Genus_MicroDiab <- read.table(paste0(directory_name,"MicroDiab_India.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Microbiome_Genus_MicroDiab <- as.data.frame(as.matrix(clr(df_Microbiome_Genus_MicroDiab+0.00001)))
df_clr_Microbiome_Genus_MicroDiab <- as.data.frame(t(apply(df_clr_Microbiome_Genus_MicroDiab,1,function(x)(x-min(x)))))

## Microbiome Core Comparison
Threshold=0.80
logmpie_core_species <- names(which(apply(df_clr_Microbiome_Species_logmpie,2,function(x)(length(x[x>0])))/nrow(df_clr_Microbiome_Species_logmpie)>=Threshold))
AIIMS_1_Delhi_Controls_core_species <- names(which(apply(df_clr_Microbiome_Species_AIIMS_1_IBD[AIIMS_1_IBD_Controls,],2,function(x)(length(x[x>0])))/length(AIIMS_1_IBD_Controls)>=Threshold))
AIIMS_1_Leh_Controls_core_species <- names(which(apply(df_clr_Microbiome_Species_AIIMS_1_IBD[AIIMS_1_IBD_Leh_Controls,],2,function(x)(length(x[x>0])))/length(AIIMS_1_IBD_Leh_Controls)>=Threshold))
AIIMS_1_IBD_Cases_core_species <- names(which(apply(df_clr_Microbiome_Species_AIIMS_1_IBD[AIIMS_1_IBD_Cases,],2,function(x)(length(x[x>0])))/length(AIIMS_1_IBD_Cases)>=Threshold))
AIIMS_ITB_Controls_core_species <- names(which(apply(df_clr_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Controls,],2,function(x)(length(x[x>0])))/length(AIIMS_ITB_Controls)>=Threshold))
AIIMS_ITB_cases_core_species <- names(which(apply(df_clr_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Cases,],2,function(x)(length(x[x>0])))/length(AIIMS_ITB_Cases)>=Threshold))
AIIMS_CD_cases_core_species <- names(which(apply(df_clr_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases,],2,function(x)(length(x[x>0])))/length(AIIMS_ITB_CD_Cases)>=Threshold))
India_Core_Reference <- intersect(logmpie_core_species,intersect(AIIMS_1_Delhi_Controls_core_species,AIIMS_1_Leh_Controls_core_species))
India_Healthy_Microbiome_Reference <- c(rownames(df_clr_Microbiome_Species_AIIMS_1_IBD[AIIMS_1_IBD_Controls,]),rownames(df_clr_Microbiome_Species_AIIMS_1_IBD[AIIMS_1_IBD_Leh_Controls,]),rownames(df_clr_Microbiome_Species_logmpie))


He_core_species <- names(which(apply(df_clr_Microbiome_Species_He,2,function(x)(length(x[x>0])))/nrow(df_clr_Microbiome_Species_He)>=Threshold))
Chinese_ITB_Cases_core_species <- names(which(apply(df_clr_Microbiome_Species_Chinese_ITB[Chinese_ITB_Cases,],2,function(x)(length(x[x>0])))/length(Chinese_ITB_Cases)>=Threshold))
Chinese_ITB_Controls_core_species <- names(which(apply(df_clr_Microbiome_Species_Chinese_ITB[Chinese_ITB_Controls,],2,function(x)(length(x[x>0])))/length(Chinese_ITB_Controls)>=Threshold))

AIIMS_ITB_All_Groups_Unique_Core <- unique(c(AIIMS_ITB_Controls_core_species,AIIMS_ITB_cases_core_species,AIIMS_CD_cases_core_species))

df_Core_Comparison <- data.frame(Control=ifelse(AIIMS_ITB_All_Groups_Unique_Core %in% AIIMS_ITB_Controls_core_species,1,0),ITB=ifelse(AIIMS_ITB_All_Groups_Unique_Core %in% AIIMS_ITB_cases_core_species,1,0),CD=ifelse(AIIMS_ITB_All_Groups_Unique_Core %in% AIIMS_CD_cases_core_species,1,0),India_Reference=ifelse(AIIMS_ITB_All_Groups_Unique_Core %in% India_Core_Reference,1,0),row.names=AIIMS_ITB_All_Groups_Unique_Core)

pdf(paste0(directory_name,"Core_Comparison.pdf"))
heatmap.2(t(df_Core_Comparison),density="none",trace="none",col=c("white","blue"),margins=c(15,15),lwid=c(1,5))
dev.off()
df_India_Core_Abundance <- data.frame(Abundance_Control_RefOverlap = rowMeans(apply(df_clr_Microbiome_Species_AIIMS_ITB[,intersect(AIIMS_ITB_Controls_core_species,India_Core_Reference)],2,rank_scale)), Abundance_ITB_CD_RefOverlap = rowMeans(apply(df_clr_Microbiome_Species_AIIMS_ITB[,intersect(unique(c(AIIMS_ITB_cases_core_species,AIIMS_CD_cases_core_species)),India_Core_Reference)],2,rank_scale)), Abundance_ITB_CD_RefDiff = rowMeans(apply(df_clr_Microbiome_Species_AIIMS_ITB[,setdiff(unique(c(AIIMS_ITB_cases_core_species,AIIMS_CD_cases_core_species)),India_Core_Reference)],2,rank_scale)), Abundance_ITB_Specific = rowMeans(apply(df_clr_Microbiome_Species_AIIMS_ITB[,setdiff(AIIMS_ITB_cases_core_species,unique(c(AIIMS_CD_cases_core_species,AIIMS_ITB_Controls_core_species,India_Core_Reference)))],2,rank_scale)), Group=factor(ifelse(rownames(df_clr_Microbiome_Species_AIIMS_ITB) %in% AIIMS_ITB_Controls,"Control",ifelse(rownames(df_clr_Microbiome_Species_AIIMS_ITB) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Control","ITB","CD")))

### Fungome Scanning Starts here
print("Fungome Scanning Starts")
print("PRJDB7616")
df_Fungome_Species_PRJDB7616 <- read.table(paste0(directory_name,"PRJDB7616.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_PRJDB7616 <- as.data.frame(clr(df_Fungome_Species_PRJDB7616+0.00001))
df_clr_Fungome_Species_PRJDB7616 <- as.data.frame(t(apply(df_clr_Fungome_Species_PRJDB7616,1,function(x)(x-min(x)))))

df_Fungome_Genus_PRJDB7616 <- read.table(paste0(directory_name,"PRJDB7616.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_PRJDB7616 <- as.data.frame(clr(df_Fungome_Genus_PRJDB7616+0.00001))
df_clr_Fungome_Genus_PRJDB7616 <- as.data.frame(t(apply(df_clr_Fungome_Genus_PRJDB7616,1,function(x)(x-min(x)))))

print("PRJEB42375")
df_Fungome_Species_PRJEB42375 <- read.table(paste0(directory_name,"PRJEB42375.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_PRJEB42375 <- as.data.frame(clr(df_Fungome_Species_PRJEB42375+0.00001))
df_clr_Fungome_Species_PRJEB42375 <- as.data.frame(t(apply(df_clr_Fungome_Species_PRJEB42375,1,function(x)(x-min(x)))))

df_Fungome_Genus_PRJEB42375 <- read.table(paste0(directory_name,"PRJEB42375.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_PRJEB42375 <- as.data.frame(clr(df_Fungome_Genus_PRJEB42375+0.00001))
df_clr_Fungome_Genus_PRJEB42375 <- as.data.frame(t(apply(df_clr_Fungome_Genus_PRJEB42375,1,function(x)(x-min(x)))))

##Reading Liguori_et_al Fungome data
print("Liguori_et_al_Fungome_data")
df_Fungome_Species_Liguori <- read.table(paste0(directory_name,"Liguori_ITS_Species_Table.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_Liguori <- as.data.frame(clr(df_Fungome_Species_Liguori+0.00001))
df_clr_Fungome_Species_Liguori <- as.data.frame(t(apply(df_clr_Fungome_Species_Liguori,1,function(x)(x-min(x)))))

df_Fungome_Genus_Liguori <- read.table(paste0(directory_name,"Liguori_ITS_Genus_Table.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_Liguori <- as.data.frame(clr(df_Fungome_Genus_Liguori+0.00001))
df_clr_Fungome_Genus_Liguori <- as.data.frame(t(apply(df_clr_Fungome_Genus_Liguori,1,function(x)(x-min(x)))))

print("PRJNA356769")
df_Fungome_Species_PRJNA356769 <- read.table(paste0(directory_name,"PRJNA356769.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_PRJNA356769 <- as.data.frame(clr(df_Fungome_Species_PRJNA356769+0.00001))
df_clr_Fungome_Species_PRJNA356769 <- as.data.frame(t(apply(df_clr_Fungome_Species_PRJNA356769,1,function(x)(x-min(x)))))

df_Fungome_Genus_PRJNA356769 <- read.table(paste0(directory_name,"PRJNA356769.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_PRJNA356769 <- as.data.frame(clr(df_Fungome_Genus_PRJNA356769+0.00001))
df_clr_Fungome_Genus_PRJNA356769 <- as.data.frame(t(apply(df_clr_Fungome_Genus_PRJNA356769,1,function(x)(x-min(x)))))

print("PRJNA662173")
df_Fungome_Species_PRJNA662173 <- read.table(paste0(directory_name,"PRJNA662173.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_PRJNA662173 <- as.data.frame(clr(df_Fungome_Species_PRJNA662173+0.00001))
df_clr_Fungome_Species_PRJNA662173 <- as.data.frame(t(apply(df_clr_Fungome_Species_PRJNA662173,1,function(x)(x-min(x)))))

df_Fungome_Genus_PRJNA662173 <- read.table(paste0(directory_name,"PRJNA662173.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_PRJNA662173 <- as.data.frame(clr(df_Fungome_Genus_PRJNA662173+0.00001))
df_clr_Fungome_Genus_PRJNA662173 <- as.data.frame(t(apply(df_clr_Fungome_Genus_PRJNA662173,1,function(x)(x-min(x)))))

print("AIIMS ITB")
df_Fungome_Species_AIIMS_ITB <- read.table(paste0(directory_name,"AIIMS_ITB_Myco.species.txt"),sep="\t",row.names=1,header=TRUE)
rownames(df_Fungome_Species_AIIMS_ITB) <- sub("sample-","itb_sample_",rownames(df_Fungome_Species_AIIMS_ITB))
df_clr_Fungome_Species_AIIMS_ITB <- as.data.frame(clr(df_Fungome_Species_AIIMS_ITB+0.00001))
df_clr_Fungome_Species_AIIMS_ITB <- as.data.frame(t(apply(df_clr_Fungome_Species_AIIMS_ITB,1,function(x)(x-min(x)))))

df_Fungome_Genus_AIIMS_ITB <- read.table(paste0(directory_name,"AIIMS_ITB_Myco.genus.txt"),sep="\t",row.names=1,header=TRUE)
rownames(df_Fungome_Genus_AIIMS_ITB) <- sub("sample-","itb_sample_",rownames(df_Fungome_Genus_AIIMS_ITB))
df_clr_Fungome_Genus_AIIMS_ITB <- as.data.frame(clr(df_Fungome_Genus_AIIMS_ITB+0.00001))
df_clr_Fungome_Genus_AIIMS_ITB <- as.data.frame(t(apply(df_clr_Fungome_Genus_AIIMS_ITB,1,function(x)(x-min(x)))))

AIIMS_ITB_Controls_Fungome <- intersect(rownames(df_clr_Fungome_Species_AIIMS_ITB),AIIMS_ITB_Controls)
AIIMS_ITB_Cases_Fungome <- intersect(rownames(df_clr_Fungome_Species_AIIMS_ITB),AIIMS_ITB_Cases)
AIIMS_ITB_CD_Cases_Fungome <- intersect(rownames(df_clr_Fungome_Species_AIIMS_ITB),AIIMS_ITB_CD_Cases)

df_clr_Fungome_Species_AIIMS_ITB <- df_clr_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome),]
df_clr_Fungome_Genus_AIIMS_ITB <- df_clr_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome),]

print("PRJNA647266")
df_Fungome_Species_PRJNA647266 <- read.table(paste0(directory_name,"PRJNA647266.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_PRJNA647266 <- as.data.frame(clr(df_Fungome_Species_PRJNA647266+0.00001))
df_clr_Fungome_Species_PRJNA647266 <- as.data.frame(t(apply(df_clr_Fungome_Species_PRJNA647266,1,function(x)(x-min(x)))))

df_Fungome_Genus_PRJNA647266 <- read.table(paste0(directory_name,"PRJNA647266.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_PRJNA647266 <- as.data.frame(clr(df_Fungome_Genus_PRJNA647266+0.00001))
df_clr_Fungome_Genus_PRJNA647266 <- as.data.frame(t(apply(df_clr_Fungome_Genus_PRJNA647266,1,function(x)(x-min(x)))))

print("PRJNA439151")
df_Fungome_Species_PRJNA439151 <- read.table(paste0(directory_name,"PRJNA439151.species.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Species_PRJNA439151 <- as.data.frame(clr(df_Fungome_Species_PRJNA439151+0.00001))
df_clr_Fungome_Species_PRJNA439151 <- as.data.frame(t(apply(df_clr_Fungome_Species_PRJNA439151,1,function(x)(x-min(x)))))

df_Fungome_Genus_PRJNA439151 <- read.table(paste0(directory_name,"PRJNA439151.genus.txt"),sep="\t",row.names=1,header=TRUE)
df_clr_Fungome_Genus_PRJNA439151 <- as.data.frame(clr(df_Fungome_Genus_PRJNA439151+0.00001))
df_clr_Fungome_Genus_PRJNA439151 <- as.data.frame(t(apply(df_clr_Fungome_Genus_PRJNA439151,1,function(x)(x-min(x)))))

print("Remove Sparse Taxa")
Selected_Fungal_Species <- which(colSums(apply(df_clr_Fungome_Species_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))>=4)
df_clr_Fungome_Species_AIIMS_ITB <- df_clr_Fungome_Species_AIIMS_ITB[,Selected_Fungal_Species]

Selected_Fungal_Genus <- which(colSums(apply(df_clr_Fungome_Genus_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))>=4)
df_clr_Fungome_Genus_AIIMS_ITB <- df_clr_Fungome_Genus_AIIMS_ITB[,Selected_Fungal_Genus]

Selected_Microbiome_Species <- which(colSums(apply(df_clr_Microbiome_Species_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))>=4)
df_clr_Microbiome_Species_AIIMS_ITB <- df_clr_Microbiome_Species_AIIMS_ITB[,Selected_Microbiome_Species]

Selected_Microbiome_Genus <- which(colSums(apply(df_clr_Microbiome_Genus_AIIMS_ITB,2,function(x)(ifelse(x>0,1,0))))>=4)
df_clr_Microbiome_Genus_AIIMS_ITB <- df_clr_Microbiome_Genus_AIIMS_ITB[,Selected_Microbiome_Genus]

df_relab_Microbiome_Species_AIIMS_ITB <- df_Microbiome_Species_AIIMS_ITB/rowSums(df_Microbiome_Species_AIIMS_ITB)
df_relab_Microbiome_Genus_AIIMS_ITB <- df_Microbiome_Genus_AIIMS_ITB/rowSums(df_Microbiome_Genus_AIIMS_ITB)

df_relab_Fungome_Species_AIIMS_ITB <- df_Fungome_Species_AIIMS_ITB/rowSums(df_Fungome_Species_AIIMS_ITB)
df_relab_Fungome_Genus_AIIMS_ITB <- df_Fungome_Genus_AIIMS_ITB/rowSums(df_Fungome_Genus_AIIMS_ITB)

Fungome_Filtered_Rows <- rownames(df_relab_Fungome_Species_AIIMS_ITB)[!is.nan(rowSums(df_relab_Fungome_Species_AIIMS_ITB))]

df_relab_Fungome_Species_AIIMS_ITB <- df_relab_Fungome_Species_AIIMS_ITB[Fungome_Filtered_Rows,]
df_relab_Fungome_Genus_AIIMS_ITB <- df_relab_Fungome_Genus_AIIMS_ITB[Fungome_Filtered_Rows,]
df_relab_ASV_AIIMS_ITB <- t(df_ASV_AIIMS_ITB)/rowSums(t(df_ASV_AIIMS_ITB))

AlphaDiversity_Microbiome <- data.frame(shannon_species = diversity(df_relab_Microbiome_Species_AIIMS_ITB),shannon_genus = diversity(df_relab_Microbiome_Genus_AIIMS_ITB),pielou_species = diversity(df_relab_Microbiome_Species_AIIMS_ITB)/log(specnumber(df_relab_Microbiome_Species_AIIMS_ITB)),pielou_genus = diversity(df_relab_Microbiome_Genus_AIIMS_ITB)/log(specnumber(df_relab_Microbiome_Genus_AIIMS_ITB)),diversity_asv = diversity(df_relab_ASV_AIIMS_ITB), pielou_asv = diversity(df_relab_ASV_AIIMS_ITB)/log(specnumber(df_relab_ASV_AIIMS_ITB)))
AlphaDiversity_Microbiome[AIIMS_ITB_Controls,"Group"] <- "Controls"
AlphaDiversity_Microbiome[AIIMS_ITB_Cases,"Group"] <- "ITB"
AlphaDiversity_Microbiome[AIIMS_ITB_CD_Cases,"Group"] <- "CD"

df_OTU_Fungome_AIIMS_ITB <- df_OTU_Fungome_AIIMS_ITB <- read.table(paste0(directory_name,"Otu_table_fungome_AIIMS_ITB.txt"),sep="\t",row.names=1,header=TRUE,check.names=FALSE)
colnames(df_OTU_Fungome_AIIMS_ITB) = paste0('itb_',sub('-','_',colnames(df_OTU_Fungome_AIIMS_ITB)))
df_OTU_Fungome_AIIMS_ITB = t(df_OTU_Fungome_AIIMS_ITB)
df_relab_ITS_AIIMS_ITB <- df_OTU_Fungome_AIIMS_ITB/rowSums(df_OTU_Fungome_AIIMS_ITB)
df_clr_Fungome_AIIMS_ITB_otu <- as.data.frame(as.matrix(clr(df_OTU_Fungome_AIIMS_ITB+0.00001)))
df_clr_Fungome_AIIMS_ITB_otu <- as.data.frame(t(apply(df_clr_Fungome_AIIMS_ITB_otu,1,function(x)(x-min(x)))))
df_relab_ITS_AIIMS_ITB <- df_relab_ITS_AIIMS_ITB[Fungome_Filtered_Rows,]

AlphaDiversity_Fungome <- data.frame(shannon_species = diversity(df_relab_Fungome_Species_AIIMS_ITB),shannon_genus = diversity(df_relab_Fungome_Genus_AIIMS_ITB),pielou_species = diversity(df_relab_Fungome_Species_AIIMS_ITB)/log(specnumber(df_relab_Fungome_Species_AIIMS_ITB)),pielou_genus = diversity(df_relab_Fungome_Genus_AIIMS_ITB)/log(specnumber(df_relab_Fungome_Genus_AIIMS_ITB)),diversity_its = diversity(df_relab_ITS_AIIMS_ITB), pielou_asv = diversity(df_relab_ITS_AIIMS_ITB)/log(specnumber(df_relab_ITS_AIIMS_ITB)))
AlphaDiversity_Fungome <- AlphaDiversity_Fungome[intersect(Fungome_Filtered_Rows,c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome)),] 
AlphaDiversity_Fungome[intersect(rownames(AlphaDiversity_Fungome),AIIMS_ITB_Controls_Fungome),"Group"] <- "Controls"
AlphaDiversity_Fungome[intersect(rownames(AlphaDiversity_Fungome),AIIMS_ITB_Cases_Fungome),"Group"] <- "ITB"
AlphaDiversity_Fungome[intersect(rownames(AlphaDiversity_Fungome),AIIMS_ITB_CD_Cases_Fungome),"Group"] <- "CD"

save.image(paste0(directory_name,"Stage1_Workflow_Output.RData"))
