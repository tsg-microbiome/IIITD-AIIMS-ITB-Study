print("Stage 2: Generating PCoA and Distance Profiles")
#dir <- "E:/Data/AIIMS"
print("Generate for Jaccard")
Dist_Jaccard_Fungome_Species_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Fungome_Species_AIIMS_ITB,method="jaccard"))
Dist_Jaccard_Fungome_Species_AIIMS_ITB <- apply(Dist_Jaccard_Fungome_Species_AIIMS_ITB,2,function(x)(ifelse(is.nan(x),1,x)))
diag(Dist_Jaccard_Fungome_Species_AIIMS_ITB) <- NA
Dist_Jaccard_Fungome_Genus_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Fungome_Genus_AIIMS_ITB,method="jaccard"))
Dist_Jaccard_Fungome_Genus_AIIMS_ITB <- apply(Dist_Jaccard_Fungome_Genus_AIIMS_ITB,2,function(x)(ifelse(is.nan(x),1,x)))
diag(Dist_Jaccard_Fungome_Genus_AIIMS_ITB) <- NA
Dist_Jaccard_Microbiome_Species_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Microbiome_Species_AIIMS_ITB,method="jaccard"))
diag(Dist_Jaccard_Microbiome_Species_AIIMS_ITB) <- NA
Dist_Jaccard_Microbiome_Genus_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Microbiome_Genus_AIIMS_ITB,method="jaccard"))
diag(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB) <- NA

DistanceProfiles_Jaccard <- data.frame(ControlCentroid_Jaccard_MicroGenus=apply(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Jaccard_MicroSpecies=apply(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Jaccard_FungGenus=apply(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Jaccard_FungSpecies=apply(Dist_Jaccard_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),WithinCohort_Jaccard_MicroGenus=c(apply(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Jaccard_MicroSpecies=c(apply(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Jaccard_FungGenus=c(apply(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Jaccard_FungSpecies=c(apply(Dist_Jaccard_Fungome_Species_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Fungome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Jaccard_Fungome_Species_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),row.names=c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome))
DistanceProfiles_Jaccard$Group <- NA
DistanceProfiles_Jaccard$Group <- factor(ifelse(rownames(DistanceProfiles_Jaccard) %in% AIIMS_ITB_Controls_Fungome,"Controls",ifelse(rownames(DistanceProfiles_Jaccard) %in% AIIMS_ITB_Cases_Fungome,"ITB","CD")),levels=c("Controls","CD","ITB"))

pco_AIIMS_ITB_Jaccard_Fungome_Species <- dudi.pco(vegdist(df_clr_Fungome_Species_AIIMS_ITB[rowSums(df_clr_Fungome_Species_AIIMS_ITB)>0,],method="jaccard"),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Jaccard_Fungome_Species_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Jaccard_Fungome_Species$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Jaccard_Fungome_Species$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)

pco_AIIMS_ITB_Jaccard_Fungome_Genus <- dudi.pco(vegdist(df_clr_Fungome_Genus_AIIMS_ITB[rowSums(df_clr_Fungome_Genus_AIIMS_ITB)>0,],method="jaccard"),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Jaccard_Fungome_Genus_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)


pco_AIIMS_ITB_Jaccard_Microbiome_Genus <- dudi.pco(vegdist(df_clr_Microbiome_Genus_AIIMS_ITB[rowSums(df_clr_Microbiome_Genus_AIIMS_ITB)>0,],method="jaccard"),scannf=FALSE,nf=20)
s.class(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()

pco_AIIMS_ITB_Jaccard_Microbiome_Species <- dudi.pco(vegdist(df_clr_Microbiome_Species_AIIMS_ITB[rowSums(df_clr_Microbiome_Species_AIIMS_ITB)>0,],method="jaccard"),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Jaccard_Microbiome_Species_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()

#######

print("Generate for Aitchison")
Dist_Aitchison_Fungome_Species_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Fungome_Species_AIIMS_ITB,method="euclidean"))
diag(Dist_Aitchison_Fungome_Species_AIIMS_ITB) <- NA
Dist_Aitchison_Fungome_Genus_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Fungome_Genus_AIIMS_ITB,method="euclidean"))
diag(Dist_Aitchison_Fungome_Genus_AIIMS_ITB) <- NA
Dist_Aitchison_Microbiome_Species_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Microbiome_Species_AIIMS_ITB,method="euclidean"))
diag(Dist_Aitchison_Microbiome_Species_AIIMS_ITB) <- NA
Dist_Aitchison_Microbiome_Genus_AIIMS_ITB <-  as.matrix(vegdist(df_clr_Microbiome_Genus_AIIMS_ITB,method="euclidean"))
diag(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB) <- NA

DistanceProfiles_Aitchison <- data.frame(ControlCentroid_Aitchison_MicroGenus=apply(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Aitchison_MicroSpecies=apply(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Aitchison_FungGenus=apply(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Aitchison_FungSpecies=apply(Dist_Aitchison_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),WithinCohort_Aitchison_MicroGenus=c(apply(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Aitchison_MicroSpecies=c(apply(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Aitchison_FungGenus=c(apply(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Aitchison_FungSpecies=c(apply(Dist_Aitchison_Fungome_Species_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Fungome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Aitchison_Fungome_Species_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),row.names=c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome))
DistanceProfiles_Aitchison$Group <- NA
DistanceProfiles_Aitchison$Group <- factor(ifelse(rownames(DistanceProfiles_Aitchison) %in% AIIMS_ITB_Controls_Fungome,"Controls",ifelse(rownames(DistanceProfiles_Aitchison) %in% AIIMS_ITB_Cases_Fungome,"ITB","CD")),levels=c("Controls","CD","ITB"))


pco_AIIMS_ITB_Aitchison_Fungome_Species <- dudi.pco(vegdist(df_clr_Fungome_Species_AIIMS_ITB[rowSums(df_clr_Fungome_Species_AIIMS_ITB)>0,],method="euclidean"),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Aitchison_Fungome_Species_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Aitchison_Fungome_Species$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Aitchison_Fungome_Species$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)

pco_AIIMS_ITB_Aitchison_Fungome_Genus <- dudi.pco(vegdist(df_clr_Fungome_Genus_AIIMS_ITB[rowSums(df_clr_Fungome_Genus_AIIMS_ITB)>0,],method="euclidean"),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Aitchison_Fungome_Genus_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)

pco_AIIMS_ITB_Aitchison_Microbiome_Genus <- dudi.pco(vegdist(df_clr_Microbiome_Genus_AIIMS_ITB[rowSums(df_clr_Microbiome_Genus_AIIMS_ITB)>0,],method="euclidean"),scannf=FALSE,nf=20)
s.class(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()

pco_AIIMS_ITB_Aitchison_Microbiome_Species <- dudi.pco(vegdist(df_clr_Microbiome_Species_AIIMS_ITB[rowSums(df_clr_Microbiome_Species_AIIMS_ITB)>0,],method="euclidean"),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Aitchison_Microbiome_Species_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)

print("Generate for Kendall")

Dist_Kendall_Fungome_Species_AIIMS_ITB <-  1-cor.fk(t(df_clr_Fungome_Species_AIIMS_ITB))/2
diag(Dist_Kendall_Fungome_Species_AIIMS_ITB) <- NA
Dist_Kendall_Fungome_Genus_AIIMS_ITB <-  1-cor.fk(t(df_clr_Fungome_Genus_AIIMS_ITB))/2
diag(Dist_Kendall_Fungome_Genus_AIIMS_ITB) <- NA
Dist_Kendall_Microbiome_Species_AIIMS_ITB <-  1-cor.fk(t(df_clr_Microbiome_Species_AIIMS_ITB))/2
diag(Dist_Kendall_Microbiome_Species_AIIMS_ITB) <- NA
Dist_Kendall_Microbiome_Genus_AIIMS_ITB <-  1-cor.fk(t(df_clr_Microbiome_Genus_AIIMS_ITB))/2
diag(Dist_Kendall_Microbiome_Genus_AIIMS_ITB) <- NA

DistanceProfiles_Kendall <- data.frame(ControlCentroid_Kendall_MicroGenus=apply(Dist_Kendall_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Kendall_MicroSpecies=apply(Dist_Kendall_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Kendall_FungGenus=apply(Dist_Kendall_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),ControlCentroid_Kendall_FungSpecies=apply(Dist_Kendall_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome),AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),WithinCohort_Kendall_MicroGenus=c(apply(Dist_Kendall_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Microbiome_Genus_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Kendall_MicroSpecies=c(apply(Dist_Kendall_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Microbiome_Species_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Kendall_FungGenus=c(apply(Dist_Kendall_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Fungome_Genus_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),WithinCohort_Kendall_FungSpecies=c(apply(Dist_Kendall_Fungome_Species_AIIMS_ITB[AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Controls_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Fungome_Species_AIIMS_ITB[AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome],1,function(x)(median(x[!is.na(x)]))),apply(Dist_Kendall_Fungome_Species_AIIMS_ITB[AIIMS_ITB_Cases_Fungome,AIIMS_ITB_Cases_Fungome],1,function(x)(median(x[!is.na(x)])))),row.names=c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome))
DistanceProfiles_Kendall$Group <- NA
DistanceProfiles_Kendall$Group <- factor(ifelse(rownames(DistanceProfiles_Kendall) %in% AIIMS_ITB_Controls_Fungome,"Controls",ifelse(rownames(DistanceProfiles_Kendall) %in% AIIMS_ITB_Cases_Fungome,"ITB","CD")),levels=c("Controls","CD","ITB"))

pco_AIIMS_ITB_Kendall_Fungome_Species <- dudi.pco(as.dist(1-cor(t(df_clr_Fungome_Species_AIIMS_ITB[rowSums(df_clr_Fungome_Species_AIIMS_ITB)>0,]),method="kendall")/2),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Kendall_Fungome_Species_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Kendall_Fungome_Species$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Kendall_Fungome_Species$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()

pco_AIIMS_ITB_Kendall_Fungome_Genus <- dudi.pco(as.dist(1-cor(t(df_clr_Fungome_Genus_AIIMS_ITB[rowSums(df_clr_Fungome_Genus_AIIMS_ITB)>0,]),method="kendall")/2),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Kendall_Fungome_Genus_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Kendall_Fungome_Genus$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Kendall_Fungome_Genus$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()

pco_AIIMS_ITB_Kendall_Microbiome_Species <- dudi.pco(as.dist(1-cor(t(df_clr_Microbiome_Species_AIIMS_ITB[rowSums(df_clr_Microbiome_Species_AIIMS_ITB)>0,]),method="kendall")/2),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Kendall_Microbiome_Species_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Kendall_Microbiome_Species$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Kendall_Microbiome_Species$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()

pco_AIIMS_ITB_Kendall_Microbiome_Genus <- dudi.pco(as.dist(1-cor(t(df_clr_Microbiome_Genus_AIIMS_ITB[rowSums(df_clr_Microbiome_Genus_AIIMS_ITB)>0,]),method="kendall")/2),scannf=FALSE,nf=20)
#pdf(paste0(dir,"\\PCo_Kendall_Microbiome_Genus_AIIMS_ITB.pdf"))
s.class(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li[,c(1:2)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
s.class(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li[,c(1,3)],factor(ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li) %in% AIIMS_ITB_Controls,"Controls",ifelse(rownames(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li) %in% AIIMS_ITB_Cases,"ITB","CD")),levels=c("Controls","CD","ITB")),col=c("Blue","Darkgoldenrod2","Pink4"),plabels.col="black",plabels.cex=1.5,plabels.boxes.alpha=0.7)
#dev.off()




#df_clr_Fungome_Species_AIIMS_ITB <- df_clr_Fungome_Species_AIIMS_ITB[rowSums(df_clr_Fungome_Species_AIIMS_ITB)>0,]
#df_clr_Fungome_Genus_AIIMS_ITB <- df_clr_Fungome_Genus_AIIMS_ITB[rowSums(df_clr_Fungome_Genus_AIIMS_ITB)>0,]

source(paste0(directory_name,"batch_dunns_function.R"))

Group_Status_Microbiome <- data.frame(Group=c(rep("Controls",length(AIIMS_ITB_Controls)),rep("CD",length(AIIMS_ITB_CD_Cases)),rep("ITB",length(AIIMS_ITB_Cases))),row.names=c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases,AIIMS_ITB_Cases))

Group_Status_Fungome <- data.frame(Group=c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("CD",length(AIIMS_ITB_CD_Cases_Fungome)),rep("ITB",length(AIIMS_ITB_Cases_Fungome))),row.names=c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome,AIIMS_ITB_Cases_Fungome))


batch_dunns_Jaccard_Microbiome_Species <- batch_dunns(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li,Group_Status_Microbiome[rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Species$li),])
batch_dunns_Aitchison_Microbiome_Species <- batch_dunns(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li,Group_Status_Microbiome[rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Species$li),])
batch_dunns_Kendall_Microbiome_Species <- batch_dunns(pco_AIIMS_ITB_Kendall_Microbiome_Species$li,Group_Status_Microbiome[rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li),])
batch_dunns_Jaccard_Microbiome_Genus <- batch_dunns(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li,Group_Status_Microbiome[rownames(pco_AIIMS_ITB_Jaccard_Microbiome_Genus$li),])
batch_dunns_Aitchison_Microbiome_Genus <- batch_dunns(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li,Group_Status_Microbiome[rownames(pco_AIIMS_ITB_Aitchison_Microbiome_Genus$li),])
batch_dunns_Kendall_Microbiome_Genus <- batch_dunns(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li,Group_Status_Microbiome[rownames(pco_AIIMS_ITB_Kendall_Microbiome_Genus$li),])
batch_dunns_Jaccard_Fungome_Species <- batch_dunns(pco_AIIMS_ITB_Jaccard_Fungome_Species$li,Group_Status_Fungome[rownames(pco_AIIMS_ITB_Jaccard_Fungome_Species$li),])
batch_dunns_Aitchison_Fungome_Species <- batch_dunns(pco_AIIMS_ITB_Aitchison_Fungome_Species$li,Group_Status_Fungome[rownames(pco_AIIMS_ITB_Aitchison_Fungome_Species$li),])
batch_dunns_Kendall_Fungome_Species <- batch_dunns(pco_AIIMS_ITB_Kendall_Fungome_Species$li,Group_Status_Fungome[rownames(pco_AIIMS_ITB_Kendall_Fungome_Species$li),])
batch_dunns_Jaccard_Fungome_Genus <- batch_dunns(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li,Group_Status_Fungome[rownames(pco_AIIMS_ITB_Jaccard_Fungome_Genus$li),])
batch_dunns_Aitchison_Fungome_Genus <- batch_dunns(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li,Group_Status_Fungome[rownames(pco_AIIMS_ITB_Aitchison_Fungome_Genus$li),])
batch_dunns_Kendall_Fungome_Genus <- batch_dunns(pco_AIIMS_ITB_Kendall_Fungome_Genus$li,Group_Status_Fungome[rownames(pco_AIIMS_ITB_Kendall_Fungome_Genus$li),])

mat <- as.matrix(batch_dunns_Kendall_Microbiome_Species$final_trends[4:20,])
heatmap.2(mat,density="none",trace="none",Rowv=FALSE,Colv=FALSE,col=c("Blue","White","Red"),sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat),rowsep=0:nrow(mat))

mat <- as.matrix(batch_dunns_Jaccard_Microbiome_Species$final_trends[4:20,])
heatmap.2(mat,density="none",trace="none",Rowv=FALSE,Colv=FALSE,col=c("Blue","White","Red"),sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat),rowsep=0:nrow(mat))

mat <- as.matrix(batch_dunns_Aitchison_Microbiome_Species$final_trends[4:20,])
heatmap.2(mat,density="none",trace="none",Rowv=FALSE,Colv=FALSE,col=c("Blue","White","Red"),sepwidth=c(0.1,0.1),sepcolor="black",colsep=0:ncol(mat),rowsep=0:nrow(mat))

df_PCoA_5 <- data.frame(Jaccard_PCo5 = pco_AIIMS_ITB_Jaccard_Microbiome_Species$li[,5],Aitchison_PCo5 = pco_AIIMS_ITB_Aitchison_Microbiome_Species$li[,5],Kendall_PCo5 = pco_AIIMS_ITB_Kendall_Microbiome_Species$li[,5],row.names=rownames(pco_AIIMS_ITB_Kendall_Microbiome_Species$li))

print("Running PERMANOVA")

#adonis(as.dist(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_Cases_Fungome))))

#adonis(as.dist(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Jaccard_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Case1",length(AIIMS_ITB_Cases_Fungome)),rep("Case2",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Jaccard_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_Cases_Fungome))))

#adonis(as.dist(Dist_Jaccard_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Jaccard_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Case1",length(AIIMS_ITB_Cases_Fungome)),rep("Case2",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_Cases))))

#adonis(as.dist(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Jaccard_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases)])~c(rep("Case1",length(AIIMS_ITB_Cases)),rep("Case2",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_Cases))))

#adonis(as.dist(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Jaccard_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases)])~c(rep("Case1",length(AIIMS_ITB_Cases)),rep("Case2",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_Cases_Fungome))))

#adonis(as.dist(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Aitchison_Fungome_Genus_AIIMS_ITB[c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Case1",length(AIIMS_ITB_Cases_Fungome)),rep("Case2",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Aitchison_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_Cases_Fungome))))

#adonis(as.dist(Dist_Aitchison_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Controls_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Controls",length(AIIMS_ITB_Controls_Fungome)),rep("Cases",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Aitchison_Fungome_Species_AIIMS_ITB[c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome),c(AIIMS_ITB_Cases_Fungome,AIIMS_ITB_CD_Cases_Fungome)])~c(rep("Case1",length(AIIMS_ITB_Cases_Fungome)),rep("Case2",length(AIIMS_ITB_CD_Cases_Fungome))))

#adonis(as.dist(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_Cases))))

#adonis(as.dist(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Aitchison_Microbiome_Species_AIIMS_ITB[c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases)])~c(rep("Case1",length(AIIMS_ITB_Cases)),rep("Case2",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_Cases))))

#adonis(as.dist(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Controls,AIIMS_ITB_CD_Cases)])~c(rep("Controls",length(AIIMS_ITB_Controls)),rep("Cases",length(AIIMS_ITB_CD_Cases))))

#adonis(as.dist(Dist_Aitchison_Microbiome_Genus_AIIMS_ITB[c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases),c(AIIMS_ITB_Cases,AIIMS_ITB_CD_Cases)])~c(rep("Case1",length(AIIMS_ITB_Cases)),rep("Case2",length(AIIMS_ITB_CD_Cases))))

save.image(paste0(directory_name,"Stage2_Workflow_Output.RData"))