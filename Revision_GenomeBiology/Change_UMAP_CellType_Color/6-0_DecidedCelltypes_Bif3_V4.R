## Make metafile & plot for estimated annotation
library(ggplot2)
library(stringr)
source("/home/sb14489/Epigenomics/scATAC-seq/0_Function/DrawFigures_QC_Annotation_forUMAP.R")

## Make new metafile bycombining all the meta
setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3")
meta <- "bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt"
loaded_meta_data <- read.table(meta)
PreAnn <- loaded_meta_data

CellOrder <- readLines("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt")
Cluster1 <- read.table("bif3_Cluster1_Recluster_Sub_res1_knear100_Partmetadata.txt")
head(Cluster1)


head(loaded_meta_data)
#levels(factor(loaded_meta_data$LouvainClusters))
loaded_meta_data[rownames(Cluster1),]$LouvainClusters <- Cluster1$LouvainClusters
levels(factor(loaded_meta_data$LouvainClusters))

NewMeta <- loaded_meta_data
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="1"),]

dim(loaded_meta_data)
dim(NewMeta)

levels(factor(NewMeta$LouvainClusters))
NewMeta$Ann_v4 <- "Ann_v4"
levels(factor(NewMeta$Ann_v4))
## New assignment ##
NewMeta[which(NewMeta$LouvainClusters =="1_1"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="1_2"),]$Ann_v4 <- "SPM-base_SM-base"
NewMeta[which(NewMeta$LouvainClusters =="1_3"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="2"),]$Ann_v4 <- "Unknown1"
NewMeta[which(NewMeta$LouvainClusters =="3"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="4"),]$Ann_v4 <- "XylemParenchyma_PithParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="5"),]$Ann_v4 <- "ProcambialMeristem_ProtoXylem_MetaXylem"
NewMeta[which(NewMeta$LouvainClusters =="6"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="7"),]$Ann_v4 <- "Unknown2"
NewMeta[which(NewMeta$LouvainClusters =="8"),]$Ann_v4 <- "L1atFloralMeristem"
NewMeta[which(NewMeta$LouvainClusters =="9"),]$Ann_v4 <- "L1"
NewMeta[which(NewMeta$LouvainClusters =="10"),]$Ann_v4 <- "Unknown_lowFRiP"
NewMeta[which(NewMeta$LouvainClusters =="11"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="12"),]$Ann_v4 <- "G2_M"
NewMeta[which(NewMeta$LouvainClusters =="13"),]$Ann_v4 <- "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="14"),]$Ann_v4 <- "IM-OC"
######################

setwd("/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3")

write.table(NewMeta, file=paste0("Bif3_AnnV4_metadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")

DrawUMAP_Ann_QC(PreAnn,NewMeta, "Ann_v4", CellOrder, "bi3_Re1", "bi3_Re2","Bif3_AnnV4")

