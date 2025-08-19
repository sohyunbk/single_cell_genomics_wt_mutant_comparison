## MRake metafile & plot for estimated annotation
library(ggplot2)
library(stringr)
source("../../scATAC-seq/functions/DrawFigures_QC_Annotation_forUMAP.R")

## Make new metafile bycombining all the meta
meta <- "Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.metadata.txt"
CellOrder <- readLines("./Ann_v4_CellType_order_forA619Bif3.txt")
loaded_meta_data <- read.table(meta)
head(loaded_meta_data)
PreAnn <- loaded_meta_data
#Cluster7 <- read.table("Cluster7_Recluster_Sub_res1_knear100_Partmetadata.txt")
#head(Cluster7)
Cluster1 <- read.table("Cluster1_Recluster_Sub_res1_knear100_Partmetadata.txt")
Cluster3 <- read.table("Cluster3_Recluster_Sub_res1_knear100_Partmetadata.txt")
#Cluster4 <- read.table("Cluster4_Recluster_Sub_res1_knear100_Partmetadata.txt")

head(loaded_meta_data)
#levels(factor(loaded_meta_data$LouvainClusters))
loaded_meta_data[rownames(Cluster1),]$LouvainClusters <- Cluster1$LouvainClusters
loaded_meta_data[rownames(Cluster3),]$LouvainClusters <- Cluster3$LouvainClusters
#loaded_meta_data[rownames(Cluster4),]$LouvainClusters <- Cluster4$LouvainClusters
#loaded_meta_data[rownames(Cluster7),]$LouvainClusters <- Cluster7$LouvainClusters
levels(factor(loaded_meta_data$LouvainClusters))

NewMeta <- loaded_meta_data
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="1"),]
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="3"),]
#NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="4"),]
#NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="7"),]
dim(loaded_meta_data)
dim(NewMeta)
levels(factor(NewMeta$LouvainClusters))
NewMeta$Ann_v4 <- "Ann_v4"
levels(factor(NewMeta$Ann_v4))
NewMeta[which(NewMeta$LouvainClusters =="1_1"),]$Ann_v4 <- "Unknown_lowFRiP"
NewMeta[which(NewMeta$LouvainClusters =="1_2"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="1_3"),]$Ann_v4 <- "Unknown_Sclerenchyma"
NewMeta[which(NewMeta$LouvainClusters =="1_4"),]$Ann_v4 <- "Unknown_lowFRiP"
NewMeta[which(NewMeta$LouvainClusters =="2"),]$Ann_v4 <- "PhloemPrecursor"
NewMeta[which(NewMeta$LouvainClusters =="3_1"),]$Ann_v4 <- "IM-OC"
NewMeta[which(NewMeta$LouvainClusters =="3_3"),]$Ann_v4 <- "IM-OC"
NewMeta[which(NewMeta$LouvainClusters =="3_2"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="3_4"),]$Ann_v4 <- "FloralMeristem_SuppressedBract"
NewMeta[which(NewMeta$LouvainClusters =="4"),]$Ann_v4 <- "XylemParenchyma_PithParenchyma"
NewMeta[which(NewMeta$LouvainClusters =="5"),]$Ann_v4 <- "Unknown1"
NewMeta[which(NewMeta$LouvainClusters =="6"),]$Ann_v4 <- "L1atFloralMeristem"
NewMeta[which(NewMeta$LouvainClusters =="7"),]$Ann_v4 <- "ProcambialMeristem_ProtoXylem_MetaXylem"
NewMeta[which(NewMeta$LouvainClusters =="8"),]$Ann_v4 <- "Unknown2"
NewMeta[which(NewMeta$LouvainClusters =="9"),]$Ann_v4 <- "G2_M"
NewMeta[which(NewMeta$LouvainClusters =="10"),]$Ann_v4 <- "SPM-base_SM-base"
NewMeta[which(NewMeta$LouvainClusters =="11"),]$Ann_v4 <- "L1"
NewMeta[which(NewMeta$LouvainClusters =="12"),]$Ann_v4 <- "ProtoPhloem_MetaPhloem_CompanionCell_PhloemParenchyma"
length(which(NewMeta$LouvainClusters =="13"))
length(which(NewMeta$LouvainClusters =="13"))
NewMeta <- NewMeta[which(NewMeta$LouvainClusters !="13"),]


write.table(NewMeta, file=paste0("Ref_AnnV4_metadata.txt"),
            quote=F, row.names=T, col.names=T, sep="\t")

## Should be sampleID
DrawUMAP_Ann_QC(PreAnn,NewMeta, "Ann_v4", CellOrder, "A619_Re1", "A619_Re2","A619_AnnV4_SmallerDot")

