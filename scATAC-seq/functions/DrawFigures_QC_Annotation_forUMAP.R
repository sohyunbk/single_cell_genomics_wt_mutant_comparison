library(ggplot2)
library(stringr)
library("RColorBrewer")
library(gridExtra)


### Getting metafile as input.
DrawUMAP_Ann_QC <- function(PreAnnMeta,Meta, Slot, CellOrder, Re1, Re2,OutfilePathName){
  print("Replicates name should be sampleID slot")
  ColorForPreAnn <- c(
    "#944a0d", # orange
    "#56B4E9", # sky blue
    "#0736b5", # yellow
    "#aec8eb", # blue
    "#D55E00", # vermillion
    "#CC79A7", # reddish purple
    "#999999", # grey
    "#FB8072", # coral
    "#FDB462", # peach
    "#BC80BD", # lavender
    "#8DD3C7", # teal-blue
    "#7d22c7", # soft purple
    "#FCCDE5", # pink
    "#FFED6F")
  
PreAnn <- ggplot(PreAnnMeta, aes(x=umap1, y=umap2, color=factor(LouvainClusters))) +
    geom_point(size=0.001) +
    scale_color_manual(values=ColorForPreAnn)+theme_minimal()+
    guides(colour = guide_legend(override.aes = list(size=12),
                                 title="LouvainCluster_OtherQCplotis_are_BasedOnPreAnnotation"))+
    labs(title = paste0("Pre-Annotation\n CellNumber: ",nrow(PreAnnMeta)),
         x = "UMAP1",
         y = "UMAP2")+
    theme(legend.key.size = unit(1, "lines"),
          legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 20))

Meta$Ann <- Meta[[Slot]]
Meta$Ann <- factor(Meta$Ann,levels=(CellOrder))
colorr <- c("#4F96C4","#84f5d9","#d4ce1e","#d62744","#FDA33F","#060878","#a97df0",
            "#876b58","#800000", "#800075","#777d7d","#fc53b6","#DE9A89","#edc7f0",
            "#deadce","#fc53b6")

All <- ggplot(Meta, aes(x=umap1, y=umap2, color=factor(Ann))) +
  geom_point(size=0.01) +
  scale_color_manual(values=colorr)+
  theme_minimal()+
  guides(colour = guide_legend(override.aes = list(size=12)))+
  labs(title = paste0("Re1+R2 \n CellNumber: ",nrow(Meta)),
       x = "UMAP1",
       y = "UMAP2")+
  theme(legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

ClustersTable_Re1 <- subset(PreAnnMeta, sampleID == Re1)
Re1_plot <- ggplot(ClustersTable_Re1, aes(x=umap1, y=umap2, color=factor(Ann))) +
  geom_point(size=0.01, color="blue") +
  theme_minimal()+
  labs(title = paste0("Re1 : ",nrow(ClustersTable_Re1)))+
  theme(axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))
ClustersTable_Re2 <- subset(PreAnnMeta, sampleID == Re2)
Re2_plot <- ggplot(ClustersTable_Re2, aes(x=umap1, y=umap2, color=factor(Ann))) +
  geom_point(size=0.01, color="red") +
  theme_minimal()+
  labs(title = paste0("Re2 : ",nrow(ClustersTable_Re2)))+
  theme(axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

## Add some plots to see the cell quality
### * Tn5 log
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(PreAnnMeta$log10nSites),
                                      max(PreAnnMeta$log10nSites)))
Q_Tn5 <- ggplot(PreAnnMeta, aes(x=umap1, y=umap2,
                                color=log10nSites)) +
  geom_point(size=0.01) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "LogTn5")+sc+
  theme(legend.key.size = unit(4, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 20))
## * Doublets
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(PreAnnMeta$doubletscore),
                                      max(PreAnnMeta$doubletscore)))
Q_doubletscore <-ggplot(PreAnnMeta, aes(x=umap1, y=umap2,
                                        color=doubletscore)) +
  geom_point(size=0.01) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "Doublet score")+sc+
  theme(legend.key.size = unit(4, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 20))
## Tss ratio
PreAnnMeta$rTSS <- PreAnnMeta$tss/PreAnnMeta$total
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(PreAnnMeta$rTSS),
                                      max(PreAnnMeta$rTSS)))
Q_rTSS <-ggplot(PreAnnMeta, aes(x=umap1, y=umap2,
                                color=rTSS)) +
  geom_point(size=0.01) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "TSS ratio")+sc+
  theme(legend.key.size = unit(4, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 20))

## Tss ratio
sc <- scale_colour_gradientn(colours = myPalette(100),
                             limits=c(min(PreAnnMeta$FRiP),
                                      max(PreAnnMeta$FRiP)))
Q_FRiP <-ggplot(PreAnnMeta, aes(x=umap1, y=umap2,
                          color=FRiP)) +
  geom_point(size=0.01) +
  theme_minimal()+
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  labs(title = "FRiP ratio")+sc+
  theme(legend.key.size = unit(4, "lines"),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 25),  # Adjust size for x-axis text
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        legend.title = element_text(size = 20))

pdf(paste0(OutfilePathName,"_AnnQCPlots.pdf"), width=45, height=20)
grid.arrange(PreAnn,Re1_plot, Q_Tn5, Q_doubletscore,
             All, Re2_plot, Q_rTSS,Q_FRiP,
             ncol=4, widths=c(2.2,1,1.2,1.2))
dev.off()
}
