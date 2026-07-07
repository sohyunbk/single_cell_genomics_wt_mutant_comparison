import scanpy as sc
import squidpy as sq
import pandas as pd
import os
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse
import os
import math
import matplotlib as mpl
import scipy.sparse

#### Designed for Pannel D - 2 WT and 2 bif3 together
#### Cluster 12 - pink WT CZ / Cluster 5 - brown Bif3 CZ
Adata_raw = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/A619_D/adata_raw.h5ad")
Adata = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/A619_D/adata_processed.h5ad")
Adata
Meta = Adata.obs

print(Meta.head())
print(Adata.obs['total_counts'])
WT_CZ = Meta[Meta['clusters'] == '12']
print(WT_CZ)
Bif3_CZ = Meta[Meta['clusters'] == '5' ]
print(Bif3_CZ)
WT_CZ_Barcodeslist = WT_CZ.index.tolist()
Bif3_CZ_Barcodeslist = Bif3_CZ.index.tolist()

cell_gene_count_matrix = Adata_raw.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_raw = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata.obs.index, columns=Adata.var.index)
print(cell_gene_count_raw)

cell_gene_count_matrix = Adata.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_normalized = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata.obs.index, columns=Adata.var.index)
print(cell_gene_count_normalized)

#### Get the CZ cells from the table ###

barcode_set = set(Bif3_CZ_Barcodeslist)
Bif3CZ_Count = cell_gene_count_normalized.loc[cell_gene_count_normalized.index.isin(barcode_set)]

barcode_set = set(WT_CZ_Barcodeslist)
WTCZ_Count = cell_gene_count_normalized.loc[cell_gene_count_normalized.index.isin(barcode_set)]

#### Load the target genes ###########
#Markergenes = pd.read_csv('/scratch/sb14489/3.scATAC/0.Data/MarkerGene/SpatialMarkerFinal.txt', sep='\t')
#filtered_genes = Markergenes[Markergenes['name'].str.startswith('arf')]
ARF_list = [
    "arftf4", "arftf30", "arftf18", "arftf3", "arftf20",  # Activator
    "arftf25", "arftf10", "arftf36",  # Repressor
    "arftf23", "arftf26"
]
ARF_geneID = Markergenes[Markergenes['name'].isin(ARF_list)]['geneID']
#len(gene_ids)
#Zm00001eb001720 = knox1
AllPlotgenes= ARF_geneID.tolist()
AllPlotgenes.append("Zm00001eb001720")

AllPlotgenes_symbol = ARF_list
AllPlotgenes_symbol.append("knox1")

Bif3CZ_Count['Zm00001eb433460']
WTCZ_Count['Zm00001eb433460']

Bif3CZ_Count['Zm00001eb066640']
WTCZ_Count['Zm00001eb066640']

Bif3CZ_Count['Zm00001eb001720']
WTCZ_Count['Zm00001eb001720']
