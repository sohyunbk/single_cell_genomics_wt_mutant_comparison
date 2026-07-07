import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--output_name', required=True)
parser.add_argument('--input_path', required=True)
parser.add_argument('--output_path', required=True)

args = parser.parse_args()
print("Current working directory:", os.getcwd())
# Create the output directory if it does not exist
adata = sc.read(f"{args.input_path}/adata.h5ad")

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)

plt.savefig(args.output_name+"_QC_Histogram.pdf") ## Save Figure

print(f'Before filtering:\n cell - {adata.n_obs}; gene - {adata.n_vars}')       # check how many genes X cells


sc.pp.filter_cells(adata, min_counts=100)
sc.pp.filter_cells(adata, min_genes=50)
adata = adata[adata.obs["total_counts_Mt"] < 20].copy()
adata = adata[adata.obs["total_counts_PT"] < 20].copy()

print(f"#cells after MT filter: {adata.n_obs}")

plt.savefig(f"{args.output_name}_QC_Histogram.pdf")

## normalization
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Manifold embedding and clustering based on transcriptional similarity
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)

# Save UMAP plot
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save="_" + args.output_name)


plt.rcParams["figure.figsize"] = (8, 8)
spatial_coords = adata.obsm['spatial'].astype(float)
adata.obsm['spatial'] = spatial_coords
sc.pl.spatial(adata, img_key="hires", color=["clusters","total_counts", "n_genes_by_counts"], wspace=0.4, save="_"+args.output_name)
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[700, 1000, 0, 600],
    alpha=0.5,
    size=1.3,
    save="Magnify_"+args.output_name
)


adata.write(args.output_path+"/adata_processed.h5ad")
