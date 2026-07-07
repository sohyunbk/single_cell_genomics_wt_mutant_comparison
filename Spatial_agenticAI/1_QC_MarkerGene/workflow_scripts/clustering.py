import scanpy as sc
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', required=True)
args = parser.parse_args()

adata = sc.read(args.input_file)

plt.rcParams["figure.figsize"] = (8, 8)
spatial_coords = adata.obsm['spatial'].astype(float)
adata.obsm['spatial'] = spatial_coords
sc.pl.spatial(adata, img_key="hires", color=["clusters","total_counts", "n_genes_by_counts"], wspace=0.4, save="_"+output_name)
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[700, 1000, 0, 600],
    alpha=0.5,
    size=1.3,
    save="Magnify_"+output_name
)

adata.write("adata_clustered.h5ad")
