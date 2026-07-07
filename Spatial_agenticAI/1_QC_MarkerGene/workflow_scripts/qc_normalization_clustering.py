import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import json
import os

parser = argparse.ArgumentParser()
parser.add_argument('--output_name', required=True)
parser.add_argument('--input_path', required=True)
parser.add_argument('--output_path', required=True)
parser.add_argument('--min_counts', type=int, default=100)
parser.add_argument('--min_genes', type=int, default=50)
parser.add_argument('--mt_pct', type=float, default=20)
parser.add_argument('--pt_pct', type=float, default=20)
parser.add_argument('--n_top_genes', type=int, default=2000)
parser.add_argument(
    '--suggest_thresholds',
    action='store_true',
    help=(
        "Dump QC distribution stats (percentiles, retained-cell counts for "
        "candidate cutoffs) as JSON and exit without filtering/clustering. "
        "Meant to be read by a human or agent to pick --min_counts/--min_genes/"
        "--mt_pct/--pt_pct before the real run."
    ),
)

args = parser.parse_args()
adata = sc.read(f"{args.input_path}/adata.h5ad")

os.chdir(args.output_path)
print("Current working directory:", os.getcwd())

if args.suggest_thresholds:
    def describe(series):
        return {
            "min": float(series.min()),
            "p1": float(series.quantile(0.01)),
            "p5": float(series.quantile(0.05)),
            "p25": float(series.quantile(0.25)),
            "median": float(series.median()),
            "p75": float(series.quantile(0.75)),
            "p95": float(series.quantile(0.95)),
            "p99": float(series.quantile(0.99)),
            "max": float(series.max()),
            "mean": float(series.mean()),
        }

    candidate_min_counts = [50, 100, 150, 200, 300, 500]
    candidate_min_genes = [25, 50, 75, 100, 150]

    stats = {
        "n_obs_raw": int(adata.n_obs),
        "n_vars_raw": int(adata.n_vars),
        "total_counts": describe(adata.obs["total_counts"]),
        "n_genes_by_counts": describe(adata.obs["n_genes_by_counts"]),
        "total_counts_Mt": describe(adata.obs["total_counts_Mt"]),
        "total_counts_PT": describe(adata.obs["total_counts_PT"]),
        "cells_retained_by_min_counts": {
            str(c): int((adata.obs["total_counts"] >= c).sum())
            for c in candidate_min_counts
        },
        "cells_retained_by_min_genes": {
            str(g): int((adata.obs["n_genes_by_counts"] >= g).sum())
            for g in candidate_min_genes
        },
    }
    stats_file = os.path.join(args.output_path, f"{args.output_name}_qc_stats.json")
    with open(stats_file, "w") as f:
        json.dump(stats, f, indent=2)
    print(json.dumps(stats, indent=2))
    raise SystemExit(0)

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


sc.pp.filter_cells(adata, min_counts=args.min_counts)
sc.pp.filter_cells(adata, min_genes=args.min_genes)
adata = adata[adata.obs["total_counts_Mt"] < args.mt_pct].copy()
adata = adata[adata.obs["total_counts_PT"] < args.pt_pct].copy()

print(f"#cells after MT filter: {adata.n_obs}")

plt.savefig(f"{args.output_name}_QC_Histogram.pdf")

## normalization
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=args.n_top_genes)

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
