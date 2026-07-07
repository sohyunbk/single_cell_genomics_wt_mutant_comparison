// QC/normalization thresholds - override via -c config or --min_counts etc on the CLI
params.min_counts  = params.min_counts  ?: 100
params.min_genes   = params.min_genes   ?: 50
params.mt_pct      = params.mt_pct      ?: 20
params.pt_pct      = params.pt_pct      ?: 20
params.n_top_genes = params.n_top_genes ?: 2000

process process_read_data {
    output:
    stdout

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import squidpy as sq
    import pandas as pd
    import os

    # Create the output directory if it does not exist
    os.makedirs("$params.output_path", exist_ok=True)

    # Read the data
    adata = sq.read.visium("$params.input_path")

    pd.set_option('display.max_columns', None)

    # Make variable names unique and calculate QC metrics
    adata.var_names_make_unique()
    adata.var.set_index('gene_ids', inplace=True)
    adata.var.index.name = 'gene_id'

    adata.var['Mt'] = adata.var_names.str.startswith('chrMt')
    adata.var['PT'] = adata.var_names.str.startswith('chrPt')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['Mt', 'PT'], inplace=True)

    # Write the output file to the specified output path
    output_file = os.path.join("$params.output_path", "adata.h5ad")
    adata.write(output_file)
    """
}


process process_qc_preprocessing {
    input:
    stdin

    output:
    stdout

    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os
    os.chdir("$params.output_path")
    print("Current working directory:", os.getcwd())
    # Create the output directory if it does not exist
    adata = sc.read("$params.output_path"+"/adata.h5ad")
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1], )
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2]) ## okay line
    sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3], )
    plt.savefig("$params.output_name"+"_QC_Histogram.pdf") ## Save Figure

    #print("Before filtering: Cell - "+ str(cell - adata.n_obs) +"gene - "+str(adata.n_vars))       # check how many genes X cells
    sc.pp.filter_cells(adata, min_counts=$params.min_counts)
    sc.pp.filter_cells(adata, min_genes=$params.min_genes)

    adata = adata[adata.obs["total_counts_Mt"] < $params.mt_pct].copy()
    adata = adata[adata.obs["total_counts_PT"] < $params.pt_pct].copy()
    #print(f"#cells after MT filter: {adata.n_obs}")
    plt.savefig("$params.output_name"+"_QC_Histogram.pdf")
    ## normalization
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=$params.n_top_genes)
    # Manifold embedding and clustering based on transcriptional similarity
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)
    # Save UMAP plot
    plt.rcParams["figure.figsize"] = (4, 4)
    sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save="_" + "$params.output_name")
    plt.rcParams["figure.figsize"] = (8, 8)
    spatial_coords = adata.obsm['spatial'].astype(float)
    adata.obsm['spatial'] = spatial_coords
    sc.pl.spatial(adata, img_key="hires", color=["clusters","total_counts", "n_genes_by_counts"], wspace=0.4, save="_"+"$params.output_name")
    sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[700, 1000, 0, 600],
    alpha=0.5,
    size=1.3,
    save="Magnify_"+"$params.output_name")
    adata.write("$params.output_path"+"/adata_processed.h5ad")
    """
}

process markergene {
    input:
    stdin

    output:
    stdout

    script:
    """
    python3 "${params.ScriptDir}marker_gene_testing.py" --output_name "$params.output_name" --input_path "$params.output_path" --markergenelist "$params.MarkerGene"
    """
    }


workflow {
    process_read_data | process_qc_preprocessing | markergene | view
    }
