import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse
import os
import math
import matplotlib as mpl

parser = argparse.ArgumentParser()
parser.add_argument('--output_name', required=True)
parser.add_argument('--input_path', required=True)
parser.add_argument('--markergenelist', required=True)
args = parser.parse_args()


adata = sc.read(f"{args.input_path}/adata_processed.h5ad")
df = pd.read_csv(args.markergenelist, sep='\t')

os.chdir(args.input_path)

gene_list = df['geneID'].tolist()
#gene_symbols = dict(zip(df['geneID'], df['geneID']+"_"+df['name']))
gene_symbols = dict(zip(df['geneID'], df['name']))


gene_ids_in_adata = adata.var.index.values
filtered_gene_list = [gene for gene in gene_list if gene in gene_ids_in_adata]

#filtered_gene_list = filtered_gene_list[1:10]

num_genes = len(filtered_gene_list)
ncols = 10  # Increase the number of columns to fit more plots in each row
nrows = math.ceil(num_genes / ncols)

# Create the spatial plot with larger figure size
#fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, nrows * 2.5))  # Adjust figsize for larger plots
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(80, nrows * 10))  # Adjust figsize for larger plots

axes = axes.flatten()  # Flatten the axes array for easy iteration

for i, gene in enumerate(filtered_gene_list):
    sc.pl.spatial(
        adata,
        color=gene,
        ax=axes[i],  # Plot on the specific subplot axis
        show=False  # Do not show immediately, we will save it manually
    )
    # Set custom title using gene_symbols dictionary with smaller font size
    if gene in gene_symbols:
        axes[i].set_title(f'{gene}\n{gene_symbols[gene]}', fontsize=8)  # Smaller font size
    else:
        axes[i].set_title(gene, fontsize=8)  # Use gene ID if symbol not found, with smaller font size
    # Customize the scale bar
    for child in axes[i].get_children():
        if isinstance(child, mpl.collections.PatchCollection):
            for path in child.get_paths():
                if path.vertices.shape[0] == 5:  # Typical of the scale bar
                    path.vertices *= 0.1  # Scale down the size
                    break

# Hide any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Adjust layout
plt.tight_layout()

# Close all figures to avoid the RuntimeWarning
plt.savefig(f"{args.input_path}MarkerGeneAll_{args.output_name}.pdf")
plt.close('all')
