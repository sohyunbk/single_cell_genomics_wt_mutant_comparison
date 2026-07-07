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
    python3 "${params.ScriptDir}read_data.py" --input_path "$params.input_path" --output_path "$params.output_path"
    """
}


process process_qc_preprocessing {
    input:
    stdin

    output:
    stdout

    script:
    """
    python3 "${params.ScriptDir}qc_normalization_clustering.py" \
        --output_name "$params.output_name" \
        --input_path "$params.output_path" \
        --output_path "$params.output_path" \
        --min_counts $params.min_counts \
        --min_genes $params.min_genes \
        --mt_pct $params.mt_pct \
        --pt_pct $params.pt_pct \
        --n_top_genes $params.n_top_genes
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
