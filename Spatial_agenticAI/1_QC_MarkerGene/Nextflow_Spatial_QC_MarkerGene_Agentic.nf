// Agentic spatial-QC pipeline.
//
// Flow:
//   process_read_data          read the Space Ranger outs/ -> adata.h5ad
//   process_suggest_thresholds dump the QC distribution as JSON
//   process_agent_pick_thresholds  Claude reads that JSON and chooses thresholds
//   process_qc_preprocessing   filter/normalize/cluster with the chosen thresholds
//   markergene                 marker-gene testing
//
// The agent only chooses numeric QC thresholds (min_counts / min_genes /
// mt_pct / pt_pct / n_top_genes); it never edits pipeline code or changes the
// analysis method. Every decision is written to <output_name>_agent_decision.json
// as an audit trail.
//
// Requires ANTHROPIC_API_KEY in the environment for the agent step.
// Set params.agent_model to pick the Claude model (default: claude-opus-4-8).

params.agent_model = params.agent_model ?: 'claude-opus-4-8'


process process_read_data {
    output:
    stdout

    script:
    """
    python3 "${params.ScriptDir}read_data.py" --input_path "$params.input_path" --output_path "$params.output_path"
    """
}


process process_suggest_thresholds {
    input:
    stdin

    output:
    stdout

    script:
    """
    python3 "${params.ScriptDir}qc_normalization_clustering.py" \
        --suggest_thresholds \
        --output_name "$params.output_name" \
        --input_path "$params.output_path" \
        --output_path "$params.output_path"
    """
}


process process_agent_pick_thresholds {
    input:
    stdin

    output:
    stdout

    script:
    """
    python3 "${params.ScriptDir}agent_pick_thresholds.py" \
        --stats_json "${params.output_path}/${params.output_name}_qc_stats.json" \
        --output_path "$params.output_path" \
        --output_name "$params.output_name" \
        --model "$params.agent_model"
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
        --thresholds_json "${params.output_path}/${params.output_name}_chosen_thresholds.json"
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
    process_read_data \
        | process_suggest_thresholds \
        | process_agent_pick_thresholds \
        | process_qc_preprocessing \
        | markergene \
        | view
}
