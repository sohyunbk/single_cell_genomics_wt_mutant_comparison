#!/usr/bin/env python3
"""Agent step: choose spatial-QC thresholds from the QC distribution.

Reads the JSON produced by `qc_normalization_clustering.py --suggest_thresholds`,
asks Claude to pick min_counts / min_genes / mt_pct / pt_pct / n_top_genes based
on the observed per-spot distribution, and writes two files:

  <output_path>/<output_name>_chosen_thresholds.json   # the numbers qc reads back
  <output_path>/<output_name>_agent_decision.json      # full audit trail

The agent only selects numeric QC thresholds -- it never edits pipeline code or
changes the analysis method. Requires ANTHROPIC_API_KEY in the environment
(or an `ant auth login` profile).
"""
import argparse
import datetime
import json
import os
import sys

import anthropic

SYSTEM = (
    "You are a spatial transcriptomics QC assistant for 10x Visium data. "
    "Given summary statistics of a sample's per-spot QC metrics, choose filtering "
    "thresholds for a Scanpy pipeline. You control only these numeric knobs:\n"
    "  - min_counts: minimum total UMI counts per spot (sc.pp.filter_cells)\n"
    "  - min_genes: minimum number of genes detected per spot\n"
    "  - mt_pct: drop spots whose mitochondrial percentage is >= this value\n"
    "  - pt_pct: drop spots whose plastid percentage is >= this value\n"
    "  - n_top_genes: number of highly variable genes kept for clustering\n\n"
    "Choose thresholds that remove low-quality / empty spots and likely-dying "
    "spots without discarding real tissue. Base every choice on the provided "
    "distribution (percentiles, retained-spot counts) rather than on generic "
    "defaults, and be conservative: when in doubt, keep more spots. Do not "
    "invent metrics you were not given."
)

# Structured-output schema. Note: structured outputs do not support numeric
# min/max constraints, so ranges are validated in Python after the call.
SCHEMA = {
    "type": "object",
    "properties": {
        "min_counts": {"type": "integer"},
        "min_genes": {"type": "integer"},
        "mt_pct": {"type": "number"},
        "pt_pct": {"type": "number"},
        "n_top_genes": {"type": "integer"},
        "reasoning": {
            "type": "string",
            "description": "Brief justification tied to the observed distribution.",
        },
    },
    "required": [
        "min_counts",
        "min_genes",
        "mt_pct",
        "pt_pct",
        "n_top_genes",
        "reasoning",
    ],
    "additionalProperties": False,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats_json", required=True,
                        help="Path to the *_qc_stats.json from --suggest_thresholds.")
    parser.add_argument("--output_path", required=True)
    parser.add_argument("--output_name", required=True)
    parser.add_argument("--model", default="claude-opus-4-8",
                        help="Claude model id (default: claude-opus-4-8).")
    args = parser.parse_args()

    with open(args.stats_json) as f:
        stats = json.load(f)

    client = anthropic.Anthropic()  # ANTHROPIC_API_KEY / ant profile from env

    prompt = (
        "Here are the QC distribution statistics for one Visium sample. "
        "Choose filtering thresholds based on them.\n\n"
        f"{json.dumps(stats, indent=2)}"
    )

    response = client.messages.create(
        model=args.model,
        max_tokens=4000,
        thinking={"type": "adaptive"},
        system=SYSTEM,
        messages=[{"role": "user", "content": prompt}],
        output_config={"format": {"type": "json_schema", "schema": SCHEMA}},
    )

    if response.stop_reason == "refusal":
        sys.exit("Claude declined the request; set thresholds manually and rerun "
                 "the non-agentic pipeline.")

    text = next((b.text for b in response.content if b.type == "text"), None)
    if text is None:
        sys.exit("No text block in the model response; cannot read thresholds.")

    text = text.strip()
    if text.startswith("```"):  # tolerate a fenced ```json block just in case
        text = text.strip("`")
        text = text[text.find("{"):]
    decision = json.loads(text)

    # Sanity-check the numbers before letting them drive a destructive filter.
    errs = []
    if not decision["min_counts"] >= 1:
        errs.append("min_counts must be >= 1")
    if not decision["min_genes"] >= 1:
        errs.append("min_genes must be >= 1")
    if not 0 < decision["mt_pct"] <= 100:
        errs.append("mt_pct must be in (0, 100]")
    if not 0 < decision["pt_pct"] <= 100:
        errs.append("pt_pct must be in (0, 100]")
    if not decision["n_top_genes"] >= 1:
        errs.append("n_top_genes must be >= 1")
    if errs:
        sys.exit("Agent returned invalid thresholds: " + "; ".join(errs))

    thresholds = {
        k: decision[k]
        for k in ("min_counts", "min_genes", "mt_pct", "pt_pct", "n_top_genes")
    }

    chosen_path = os.path.join(
        args.output_path, f"{args.output_name}_chosen_thresholds.json"
    )
    with open(chosen_path, "w") as f:
        json.dump(thresholds, f, indent=2)

    try:
        usage = response.usage.model_dump()
    except Exception:
        usage = None
    audit = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "model": response.model,
        "request_id": getattr(response, "_request_id", None),
        "usage": usage,
        "stats_json": os.path.abspath(args.stats_json),
        "input_stats": stats,
        "chosen_thresholds": thresholds,
        "reasoning": decision["reasoning"],
    }
    audit_path = os.path.join(
        args.output_path, f"{args.output_name}_agent_decision.json"
    )
    with open(audit_path, "w") as f:
        json.dump(audit, f, indent=2)

    print(f"Agent ({response.model}) chose thresholds: {json.dumps(thresholds)}")
    print(f"Reasoning: {decision['reasoning']}")
    print(f"Wrote {chosen_path}")
    print(f"Wrote audit trail {audit_path}")


if __name__ == "__main__":
    main()
