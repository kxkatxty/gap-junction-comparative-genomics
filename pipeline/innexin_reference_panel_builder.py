#!/usr/bin/env python3
"""Build metadata and SynVoy commands for the 6 innexin reference species panel."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from pipeline.common import (
    METADATA_DIR,
    REFERENCES_DIR,
    read_csv_rows,
    write_csv,
)

DEFAULT_SPECIES_CONFIG = METADATA_DIR / "species_config.csv"
DEFAULT_OUT = METADATA_DIR / "innexin_reference_panel.csv"
SYNVOY_CMDS = METADATA_DIR / "innexin_reference_synvoy_commands.sh"

INSECT_SPECIES = (
    "Drosophila melanogaster",
    "Anopheles gambiae",
    "Aedes aegypti",
    "Schistocerca americana",
)
NEMATODE_SPECIES = (
    "Caenorhabditis elegans",
    "Caenorhabditis briggsae",
)

# Representative ortholog queries per clade (one synteny run each)
INSECT_SYNVOY_GENES = (
    {"gene": "Inx2", "accession": "Q9V427", "home": "Drosophila melanogaster", "outdir": "refpanel_dmel_Inx2"},
    {"gene": "shakB", "accession": "P33085", "home": "Drosophila melanogaster", "outdir": "refpanel_dmel_shakB"},
)
NEMATODE_SYNVOY_GENES = (
    {"gene": "inx-2", "accession": "Q9U3K5", "home": "Caenorhabditis elegans", "outdir": "refpanel_cele_inx-2"},
)


def load_innexin_species(config_path: Path) -> list[str]:
    return sorted(
        {
            row["organism"].strip()
            for row in read_csv_rows(config_path)
            if row.get("family") == "innexin" and row.get("organism", "").strip()
        }
    )


def count_reference_proteins(organism: str) -> tuple[int, str]:
    slug = organism.replace(" ", "_")
    ref_dir = REFERENCES_DIR / "innexins" / slug
    if not ref_dir.exists():
        return 0, ""
    fastas = sorted(ref_dir.glob("*.fasta"))
    genes = [f.stem.split("__")[0] for f in fastas]
    return len(fastas), ";".join(genes[:15]) + ("..." if len(genes) > 15 else "")


def panel_rows(species: list[str]) -> list[dict[str, str]]:
    quality = {row.get("organism", ""): row for row in read_csv_rows(METADATA_DIR / "gap_junction_quality_table.csv")}
    annotation = {
        (row.get("family", ""), row.get("organism", "")): row
        for row in read_csv_rows(METADATA_DIR / "family_annotation_status_connexin_vertebrate.csv")
    }
    for path in sorted(METADATA_DIR.glob("family_annotation_status*.csv")):
        for row in read_csv_rows(path):
            key = (row.get("family", ""), row.get("organism", ""))
            if key[0] == "innexin" and key not in annotation:
                annotation[key] = row

    rows: list[dict[str, str]] = []
    for organism in species:
        n_prot, gene_list = count_reference_proteins(organism)
        qual = quality.get(organism, {})
        ann = annotation.get(("innexin", organism), {})
        if organism in INSECT_SPECIES:
            clade = "insect"
        elif organism in NEMATODE_SPECIES:
            clade = "nematode"
        else:
            clade = "other"
        rows.append(
            {
                "organism": organism,
                "clade": clade,
                "reference_protein_count": str(n_prot),
                "reference_genes": gene_list,
                "annotation_category": ann.get("final_category", ""),
                "assembly_level": qual.get("assembly_level", ""),
                "annotation_available": qual.get("annotation_available", ""),
            }
        )
    return rows


def target_species(home: str, panel: list[str]) -> str:
    return ",".join(s for s in panel if s != home)


def build_synvoy_script(panel: list[str]) -> str:
    insects = [s for s in panel if s in INSECT_SPECIES]
    nematodes = [s for s in panel if s in NEMATODE_SPECIES]
    lines = [
        "#!/usr/bin/env bash",
        "# SynVoy runs: innexin reference panel only (species_config.csv).",
        "# Insect and nematode clades are run separately — they do not share microsynteny.",
        "set -euo pipefail",
        'cd "$(dirname "$0")/../../SynVoy"',
        'source "${HOME}/miniconda3/etc/profile.d/conda.sh" && conda activate synvoy_env',
        "",
    ]
    # Schistocerca excluded: large genome causes MMseqs OOM on typical laptops.
    insect_targets = "Aedes aegypti,Anopheles gambiae"
    insect_max = 2
    lines.extend(
        [
            "SYNVOY_ARGS=(--mode easy --auto_params false --multi_profile false",
            "  --n_flanking_genes 5 --adaptive_max_regions 3",
            "  --enable_smith_waterman false --exon_level_search false",
            "  --mmseqs_split_memory_limit 1G",
            "  -resume)",
            "",
            f'INSECT_TARGETS="{insect_targets}"',
            f"INSECT_MAX_GENOMES={insect_max}",
            "",
        ]
    )

    for job in INSECT_SYNVOY_GENES:
        lines.extend(
            [
                f"echo '=== insect panel: {job['gene']} ({job['accession']}) ==='",
                "nextflow run main.nf -profile standard \\",
                f"  --query_id {job['accession']} \\",
                '  --target_species "${INSECT_TARGETS}" \\',
                "  --max_genomes ${INSECT_MAX_GENOMES} \\",
                f"  --outdir results/innexin_synvoy/{job['outdir']} \\",
                '  "${SYNVOY_ARGS[@]}"',
                "",
            ]
        )

    nematode_targets = target_species("Caenorhabditis elegans", nematodes)
    for job in NEMATODE_SYNVOY_GENES:
        lines.extend(
            [
                f"echo '=== nematode panel: {job['gene']} ({job['accession']}) ==='",
                "nextflow run main.nf -profile standard \\",
                f"  --query_id {job['accession']} \\",
                f'  --target_species "{nematode_targets}" \\',
                f"  --max_genomes {max(len(nematodes) - 1, 1)} \\",
                f"  --outdir results/innexin_synvoy/{job['outdir']} \\",
                '  "${SYNVOY_ARGS[@]}"',
                "",
            ]
        )

    lines.append("echo 'Done. Plots in SynVoy/results/innexin_synvoy/refpanel_*/'")
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build innexin reference panel comparison metadata.")
    parser.add_argument("--species-config", type=Path, default=DEFAULT_SPECIES_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--write-synvoy-script", type=Path, default=SYNVOY_CMDS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    species = load_innexin_species(args.species_config)
    if len(species) < 2:
        print("Need at least 2 innexin species in species_config.csv", file=sys.stderr)
        return 1

    rows = panel_rows(species)
    write_csv(args.output, rows)

    script = build_synvoy_script(species)
    args.write_synvoy_script.parent.mkdir(parents=True, exist_ok=True)
    args.write_synvoy_script.write_text(script, encoding="utf-8")
    args.write_synvoy_script.chmod(0o755)

    print(f"Wrote {len(rows)} species -> {args.output}")
    print(f"Wrote SynVoy runner -> {args.write_synvoy_script}")
    print("\nPanel clades:")
    for clade in ("insect", "nematode"):
        members = [r["organism"] for r in rows if r["clade"] == clade]
        print(f"  {clade}: {', '.join(members)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
