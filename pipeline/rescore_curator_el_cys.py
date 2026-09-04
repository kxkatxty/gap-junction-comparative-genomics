#!/usr/bin/env python3
"""Re-score existing curator hits with EL Cys 2+2 motif (no remapping).

For each species with hits.csv + miniprot.gff (or locate region GFFs), re-translate
overlapping models, apply the updated classify_product (EL motif aware), rewrite
hits.csv, then rebuild innexin_search_summary.
"""

from __future__ import annotations

import csv
from pathlib import Path

from pipeline.batch_innexin_curator_probe import (
    GENOME_ROOT,
    OUT_ROOT,
    Hit,
    assess_species,
    classify_product,
    find_genome_fasta,
    parse_miniprot_models,
    prior_discovery,
    summarize,
    translate_model,
    write_species_hits,
)
from pipeline.common import PROJECT_ROOT
from pipeline.innexin_topology import evaluate_el_cys_motif


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def _rescore_from_gff(genome: Path, gff: Path) -> list[Hit]:
    return assess_species(genome, gff)


def _annotate_existing_hits(genome: Path, gff: Path, old_hits: list[dict]) -> list[Hit]:
    """Keep hit coordinates; re-translate best overlapping GFF model for motif + verdict."""
    models = parse_miniprot_models(gff)
    out: list[Hit] = []
    for row in old_hits:
        seqid = row["seqid"]
        start, end = int(row["start"]), int(row["end"])
        identity = float(row.get("identity") or 0)
        best_m = None
        best_ov = 0
        for m in models:
            if m["seqid"] != seqid:
                continue
            ov = _overlap(start, end, int(m["start"]), int(m["end"]))
            if ov > best_ov:
                best_ov = ov
                best_m = m
        if best_m is None or best_ov < 50:
            # keep row, mark motif unresolved
            motif = evaluate_el_cys_motif("")  # unresolved empty
            out.append(
                Hit(
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=row.get("strand") or "+",
                    identity=identity,
                    target=row.get("target") or "",
                    aa_length=int(float(row.get("aa_length") or 0)),
                    cys=int(float(row.get("cys") or 0)),
                    stops=int(float(row.get("stops") or 0)),
                    verdict=row.get("verdict") or "weak",
                    el_motif="unresolved",
                    tm_pred="0",
                )
            )
            continue
        try:
            prot = translate_model(genome, best_m)
        except Exception:
            out.append(
                Hit(
                    seqid=seqid,
                    start=start,
                    end=end,
                    strand=row.get("strand") or "+",
                    identity=identity,
                    target=row.get("target") or best_m.get("target") or "",
                    aa_length=int(float(row.get("aa_length") or 0)),
                    cys=int(float(row.get("cys") or 0)),
                    stops=int(float(row.get("stops") or 0)),
                    verdict=row.get("verdict") or "weak",
                    el_motif="unresolved",
                    tm_pred="0",
                )
            )
            continue
        clean = prot.split("*")[0]
        verdict, motif = classify_product(float(best_m["identity"]), prot)
        if verdict == "noise":
            continue
        out.append(
            Hit(
                seqid=seqid,
                start=start,
                end=end,
                strand=str(best_m["strand"]),
                identity=float(best_m["identity"]),
                target=str(best_m.get("target") or ""),
                aa_length=len(clean),
                cys=clean.count("C"),
                stops=prot.count("*"),
                verdict=verdict,
                el1_cys="" if motif.el1_cys is None else str(motif.el1_cys),
                el2_cys="" if motif.el2_cys is None else str(motif.el2_cys),
                el_motif=motif.status,
                tm_pred=str(motif.tm_count),
            )
        )
    return out


def _pick_gff(sp_dir: Path) -> Path | None:
    gff = sp_dir / "miniprot.gff"
    if gff.exists() and gff.stat().st_size > 0:
        return gff
    # Drosophila / locate-only: concatenate region GFFs if present
    regions = sp_dir / "locate_tblastn_regions" / "regions"
    if regions.is_dir():
        parts = sorted(regions.glob("*/miniprot.gff"))
        if parts:
            merged = sp_dir / "locate_tblastn_regions" / "merged_miniprot.gff"
            with merged.open("w", encoding="utf-8") as out:
                for p in parts:
                    out.write(p.read_text(encoding="utf-8"))
            return merged
    return None


def main() -> int:
    summary_rows: list[dict[str, str]] = []
    existing: dict[str, dict[str, str]] = {}
    csv_path = OUT_ROOT / "innexin_search_summary.csv"
    if csv_path.exists():
        with csv_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                existing[row["species"]] = row

    present_species = [
        r["species"]
        for r in existing.values()
        if int(r.get("present_loci") or 0) > 0 or r.get("status") == "searched"
    ]
    # Prefer all searched with hits
    targets = sorted(
        {
            p.name
            for p in OUT_ROOT.iterdir()
            if p.is_dir() and (p / "hits.csv").exists()
        }
    )

    motif_counts = {"pass": 0, "fail": 0, "unresolved": 0}

    for species in targets:
        sp_dir = OUT_ROOT / species
        genome = find_genome_fasta(GENOME_ROOT / species)
        gff = _pick_gff(sp_dir)
        prior, prior_accepted = prior_discovery(species)
        old = list(csv.DictReader((sp_dir / "hits.csv").open(encoding="utf-8", newline="")))
        if genome is None or gff is None:
            # keep existing; no motif columns
            hits = [
                Hit(
                    seqid=r["seqid"],
                    start=int(r["start"]),
                    end=int(r["end"]),
                    strand=r.get("strand") or "+",
                    identity=float(r.get("identity") or 0),
                    target=r.get("target") or "",
                    aa_length=int(float(r.get("aa_length") or 0)),
                    cys=int(float(r.get("cys") or 0)),
                    stops=int(float(r.get("stops") or 0)),
                    verdict=r.get("verdict") or "weak",
                    el_motif="unresolved",
                    tm_pred="",
                )
                for r in old
                if r.get("verdict") != "noise"
            ]
        else:
            hits = _annotate_existing_hits(genome, gff, old)

        write_species_hits(sp_dir, hits)
        for h in hits:
            if h.el_motif in motif_counts:
                motif_counts[h.el_motif] += 1

        present = [h for h in hits if h.verdict == "present"]
        frag = [h for h in hits if h.verdict == "fragmentary"]
        weak = [h for h in hits if h.verdict == "weak"]
        best = max(hits, key=lambda h: (h.aa_length, h.identity), default=None)
        el_pass = sum(1 for h in present if h.el_motif == "pass")
        row = {
            "species": species,
            "status": "searched",
            "genome": str(genome.relative_to(PROJECT_ROOT)) if genome else existing.get(species, {}).get("genome", ""),
            "prior_discovery": prior,
            "prior_accepted": prior_accepted,
            "miniprot_mrna_count": existing.get(species, {}).get("miniprot_mrna_count", str(len(hits))),
            "present_loci": str(len(present)),
            "fragmentary_loci": str(len(frag)),
            "weak_loci": str(len(weak)),
            "best_identity": f"{best.identity:.4f}" if best else "",
            "best_aa": str(best.aa_length) if best else "",
            "best_locus": f"{best.seqid}:{best.start}-{best.end}" if best else "",
            "locate_method": existing.get(species, {}).get("locate_method", ""),
            "notes": (
                f"present_loci floor; el_motif_pass_among_present={el_pass}/"
                f"{len(present)}; rescored_el_cys"
            ),
        }
        summary_rows.append(row)
        print(
            f"{species}: present=≥{len(present)} frag={len(frag)} "
            f"el_pass={el_pass} motif_hits={sum(1 for h in hits if h.el_motif=='pass')}"
        )

    summarize(summary_rows)
    print("EL motif tallies across hits:", motif_counts)
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
