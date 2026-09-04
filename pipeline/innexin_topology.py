"""Innexin topology helpers: TM spans + extracellular-loop cysteine motif.

Innexin / pannexin monomers have four transmembrane helices with N-cytoplasmic
topology. The conserved cysteine ladder is **2 Cys in EL1** (between TM1–TM2)
and **2 Cys in EL2** (between TM3–TM4), not merely Cys≥4 anywhere in the chain.

TM spans are predicted with the same Kyte–Doolittle sliding-window heuristic used
in ``tools/discover_innexins.py`` (no Phobius/TMHMM required). Predictions on
truncated / non-Met products can be wrong — callers should treat
``el_motif_status == "unresolved"`` as "could not test", not "failed".
"""

from __future__ import annotations

from dataclasses import dataclass

KYTE_DOOLITTLE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}


@dataclass(frozen=True)
class TmSpan:
    start: int  # inclusive, 0-based
    end: int  # exclusive


@dataclass(frozen=True)
class ElCysMotif:
    tm_count: int
    el1_cys: int | None
    el2_cys: int | None
    total_cys: int
    status: str  # "pass" | "fail" | "unresolved"
    detail: str

    @property
    def ok(self) -> bool:
        return self.status == "pass"


def predict_tm_spans(
    protein: str,
    *,
    window: int = 19,
    threshold: float = 1.15,
    min_gap: int = 5,
) -> list[TmSpan]:
    """Return non-overlapping hydrophobic spans (approximate TM helices)."""
    clean = protein.replace("*", "")
    if len(clean) < window:
        return []
    scores = [
        sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in clean[i : i + window]) / window
        for i in range(len(clean) - window + 1)
    ]
    raw: list[TmSpan] = []
    i = 0
    while i < len(scores):
        if scores[i] >= threshold:
            j = i
            while j < len(scores) and scores[j] >= threshold:
                j += 1
            # helix covers residues i .. j+window-2 inclusive → end exclusive j+window-1
            raw.append(TmSpan(start=i, end=min(len(clean), j + window - 1)))
            i = j
        else:
            i += 1
    # merge near-adjacent peaks
    if not raw:
        return []
    merged: list[TmSpan] = [raw[0]]
    for span in raw[1:]:
        prev = merged[-1]
        if span.start <= prev.end + min_gap:
            merged[-1] = TmSpan(start=prev.start, end=max(prev.end, span.end))
        else:
            merged.append(span)
    return merged


def evaluate_el_cys_motif(
    protein: str,
    *,
    min_per_loop: int = 2,
    loop_pad: int = 6,
) -> ElCysMotif:
    """Test 2+2 extracellular Cys using predicted TM1–TM4 topology.

    When Kyte–Doolittle over-segments (>4 peaks), try every consecutive
    4-helix window. Loop windows are padded by ``loop_pad`` aa into adjacent
    TM edges so boundary cysteines are not swallowed by hydrophobic spans.
    """
    clean = (protein[:-1] if protein.endswith("*") else protein).split("*")[0]
    total = clean.count("C")
    spans = predict_tm_spans(clean)
    n_tm = len(spans)
    if n_tm < 4:
        return ElCysMotif(
            tm_count=n_tm,
            el1_cys=None,
            el2_cys=None,
            total_cys=total,
            status="unresolved",
            detail=f"predicted_tm={n_tm}<4; cannot place EL1/EL2",
        )

    n = len(clean)
    best: tuple[int, int, str] | None = None
    for i in range(n_tm - 3):
        tm = spans[i : i + 4]
        el1_a = max(0, tm[0].end - loop_pad)
        el1_b = min(n, tm[1].start + loop_pad)
        el2_a = max(0, tm[2].end - loop_pad)
        el2_b = min(n, tm[3].start + loop_pad)
        if el1_b <= el1_a or el2_b <= el2_a:
            continue
        el1_c = clean[el1_a:el1_b].count("C")
        el2_c = clean[el2_a:el2_b].count("C")
        detail = (
            f"EL1={el1_c}cys EL2={el2_c}cys window_tm={i}:{i+3} "
            f"pad={loop_pad} (need ≥{min_per_loop}+≥{min_per_loop})"
        )
        if el1_c >= min_per_loop and el2_c >= min_per_loop:
            return ElCysMotif(
                tm_count=n_tm,
                el1_cys=el1_c,
                el2_cys=el2_c,
                total_cys=total,
                status="pass",
                detail=detail,
            )
        if best is None or (el1_c + el2_c) > (best[0] + best[1]):
            best = (el1_c, el2_c, detail)

    el1_c, el2_c, detail = best if best else (0, 0, "no_valid_tm_window")
    return ElCysMotif(
        tm_count=n_tm,
        el1_cys=el1_c,
        el2_cys=el2_c,
        total_cys=total,
        status="fail",
        detail=detail,
    )
