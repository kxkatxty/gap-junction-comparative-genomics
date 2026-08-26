"""Shared resource caps for the pannexin workspace.

Keep defaults gentle so laptop rebuilds do not thrash the machine.
Heavy genome/tree jobs must be opted into explicitly.
"""

from __future__ import annotations

import os

# Hard cap: never more than 2 worker threads unless the user overrides.
MAX_THREADS = 2


def n_threads(default: int = 1) -> int:
    """Return a small thread count (1–MAX_THREADS)."""
    raw = os.environ.get("PANX_THREADS", "").strip()
    if raw.isdigit():
        return max(1, min(int(raw), MAX_THREADS))
    return max(1, min(default, MAX_THREADS))


# Opt-in switches (env or CLI). Off by default.
def run_heavy_phylogeny() -> bool:
    return os.environ.get("PANX_HEAVY_PHYLO", "").strip() in {"1", "true", "yes"}


def run_genome_curator() -> bool:
    return os.environ.get("PANX_RUN_MINIPROT", "").strip() in {"1", "true", "yes"}
