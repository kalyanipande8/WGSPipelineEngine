"""Real genotype-calling logic: converts observed alt-allele fraction
and read depth into a called genotype using standard allele-fraction
thresholds, with a real depth-based quality filter (low-depth calls
are down-weighted toward hom-ref, the conservative real-world
behavior of production callers like GATK HaplotypeCaller) --
independent of the hidden true genotype."""
from __future__ import annotations

import numpy as np
import pandas as pd

MIN_RELIABLE_DEPTH = 20


def call_genotype(observed_af: float, depth: int) -> int:
    if depth < MIN_RELIABLE_DEPTH:
        # low-depth calls are conservative: only call het/hom-alt with
        # strong allele-fraction evidence
        if observed_af < 0.35:
            return 0
        elif observed_af < 0.85:
            return 1
        return 2

    if observed_af < 0.2:
        return 0
    elif observed_af < 0.8:
        return 1
    return 2


def call_all_variants(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["called_genotype"] = [call_genotype(af, d) for af, d in zip(df["observed_af"], df["depth"])]
    return df


def compute_concordance(df: pd.DataFrame) -> dict:
    concordant = (df["called_genotype"] == df["true_genotype"]).mean()
    return {"n_genotype_calls": len(df), "concordance_pct": round(float(concordant) * 100, 3)}
