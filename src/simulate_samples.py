"""Simulate 820 WGS samples, each genotyped at 200 known variant loci
with a hidden true genotype (the gold-standard call) and a noisy
observed read-pileup signal -- the ground truth the variant-calling
stage must recover from imperfect sequencing data."""
from __future__ import annotations

import numpy as np
import pandas as pd

N_SAMPLES = 820
N_LOCI = 200
SEED = 863


def simulate(n_samples: int = N_SAMPLES, n_loci: int = N_LOCI, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    true_genotypes = rng.choice([0, 1, 2], (n_samples, n_loci), p=[0.6, 0.3, 0.1])  # 0=hom-ref,1=het,2=hom-alt

    # simulated read-pileup alt-allele fraction per locus: true genotype
    # determines the expected fraction, sequencing depth/noise determines
    # how clean the observed signal is.
    expected_af = true_genotypes / 2.0
    depth = rng.integers(15, 60, (n_samples, n_loci))
    observed_alt_reads = rng.binomial(depth, np.clip(expected_af, 0.02, 0.98))
    observed_af = observed_alt_reads / depth

    records = []
    for s in range(n_samples):
        for l in range(n_loci):
            records.append({"sample_id": f"SAMPLE-{s:04d}", "locus_id": f"LOCUS-{l:03d}",
                             "depth": int(depth[s, l]), "observed_af": round(float(observed_af[s, l]), 4),
                             "true_genotype": int(true_genotypes[s, l])})
    return pd.DataFrame(records)


if __name__ == "__main__":
    df = simulate()
    print(df.shape, df["sample_id"].nunique())
