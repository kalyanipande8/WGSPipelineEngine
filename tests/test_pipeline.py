import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.simulate_samples import simulate
from src.variant_caller import call_all_variants, compute_concordance


def test_at_least_800_samples():
    df = simulate()
    assert df["sample_id"].nunique() >= 800


def test_concordance_is_high():
    df = simulate()
    called = call_all_variants(df)
    concordance = compute_concordance(called)
    assert concordance["concordance_pct"] >= 85.0


def test_low_depth_calls_are_conservative():
    from src.variant_caller import call_genotype
    assert call_genotype(0.5, 10) == 1  # moderate AF at low depth calls het, not hom-alt
    assert call_genotype(0.3, 10) == 0


def test_high_depth_calls_follow_standard_thresholds():
    from src.variant_caller import call_genotype
    assert call_genotype(0.05, 50) == 0
    assert call_genotype(0.5, 50) == 1
    assert call_genotype(0.95, 50) == 2


def test_genotype_calls_are_valid_values():
    df = simulate()
    called = call_all_variants(df)
    assert called["called_genotype"].isin([0, 1, 2]).all()
