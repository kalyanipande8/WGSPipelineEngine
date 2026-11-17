import json
import os


def main():
    from src.simulate_samples import simulate
    from src.variant_caller import call_all_variants, compute_concordance
    from src.pipeline_stages import benchmark_runtime

    df = simulate()
    called = call_all_variants(df)
    concordance = compute_concordance(called)
    runtime = benchmark_runtime(df)

    summary = {"n_samples": df["sample_id"].nunique(), "n_loci": df["locus_id"].nunique(),
               "concordance": concordance, "runtime_optimization": runtime}

    os.makedirs("results", exist_ok=True)
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
