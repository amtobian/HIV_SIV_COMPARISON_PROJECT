import argparse
import csv
import random
from collections import Counter
from itertools import combinations
from math import ceil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REGIONS = ("env", "pol", "gag")
BASES = frozenset("ACGT")
OUTPUT_FIELDS = (
    "comparison_type",
    "gene_a",
    "group_a",
    "gene_b",
    "group_b",
    "pi_per_site_a",
    "pi_per_site_b",
    "difference_a_minus_b",
    "ci95_lower",
    "ci95_upper",
    "p_value",
    "holm_adjusted_p_value",
    "n_sequences_a",
    "n_sequences_b",
    "callable_columns_a",
    "callable_columns_b",
    "alignment_blocks_a",
    "alignment_blocks_b",
    "valid_bootstrap_replicates",
    "bootstrap_replicates",
    "block_size_columns",
    "random_seed",
)


def read_alignment(path):
    sequences_by_group = {}
    seen_identifiers = set()
    identifier = None
    sequence_parts = []

    def save_record():
        if identifier is None:
            return
        sequence = "".join(sequence_parts).upper().replace("U", "T")
        if not sequence:
            raise ValueError(f"Empty sequence for {identifier!r} in {path}")
        group, separator, _ = identifier.partition("__")
        if not separator or not group:
            raise ValueError(
                f"Expected a group-prefixed FASTA identifier in {path}: "
                f"{identifier!r}"
            )
        if group == "ALL":
            raise ValueError("'ALL' is reserved for the pooled profile")
        sequences_by_group.setdefault(group, []).append(sequence)

    with path.open(encoding="utf-8") as alignment_file:
        for line_number, line in enumerate(alignment_file, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_record()
                identifier = line[1:].split(None, 1)[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                if identifier in seen_identifiers:
                    raise ValueError(
                        f"Duplicate FASTA identifier in {path}: {identifier}"
                    )
                seen_identifiers.add(identifier)
                sequence_parts = []
            elif identifier is None:
                raise ValueError(f"Sequence before FASTA header at {path}:{line_number}")
            else:
                sequence_parts.append(line)
    save_record()

    if not sequences_by_group:
        raise ValueError(f"No FASTA records found in {path}")

    lengths = {
        len(sequence)
        for sequences in sequences_by_group.values()
        for sequence in sequences
    }
    if len(lengths) != 1:
        raise ValueError(f"Aligned sequences do not all have the same length in {path}")
    alignment_length = lengths.pop()
    if alignment_length == 0:
        raise ValueError(f"Alignment contains zero-length sequences in {path}")

    return sequences_by_group, alignment_length


def diversity_by_column(sequences):
    profile = []
    for column in zip(*sequences):
        counts = Counter(base for base in column if base in BASES)
        callable_count = sum(counts.values())
        pair_count = callable_count * (callable_count - 1) // 2
        if pair_count == 0:
            profile.append(None)
            continue
        identical_pairs = sum(
            count * (count - 1) // 2 for count in counts.values()
        )
        profile.append((pair_count - identical_pairs) / pair_count)
    return profile


def summarize_blocks(profile, block_size):
    block_count = ceil(len(profile) / block_size)
    if block_count < 2:
        raise ValueError(
            f"Alignment has {len(profile)} columns, too few for a block bootstrap "
            f"with block size {block_size}"
        )

    blocks = []
    for index in range(block_count):
        start = index * len(profile) // block_count
        end = (index + 1) * len(profile) // block_count
        values = [value for value in profile[start:end] if value is not None]
        blocks.append((sum(values), len(values)))
    callable_columns = sum(count for _, count in blocks)
    if callable_columns == 0:
        return None

    callable_blocks = sum(count > 0 for _, count in blocks)
    if callable_blocks < 2:
        return None
    return {
        "block_sums": [total for total, _ in blocks],
        "block_counts": [count for _, count in blocks],
        "mean_pi": sum(total for total, _ in blocks) / callable_columns,
        "callable_columns": callable_columns,
        "alignment_blocks": block_count,
    }


def bootstrap_profiles(group_summaries, replicates, rng):
    group_names = sorted(group_summaries)
    block_count = group_summaries[group_names[0]]["alignment_blocks"]
    profiles = {group: [] for group in group_names}

    for _ in range(replicates):
        sampled_blocks = rng.choices(range(block_count), k=block_count)
        for group in group_names:
            summary = group_summaries[group]
            numerator = sum(summary["block_sums"][i] for i in sampled_blocks)
            denominator = sum(summary["block_counts"][i] for i in sampled_blocks)
            profiles[group].append(
                numerator / denominator if denominator else None
            )
    return profiles


def percentile(values, probability):
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def compare_profiles(observed_a, observed_b, bootstrap_a, bootstrap_b, replicates):
    difference = observed_a - observed_b
    draws = [
        value_a - value_b
        for value_a, value_b in zip(bootstrap_a, bootstrap_b)
        if value_a is not None and value_b is not None
    ]
    minimum_valid = ceil(replicates * 0.95)
    if len(draws) < minimum_valid:
        raise ValueError(
            f"Only {len(draws)} of {replicates} bootstrap replicates had callable "
            f"sites for both profiles; at least {minimum_valid} are required"
        )

    centered_null = [draw - difference for draw in draws]
    extreme = sum(abs(value) >= abs(difference) for value in centered_null)
    p_value = (extreme + 1) / (len(centered_null) + 1)
    return {
        "difference": difference,
        "ci95_lower": percentile(draws, 0.025),
        "ci95_upper": percentile(draws, 0.975),
        "p_value": p_value,
        "valid_replicates": len(draws),
    }


def holm_adjust(p_values):
    adjusted = [0.0] * len(p_values)
    ordered = sorted(enumerate(p_values), key=lambda item: item[1])
    running_max = 0.0
    count = len(p_values)
    for rank, (index, p_value) in enumerate(ordered):
        running_max = max(running_max, (count - rank) * p_value)
        adjusted[index] = min(1.0, running_max)
    return adjusted


def analyze(results_dir, block_size, replicates, seed):
    summaries_by_gene = {}
    bootstrap_by_gene = {}
    skipped = []

    for region_index, region in enumerate(REGIONS):
        alignment_path = results_dir / f"{region}.aln.fasta"
        if not alignment_path.is_file():
            raise FileNotFoundError(f"Required alignment not found: {alignment_path}")
        sequences_by_group, _ = read_alignment(alignment_path)
        group_summaries = {}
        for group, sequences in sequences_by_group.items():
            if len(sequences) < 2:
                skipped.append(
                    f"{region}/{group}: fewer than two sequences; diversity is undefined"
                )
                continue
            summary = summarize_blocks(diversity_by_column(sequences), block_size)
            if summary is None:
                skipped.append(
                    f"{region}/{group}: fewer than two callable alignment blocks"
                )
                continue
            summary["n_sequences"] = len(sequences)
            group_summaries[group] = summary

        if group_summaries:
            summaries_by_gene[region] = group_summaries
            bootstrap_by_gene[region] = bootstrap_profiles(
                group_summaries,
                replicates,
                random.Random(seed + region_index),
            )

    comparisons = []

    def add_comparison(kind, gene_a, group_a, gene_b, group_b):
        summary_a = summaries_by_gene[gene_a][group_a]
        summary_b = summaries_by_gene[gene_b][group_b]
        result = compare_profiles(
            summary_a["mean_pi"],
            summary_b["mean_pi"],
            bootstrap_by_gene[gene_a][group_a],
            bootstrap_by_gene[gene_b][group_b],
            replicates,
        )
        comparisons.append(
            {
                "comparison_type": kind,
                "gene_a": gene_a,
                "group_a": group_a,
                "gene_b": gene_b,
                "group_b": group_b,
                "pi_per_site_a": summary_a["mean_pi"],
                "pi_per_site_b": summary_b["mean_pi"],
                "difference_a_minus_b": result["difference"],
                "ci95_lower": result["ci95_lower"],
                "ci95_upper": result["ci95_upper"],
                "p_value": result["p_value"],
                "n_sequences_a": summary_a["n_sequences"],
                "n_sequences_b": summary_b["n_sequences"],
                "callable_columns_a": summary_a["callable_columns"],
                "callable_columns_b": summary_b["callable_columns"],
                "alignment_blocks_a": summary_a["alignment_blocks"],
                "alignment_blocks_b": summary_b["alignment_blocks"],
                "valid_bootstrap_replicates": result["valid_replicates"],
                "bootstrap_replicates": replicates,
                "block_size_columns": block_size,
                "random_seed": seed,
            }
        )

    for region, group_summaries in summaries_by_gene.items():
        for group_a, group_b in combinations(sorted(group_summaries), 2):
            add_comparison(
                "viral_groups_within_gene",
                region,
                group_a,
                region,
                group_b,
            )

    for gene_a, gene_b in combinations(REGIONS, 2):
        if gene_a not in summaries_by_gene or gene_b not in summaries_by_gene:
            continue
        shared_groups = (
            set(summaries_by_gene[gene_a]) & set(summaries_by_gene[gene_b])
        )
        for group in sorted(shared_groups):
            add_comparison(
                "genes_within_viral_group",
                gene_a,
                group,
                gene_b,
                group,
            )

    if not comparisons:
        raise ValueError("No valid diversity comparisons could be calculated")

    adjusted_p_values = holm_adjust(
        [comparison["p_value"] for comparison in comparisons]
    )
    for comparison, adjusted_p in zip(comparisons, adjusted_p_values):
        comparison["holm_adjusted_p_value"] = adjusted_p

    return comparisons, skipped


def write_results(path, comparisons):
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        for comparison in comparisons:
            writer.writerow(
                {
                    field: (
                        f"{comparison[field]:.8f}"
                        if field in {
                            "pi_per_site_a",
                            "pi_per_site_b",
                            "difference_a_minus_b",
                            "ci95_lower",
                            "ci95_upper",
                            "p_value",
                            "holm_adjusted_p_value",
                        }
                        else comparison[field]
                    )
                    for field in OUTPUT_FIELDS
                }
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare mean per-site nucleotide diversity among genes and viral groups "
            "using an alignment-block bootstrap."
        )
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=PROJECT_ROOT / "results",
        help="Directory containing env, pol, and gag aligned FASTA files.",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=100,
        help="Approximate bootstrap block width in alignment columns (default: 100).",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=10000,
        help="Number of bootstrap replicates (default: 10000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=314159,
        help="Random seed for reproducibility (default: 314159).",
    )
    args = parser.parse_args()
    if args.block_size < 1 or args.bootstrap_replicates < 1:
        parser.error("--block-size and --bootstrap-replicates must be positive")

    comparisons, skipped = analyze(
        args.results_dir,
        args.block_size,
        args.bootstrap_replicates,
        args.seed,
    )
    output_path = args.results_dir / "statistical_diversity_comparisons.csv"
    write_results(output_path, comparisons)
    for message in skipped:
        print(f"Skipped {message}")
    print(f"Wrote {len(comparisons)} comparisons to {output_path}")


if __name__ == "__main__":
    main()
