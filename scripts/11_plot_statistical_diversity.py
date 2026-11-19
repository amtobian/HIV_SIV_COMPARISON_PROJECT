import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "results" / "statistical_diversity_comparisons.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "statistical_diversity_forest_plot.svg"
REQUIRED_COLUMNS = {
    "comparison_type",
    "gene_a",
    "group_a",
    "gene_b",
    "group_b",
    "difference_a_minus_b",
    "ci95_lower",
    "ci95_upper",
    "holm_adjusted_p_value",
}
COMPARISON_LABELS = {
    "viral_groups_within_gene": "Viral groups within each gene",
    "genes_within_viral_group": "Genes within each viral group",
}
SIGNIFICANT_COLOR = "#b91c1c"
NOT_SIGNIFICANT_COLOR = "#2563eb"


def read_comparisons(path):
    comparisons = []
    with path.open(newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must contain columns: "
                + ", ".join(sorted(REQUIRED_COLUMNS))
            )

        for line_number, row in enumerate(reader, start=2):
            comparison_type = row["comparison_type"]
            if comparison_type not in COMPARISON_LABELS:
                raise ValueError(
                    f"Unexpected comparison type in {path}:{line_number}: "
                    f"{comparison_type!r}"
                )
            if not all(row[field].strip() for field in (
                "gene_a", "group_a", "gene_b", "group_b"
            )):
                raise ValueError(f"Missing comparison label in {path}:{line_number}")
            try:
                difference = float(row["difference_a_minus_b"])
                lower = float(row["ci95_lower"])
                upper = float(row["ci95_upper"])
                adjusted_p = float(row["holm_adjusted_p_value"])
            except ValueError as error:
                raise ValueError(
                    f"Invalid statistical value in {path}:{line_number}"
                ) from error
            if lower > upper or not 0 <= adjusted_p <= 1:
                raise ValueError(
                    f"Invalid confidence interval or adjusted p-value "
                    f"in {path}:{line_number}"
                )
            comparisons.append(
                {
                    "comparison_type": comparison_type,
                    "gene_a": row["gene_a"],
                    "group_a": row["group_a"],
                    "gene_b": row["gene_b"],
                    "group_b": row["group_b"],
                    "difference": difference,
                    "lower": lower,
                    "upper": upper,
                    "adjusted_p": adjusted_p,
                }
            )

    if not comparisons:
        raise ValueError(f"No statistical comparisons found in {path}")
    return comparisons


def format_p_value(value):
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def render_plot(comparisons):
    by_type = defaultdict(list)
    for comparison in comparisons:
        by_type[comparison["comparison_type"]].append(comparison)
    ordered_types = (
        "viral_groups_within_gene",
        "genes_within_viral_group",
    )
    max_abs_value = max(
        abs(value)
        for comparison in comparisons
        for value in (comparison["lower"], comparison["upper"])
    )
    x_limit = max_abs_value * 1.12 if max_abs_value else 0.01

    width = 1120
    left = 300
    right = 920
    top = 72
    row_height = 25
    section_height = 25
    subsection_height = 22
    rows = []
    for comparison_type in ordered_types:
        items = by_type.get(comparison_type, [])
        if not items:
            continue
        rows.append(("section", COMPARISON_LABELS[comparison_type], None))
        subgroup_field = (
            "gene_a"
            if comparison_type == "viral_groups_within_gene"
            else "group_a"
        )
        subgroups = []
        for item in items:
            subgroup = item[subgroup_field]
            if subgroup not in subgroups:
                subgroups.append(subgroup)
        for subgroup in subgroups:
            rows.append(("subsection", subgroup, None))
            rows.extend(
                ("comparison", "", item)
                for item in items
                if item[subgroup_field] == subgroup
            )

    height = (
        top
        + sum(
            section_height if kind == "section"
            else subsection_height if kind == "subsection"
            else row_height
            for kind, _, _ in rows
        )
        + 95
    )

    def x_position(value):
        return left + (value + x_limit) / (2 * x_limit) * (right - left)

    parts = [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}" viewBox="0 0 {width} {height}" role="img" '
            'aria-labelledby="chart-title chart-description">'
        ),
        "<title id=\"chart-title\">Statistical diversity comparisons</title>",
        (
            "<desc id=\"chart-description\">Forest plot of estimated differences "
            "in mean nucleotide diversity with 95% block-bootstrap confidence "
            "intervals. Red points have Holm-adjusted p-values below 0.05.</desc>"
        ),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        (
            f'<text x="{width / 2}" y="32" text-anchor="middle" '
            'font-family="sans-serif" font-size="21" font-weight="bold">'
            "Statistical diversity comparisons</text>"
        ),
        (
            f'<text x="{(left + right) / 2}" y="55" text-anchor="middle" '
            'font-family="sans-serif" font-size="12">'
            "Difference in mean nucleotide diversity (A − B; π/site)</text>"
        ),
        (
            f'<text x="{right + 30}" y="55" text-anchor="middle" '
            'font-family="sans-serif" font-size="12">Holm-adjusted p</text>'
        ),
    ]

    y = top
    plot_bottom = height - 79
    tick_values = [
        -x_limit,
        -x_limit / 2,
        0,
        x_limit / 2,
        x_limit,
    ]
    for value in tick_values:
        x = x_position(value)
        grid_color = "#9ca3af" if value == 0 else "#e5e7eb"
        grid_width = ' stroke-width="1.5"' if value == 0 else ""
        parts.extend(
            [
                (
                    f'<line x1="{x:.2f}" y1="{top - 4}" x2="{x:.2f}" '
                    f'y2="{plot_bottom}" stroke="{grid_color}"{grid_width}/>'
                ),
                (
                    f'<text x="{x:.2f}" y="{plot_bottom + 20}" '
                    'text-anchor="middle" font-family="sans-serif" font-size="11">'
                    f"{value:.3f}</text>"
                ),
            ]
        )

    for kind, label, comparison in rows:
        if kind == "section":
            y += section_height
            parts.append(
                f'<text x="20" y="{y - 6}" font-family="sans-serif" '
                f'font-size="15" font-weight="bold">{html.escape(label)}</text>'
            )
            continue
        if kind == "subsection":
            y += subsection_height
            parts.append(
                f'<text x="38" y="{y - 5}" font-family="sans-serif" '
                f'font-size="12" font-weight="bold">{html.escape(label.upper())}</text>'
            )
            continue

        y += row_height
        if comparison["comparison_type"] == "viral_groups_within_gene":
            row_label = f'{comparison["group_a"]} vs {comparison["group_b"]}'
        else:
            row_label = f'{comparison["gene_a"]} vs {comparison["gene_b"]}'
        color = (
            SIGNIFICANT_COLOR
            if comparison["adjusted_p"] < 0.05
            else NOT_SIGNIFICANT_COLOR
        )
        low_x = x_position(comparison["lower"])
        high_x = x_position(comparison["upper"])
        center_x = x_position(comparison["difference"])
        parts.extend(
            [
                (
                    f'<text x="54" y="{y}" dominant-baseline="middle" '
                    'font-family="sans-serif" font-size="11">'
                    f"{html.escape(row_label)}</text>"
                ),
                (
                    f'<line x1="{low_x:.2f}" y1="{y}" x2="{high_x:.2f}" y2="{y}" '
                    f'stroke="{color}" stroke-width="2"/>'
                ),
                (
                    f'<line x1="{low_x:.2f}" y1="{y - 4}" x2="{low_x:.2f}" '
                    f'y2="{y + 4}" stroke="{color}" stroke-width="1.5"/>'
                ),
                (
                    f'<line x1="{high_x:.2f}" y1="{y - 4}" x2="{high_x:.2f}" '
                    f'y2="{y + 4}" stroke="{color}" stroke-width="1.5"/>'
                ),
                (
                    f'<circle cx="{center_x:.2f}" cy="{y}" r="4" fill="{color}">'
                    f'<title>Difference {comparison["difference"]:.5f}; 95% CI '
                    f'[{comparison["lower"]:.5f}, {comparison["upper"]:.5f}]</title>'
                    "</circle>"
                ),
                (
                    f'<text x="{right + 30}" y="{y}" text-anchor="middle" '
                    'dominant-baseline="middle" font-family="sans-serif" '
                    f'font-size="11">{format_p_value(comparison["adjusted_p"])}</text>'
                ),
            ]
        )

    legend_y = height - 42
    parts.extend(
        [
            (
                f'<circle cx="{left + 4}" cy="{legend_y}" r="4" '
                f'fill="{SIGNIFICANT_COLOR}"/>'
            ),
            (
                f'<text x="{left + 14}" y="{legend_y}" dominant-baseline="middle" '
                'font-family="sans-serif" font-size="11">'
                "Holm-adjusted p &lt; 0.05</text>"
            ),
            (
                f'<circle cx="{left + 200}" cy="{legend_y}" r="4" '
                f'fill="{NOT_SIGNIFICANT_COLOR}"/>'
            ),
            (
                f'<text x="{left + 210}" y="{legend_y}" dominant-baseline="middle" '
                'font-family="sans-serif" font-size="11">'
                "Holm-adjusted p ≥ 0.05</text>"
            ),
            (
                f'<text x="{width / 2}" y="{height - 13}" text-anchor="middle" '
                'font-family="sans-serif" font-size="10" fill="#4b5563">'
                "Positive values indicate higher diversity for A; negative values "
                "indicate higher diversity for B. Intervals reflect alignment-block "
                "bootstrap uncertainty.</text>"
            ),
            "</svg>",
        ]
    )
    return "\n".join(parts) + "\n"


def main():
    parser = argparse.ArgumentParser(
        description="Plot statistical diversity comparisons from their CSV output."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    comparisons = read_comparisons(args.input)
    args.output.write_text(render_plot(comparisons), encoding="utf-8")
    print(f"Wrote statistical diversity plot to {args.output}")


if __name__ == "__main__":
    main()
