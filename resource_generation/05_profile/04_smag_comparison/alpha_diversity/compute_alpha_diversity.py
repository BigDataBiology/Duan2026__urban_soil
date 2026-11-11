import sys
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from os import makedirs

import polars as pl
import polars.selectors as cs

# Run from this directory. The large OTU tables are not tracked by git; they live
# in the repository-level intermediate_data/ tree (mirroring resource_generation/).
OTU_DIR = "../../../../intermediate_data/05_profile/04_smag_comparison/alpha_diversity/OTU"


hmms_and_names = pl.read_csv("hmms_and_names", separator="\t")
hmms_and_names = hmms_and_names.with_columns(pl.col("type").alias("domain"))

# Note that one file is comma-separated and the other is tab-separated

otu_sm = f"{OTU_DIR}/PRJNA983538.csv"
otu_sm = pl.read_csv(otu_sm, separator=",")
otu_sm = otu_sm.join(hmms_and_names["name", "domain"], left_on="gene", right_on="name")

otu_sh = f"{OTU_DIR}/shanghai_soil_otu.tsv"
otu_sh = pl.read_csv(otu_sh, separator="\t")
otu_sh = otu_sh.join(hmms_and_names["name", "domain"], left_on="gene", right_on="name")

otu_dog = f"{OTU_DIR}/Dog_otu.tsv"
otu_dog = pl.read_csv(otu_dog, separator="\t")
otu_dog = otu_dog.join(
    hmms_and_names["name", "domain"], left_on="gene", right_on="name"
)

# Use only bacterial and both domains for alpha diversity computation
KEEP_DOMAINS = ["both", "bacteria"]
otu_sm = otu_sm.filter(pl.col("domain").is_in(KEEP_DOMAINS))
otu_sh = otu_sh.filter(pl.col("domain").is_in(KEEP_DOMAINS))
otu_dog = otu_dog.filter(pl.col("domain").is_in(KEEP_DOMAINS))

AGGREGATIONS = ["mean", "tmean10", "median"]
METRICS = ["Shannon", "Chao1", "Simpson"]
GROUP_KEYS = ["gene", "sample"]


def compute_gene_sample_metrics(otu):
    metrics = otu.group_by(GROUP_KEYS).agg(
        pl.col("num_hits").entropy().alias("Shannon"),
        (1 - ((pl.col("num_hits") / pl.col("num_hits").sum()) ** 2).sum()).alias(
            "Simpson"
        ),
        pl.len().alias("S_obs"),
        pl.col("num_hits").eq(1).sum().alias("F1"),
        pl.col("num_hits").eq(2).sum().alias("F2"),
    )
    return metrics.with_columns(
        pl.when(pl.col("F2") > 0)
        .then(pl.col("S_obs") + (pl.col("F1") ** 2) / (2 * pl.col("F2")))
        .otherwise(pl.col("S_obs") + (pl.col("F1") * (pl.col("F1") - 1)) / 2)
        .alias("Chao1")
    ).select(GROUP_KEYS + METRICS)


def aggregate_metric(metric_values, metric, aggregation):
    values = pl.col(metric)

    if aggregation == "mean":
        agg_f = values.mean()
    elif aggregation == "tmean10":
        agg_f = values.clip(values.quantile(0.1), values.quantile(0.9)).mean()
    elif aggregation == "median":
        agg_f = values.median()
    else:
        raise ValueError(f"Unknown aggregation method: {aggregation}")
    return metric_values.group_by("sample").agg(agg_f.alias(metric))


makedirs("figures", exist_ok=True)
makedirs("results", exist_ok=True)
sh_metrics = compute_gene_sample_metrics(otu_sh)
smag_metrics = compute_gene_sample_metrics(otu_sm)
dog_metrics = compute_gene_sample_metrics(otu_dog)

# Group label -> per gene x sample metrics (order used for plots)
GROUPS = {
    "Urban soil": sh_metrics,
    "Smag": smag_metrics,
    "Dog": dog_metrics,
}
GROUP_ORDER = list(GROUPS)
COMPARISONS = [("Urban soil", "Smag"), ("Smag", "Dog"), ("Urban soil", "Dog")]

per_sample_means = {}
for metric in METRICS:
    fig, axes = plt.subplots(1, len(AGGREGATIONS), figsize=(10, 4), sharey=True)
    data = []
    for ax, agg in zip(axes, AGGREGATIONS):
        ax.clear()
        print(f"{metric} aggregation method: {agg}")
        aggregated = {
            group: aggregate_metric(metrics, metric, agg)
            for group, metrics in GROUPS.items()
        }
        for g1, g2 in COMPARISONS:
            stat, p_value = stats.mannwhitneyu(
                aggregated[g1][metric], aggregated[g2][metric], alternative="two-sided"
            )
            print(
                f"{g1} vs {g2}: Mann-Whitney U test: statistic={stat}, p-value={p_value}"
            )

        cur = pl.concat(
            [
                values.with_columns(pl.lit(group).alias("Group"))
                for group, values in aggregated.items()
            ]
        ).with_columns(pl.lit(agg).alias("Aggregation"))
        sns.boxplot(
            x="Group",
            y=metric,
            ax=ax,
            data=cur,
            fill=False,
            color="black",
            fliersize=0,
            order=GROUP_ORDER,
        )
        sns.stripplot(
            x="Group",
            y=metric,
            ax=ax,
            data=cur,
            hue="Group",
            size=2,
            jitter=True,
            order=GROUP_ORDER,
            alpha=0.7,
        )
        if ax.legend_ is not None:
            ax.legend_.remove()

        ax.set_title(f"{agg}")
        data.append(cur)
    sns.despine(fig, trim=True)
    fig.tight_layout()
    fig.savefig(f"figures/{metric.lower()}.svg")

    data = pl.concat(data)
    data = data.pivot(index="sample", on="Aggregation", values=metric).sort("sample")
    data.write_csv(f"results/{metric.lower()}.tsv", separator="\t", float_precision=8)
    per_sample_means[metric] = data.select("sample", pl.col("mean").alias(metric))

    print(
        f"\n# Spearman correlation between different aggregation methods for {metric}\n"
    )
    print("Agg1\tAgg2\tSpearman correlation")
    for col in AGGREGATIONS:
        for col2 in AGGREGATIONS:
            if col > col2:
                print(
                    f"{col}\t{col2}\t{stats.spearmanr(data[col], data[col2]).statistic:.4f}"
                )


# Combined table (per-sample mean across genes) used by analysis/figure1/figure1.ipynb
alpha_diversity = per_sample_means["Chao1"]
for metric in ["Simpson", "Shannon"]:
    alpha_diversity = alpha_diversity.join(per_sample_means[metric], on="sample")
alpha_diversity.sort("sample").write_csv(
    "alpha_diversity.tsv", separator="\t", float_precision=8
)


print("# Spot check computation")

SPOT_SAMPLE = "CPSNJ01_350"
spot_metrics = []
sel_1sample = otu_sh.filter(pl.col("sample") == SPOT_SAMPLE)
for g in sel_1sample["gene"].unique():
    sel = sel_1sample.filter(pl.col("gene") == g)
    counts = sel["num_hits"]
    n_singletons = sum(v == 1 for v in counts)
    n_doubletons = sum(v == 2 for v in counts)
    s_obs = len(counts)
    if n_doubletons > 0:
        chao1 = s_obs + n_singletons**2 / (2 * n_doubletons)
    else:
        chao1 = s_obs + n_singletons * (n_singletons - 1) / 2
    total = sum(counts)
    simpson = 1 - sum((v / total) ** 2 for v in counts)
    spot_metrics.append((g, stats.entropy(counts), chao1, simpson))
spot_metrics = pl.DataFrame(
    spot_metrics, schema=["gene", "Shannon", "Chao1", "Simpson"], orient="row"
)
orig_metrics = sh_metrics.filter(pl.col("sample") == SPOT_SAMPLE).sort("gene")
recomp_metrics = spot_metrics.sort("gene")

print("Difference between original and recomputed per-gene values:")
diff = orig_metrics.select(cs.float()) - recomp_metrics.select(cs.float())
print(diff)
print("Max absolute difference:")
print(diff.select(pl.col("*").abs().max()))
