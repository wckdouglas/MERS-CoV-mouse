import logging
from enum import Enum
from pathlib import Path
from typing import List, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandera as pa
import seaborn as sns
from pandera.typing import DataFrame, Series
from pydantic.dataclasses import dataclass
from ridgeplot.colors import ColorEncoder, ColorPalette

from common import COMPARISONS, comparison_genic_expression, read_meta

logging.basicConfig(level=logging.INFO)

LOG = logging.getLogger(__file__)

plt.rc("xtick", labelsize=15)
plt.rc("ytick", labelsize=15)
plt.rc("axes", labelsize=15)
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 1

DATA_PATH = Path(__file__).parents[1]


class DESeqDataFrame(pa.SchemaModel):
    gene_id: Series[str]
    baseMean: Series[float]
    label: Series[str]
    gene_name: Series[str]
    gene_type: Series[str]
    log2FoldChange: Series[float] = pa.Field(nullable=True)
    lfcSE: Series[float] = pa.Field(nullable=True)
    stat: Series[float] = pa.Field(nullable=True)
    pvalue: Series[float] = pa.Field(nullable=True)
    padj: Series[float] = pa.Field(nullable=True)


@pa.check_types
def get_significant_protein_gene(
    df: DataFrame[DESeqDataFrame],
    lfc_threshold: Union[int, float] = LFC_THRESHOLD,
    padj_threshold: float = PADJ_THRESHOLD,
) -> DataFrame[DESeqDataFrame]:
    """
    Picking out significant protein genes

    :param df: DESeq2 result dataframe
    :param lfc_threshold: threshold to cutoff absolute log fold change values
    :param padj_threshold: threshold to cutoff adjusted p value
    """
    data = (
        df.query(f"padj < {padj_threshold}")
        .query("gene_type=='protein_coding'")
        .query(f"log2FoldChange > {lfc_threshold} | log2FoldChange < {-lfc_threshold}")
    )
    return data


def plot_volcano(ax, protein_genes_df: DataFrame[DESeqDataFrame]):
    for sig, sig_df in protein_genes_df.assign(
        Significant=lambda d: list(map(label_significant, d.itertuples()))
    ).groupby("Significant"):
        color, alpha = Significant.__members__[sig].value
        ax.scatter(
            sig_df["log2FoldChange"],
            -np.log2(sig_df["padj"]),
            color=color,
            alpha=alpha,
            label=sig,
        )
    ax.set_title(size=15, label=sig_df["label"].tolist()[0])
    ax.set_ylabel("-log2 adjusted p-value ")
    ax.set_xlabel("Fold change (log2)")
    sns.despine()


def plot_ma(ax, protein_genes_df: DataFrame[DESeqDataFrame]):
    for sig, sig_df in protein_genes_df.assign(
        Significant=lambda d: list(map(label_significant, d.itertuples()))
    ).groupby("Significant"):
        color, alpha = Significant.__members__[sig].value
        ax.scatter(
            np.log2(sig_df["baseMean"] + 1),
            sig_df["log2FoldChange"],
            color=color,
            alpha=alpha,
            label=sig,
        )
    ax.set_title(size=15, label=sig_df["label"].tolist()[0])
    ax.hlines(y=0, xmin=-1, xmax=22, color="red")
    ax.set_xlabel("Expression level (log2 base mean)")
    ax.set_ylabel("Fold change (log2)")
    sns.despine()


class Significant(Enum):
    Significant = ("red", 0.3)
    No_change = ("lightgray", 0.3)


def label_significant(row) -> str:
    if abs(row.log2FoldChange) > LFC_THRESHOLD and row.padj < PADJ_THRESHOLD:
        return Significant.Significant.name
    else:
        return Significant.No_change.name


@dataclass
class Pathway:
    name: str
    description_url: str
    genes: list[str]


def index_gmt_pathways(gmt_file):
    pathway_index: Dict[str, Pathway] = dict()
    with Path(gmt_file).open("r") as f:
        for pathway_count, line in enumerate(f):
            fields = line.strip().split("\t")
            pathway = Pathway(
                name=fields[0], description_url=fields[1], genes=fields[2:]
            )
            pathway_index[pathway.name] = pathway
    LOG.info("Indexed %i pathways", pathway_count + 1)
    return pathway_index


def plot_heatmap(
    de_table,
    comparison,
    gene_list=None,
    figsize=(10, 20),
    sample_coloring="infection",
    diff_expr_gene_only=False,
    zscore_normalize=False,
):
    sample_data = read_meta()
    color_palette = ColorPalette["okabeito"]
    if diff_expr_gene_only:
        gene_list_diff = (
            de_table.pipe(lambda d: d[d["label"] == COMPARISONS[comparison].de_label])
            .pipe(
                lambda d: get_significant_protein_gene(
                    d, lfc_threshold=1, padj_threshold=0.05
                )
            )["gene_name"]
            .tolist()
        )

        if gene_list is not None:
            gene_list = set(gene_list_diff).intersection(gene_list)
        else:
            gene_list = gene_list_diff

    norm_df = comparison_genic_expression(gene_list, comparison, allow_missing=True)
    plot_mat = norm_df.pipe(
        pd.pivot, index="gene_name", columns="sample_id", values="value"
    )

    labels = [sample_data[col].__dict__[sample_coloring] for col in plot_mat.columns]
    color_encoder = ColorEncoder()
    color_encoder.fit(labels, color_palette)

    p = sns.clustermap(
        np.log10(plot_mat + 1),
        col_colors=color_encoder.transform(labels),
        figsize=figsize,
        cbar_kws={
            "label": "$log_{10}$ normalized expression",
        },
        cmap="viridis",
        col_cluster=False,
        colors_ratio=0.02,
        cbar_pos=(1, 0.2, 0.1, 0.6),
        z_score=0 if zscore_normalize else None,
    )
    color_encoder.show_legend(
        ax=p.ax_heatmap, bbox_to_anchor=(1.05, 1.1), frameon=False, fontsize=15
    )
    yt = p.ax_heatmap.set_yticks(np.arange(plot_mat.shape[0]) + 0.5)
    row_order = p.dendrogram_row.reordered_ind
    yt = p.ax_heatmap.set_yticklabels(plot_mat.index.values[row_order])
    # p.ax_heatmap.yaxis.set_visible(False)
    p.ax_heatmap.set_ylabel("")
    p.ax_heatmap.set_xlabel("")

    if zscore_normalize:
        cbar_label = "Z-score, $log_{10}$ normalized expression"
    else:
        cbar_label = "$log_{10}$ normalized expression"
    p.ax_cbar.set_ylabel(cbar_label, va="bottom", rotation=270, fontsize=15)
    return p
