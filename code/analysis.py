from enum import Enum
from pathlib import Path
from pydantic.dataclasses import dataclass
from typing import Union, List
import matplotlib.pyplot as plt
import numpy as np
import pandera as pa
import seaborn as sns
from pandera.typing import DataFrame, Series
import logging

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
