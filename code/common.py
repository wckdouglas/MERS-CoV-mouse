import csv
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Tuple, Optional, Dict

import pandas as pd

DATA_PATH = Path(__file__).resolve().parents[1] / "data"
RESULT_PATH = Path(__file__).resolve().parents[1] / "result"
REFERENCE_PATH = Path(__file__).resolve().parents[1] / "ref"
FIGURE_PATH = Path(__file__).resolve().parents[1] / "figures"

ISG_GENES = ["Isg15", "Isg20", "lfit1", "lfit2", "lfit3", "Bst2", "Oas3", "Oasl1", "Mx1", "Rsad2"]


@dataclass(frozen=True)
class Sample:
    sample_id: str
    infection: Optional[str]
    antibody: Optional[str]
    time_point: str


@dataclass
class Comparison:
    filter: str
    label: str
    comparison: Tuple[Tuple[str, str, str], str]
    de_label: str


COMPARISONS = {
    "comparison_1": Comparison(
        filter="Antibody == 'α-CCR2' | Antibody == 'control'",
        label="α-CCR vs control antibody",
        comparison=(("Antibody", "α-CCR2", "control"), "5d"),
        de_label="Antibody: α-CCR2 vs control at 5d",
    ),
    "comparison_2": Comparison(
        filter="(Infection == 'MERS-CoV wt' | Infection == 'MERS-CoV mutant') & time_point == '2d'",
        label="MERS-Cov mutant vs WT at 2d",
        comparison=(("Infection", "MERS-CoV mutant", "MERS-CoV wt"), "2d"),
        de_label="Infection: MERS-CoV mutant vs MERS-CoV wt at 2d",
    ),
    "comparison_3": Comparison(
        filter="(Infection == 'MERS-CoV wt' | Infection == 'MERS-CoV mutant') & time_point == '5d'",
        label="MERS-Cov mutant vs WT at 5d",
        comparison=(("Infection", "MERS-CoV mutant", "MERS-CoV wt"), "5d"),
        de_label="Infection: MERS-CoV mutant vs MERS-CoV wt at 5d",
    ),
}


@cache
def get_tx2gene():
    return pd.read_csv(REFERENCE_PATH / "tx2gene.tsv", sep="\t", names=["tid", "gid"])


def get_sample_table():
    return (
        pd.read_csv(DATA_PATH / "2023-01-27_metadata.csv")
        .query('species=="Mouse"')
        .query('species=="Mouse"')
        .assign(id=lambda d: d["sample_id"])
        .set_index("id")
    )


@cache
def get_expression_data():
    return (
        pd.read_csv(RESULT_PATH / "norm_count.csv")
        .rename(columns={"variable": "sample_id", "index": "gene_id"})
        .merge(get_sample_table())
        .merge(get_gene_table())
    )


@cache
def get_gene_table():
    return pd.read_csv(DATA_PATH / "gene_table.csv").rename(
        columns={
            "Gene type": "gene_type",
            "Gene stable ID": "gene_id",
            "Gene name": "gene_name",
        }
    )


def genic_express(gene_list):
    gene_df = get_expression_data()
    assert set(gene_list).issubset(set(gene_df["gene_name"]))
    return gene_df.pipe(lambda d: d[d["gene_name"].isin(gene_list)])


def comparison_genic_expression(gene_list, comparison):
    if comparison not in COMPARISONS.keys():
        raise ValueError(f"comparison must be one of {COMPARISONS.keys()}")
    gene_data = genic_express(gene_list)
    return gene_data.query(COMPARISONS[comparison].filter)


def read_meta() -> Dict[str, Sample]:
    sample_data = {}
    with open(DATA_PATH  /  "2023-01-27_metadata.csv") as metadata:
        for row in csv.DictReader(metadata):
            sample_data[row['sample_id']] = Sample(
                sample_id=row['sample_id'],
                infection=row['Infection'],
                antibody=row['Antibody'],
                time_point=row['time_point'],
            )
    return sample_data