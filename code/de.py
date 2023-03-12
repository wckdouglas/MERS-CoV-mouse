import logging

import pandas as pd
from diffexpr.py_deseq import py_DESeq2

from common import COMPARISONS, RESULT_PATH, get_sample_table, get_tx2gene

logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

REMOVE_SAMPLES = ["CCR2_2"]


def run_deseq(h5_dict, subsampled_metadata, formula, contrast):
    h5_dict_subset = {k: h5_dict[k].as_posix() for k in subsampled_metadata.sample_id}
    tx2gene = get_tx2gene()
    dds = py_DESeq2(
        count_matrix=h5_dict_subset,
        design_matrix=subsampled_metadata,
        design_formula=formula,
        kallisto=True,
        tx2gene=tx2gene,
        threads=4,
    )

    dds.run_deseq()
    dds.get_deseq_result(contrast=list(contrast))
    res = dds.deseq_result
    norm_count = dds.normalized_count().reset_index().pipe(pd.melt, id_vars="index")
    print(norm_count.head())
    return res, norm_count


def make_deseq(h5_dict, sample_table, contrast, day):
    if len(contrast) != 3:
        raise ValueError(
            "contrast must be a 3 value tuple: [column name, factor 1, factor 2]"
        )

    column = contrast[0]
    factors = contrast[1:]
    subsampled_metadata = (
        sample_table.pipe(
            lambda d: d[(d[column].isin(factors)) & (d["time_point"] == day)]
        )
        .pipe(lambda d: d[~d["sample_id"].isin(REMOVE_SAMPLES)])
        .filter(["sample_id", column])
    )
    return run_deseq(h5_dict, subsampled_metadata, f"~ {column}", contrast)


if __name__ == "__main__":
    LOG.info("Getting data from %s", RESULT_PATH)

    # sample metadata table
    sample_table = get_sample_table()

    # collect all kallisto results
    h5_dict = {}
    for h5 in RESULT_PATH.glob("*/kallisto/abundance.h5"):
        h5_dict[h5.parents[1].name] = h5
    LOG.info("Collected %i h5 files from kallisto", len(h5_dict))

    # run deseq2 for each comparison and concat
    de_df_list, norm_count_df_list = [], []
    for comparison_item in COMPARISONS.values():
        comparison = list(comparison_item.comparison)
        LOG.info(f"Running deseq2 for comparison {comparison}")
        label = "{}: {} vs {} at {}".format(*comparison[0], comparison[1])
        de_df, norm_count_df = make_deseq(
            h5_dict, sample_table, comparison[0], comparison[1]
        )
        de_df_list.append(de_df.assign(label=label))
        norm_count_df_list.append(norm_count_df)

    de_df = pd.concat(de_df_list).reset_index()
    de_table_name = RESULT_PATH / "de_result.csv"
    de_df.to_csv(de_table_name, index=False)
    LOG.info("Written %s" % de_table_name)

    norm_count_df = pd.concat(norm_count_df_list).reset_index()
    norm_count_table_name = RESULT_PATH / "norm_count.csv"
    norm_count_df.to_csv(norm_count_table_name, index=False)
    LOG.info("Written %s" % norm_count_table_name)
