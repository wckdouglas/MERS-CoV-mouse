# RNA-seq analysis for MERS-CoV infected mouse cells #

## Make reference ##

The `makefile` contains all procedures to create the mm10 reference.

```
cd code
make -j 4
```


## Calculate gene expression with kallisto ## 


The gene expression calculation is done with the nextflow pipeline `gene_expr_workflow.nf`:

```
cd code
bash run_pipeline.sh
```

##  Differential expression with DESeq2 ##

Differential expression is calculated by DESeq2 (via [diffexpr](https://github.com/wckdouglas/diffexpr))


```
singularity pull -F $HOME/singularity_image/diffexpr.sif docker://ghcr.io/wckdouglas/diffexpr/diffexpr-dev:60638aee31131566dd18b0601ce0c3433eeb5e47
cd code
bash run_de.sh
```
