
IMAGE=$HOME/singularity_images/diffexpr.sif
WORKDIR=/raid/home/wu58/projects/rna-seq/mouse/
singularity exec -B $WORKDIR $IMAGE python de.py
