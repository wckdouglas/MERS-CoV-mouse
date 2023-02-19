
IMAGE=$HOME/singularity_images/diffexpr.sif
WORKDIR=$(dirname $(pwd))
singularity exec -B $WORKDIR $IMAGE python de.py
