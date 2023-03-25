
process map_to_db {
    conda params.CONDA_ENV
    cpus params.N_CPU
    
    input:
        tuple val(id), path(r1), path(r2), val(label), val(bt2_index)

    publishDir "${params.PROJECT_PATH}/result/${id}/${label}", pattern: "*"

    output:
        tuple val(id), path(r1), path(r2), val(label), path("mapped_id.txt"), path("*.bam")

    script:
        """
        bowtie2 \
            -p ${task.cpus} -x ${bt2_index} \
            --local --mm \
            --met-file bowtie2.metrics.txt \
            --no-mixed  --no-discordant \
            -1 ${r1} -2 ${r2} \
        | samtools view -b@ ${task.cpus} \
        | filter-clipped --in-bam - --unalign --right-side 0.1 --left-side 0.1 --both-end 0.1 \
        > aligned.bam  \
        && samtools view -F4 aligned.bam | cut -f1 |sort | uniq > mapped_id.txt 
        """
}

process filter_aligned {
    // this is a tool from https://github.com/wckdouglas/fq-filter-reads
    conda params.CONDA_ENV
    cpus params.N_CPU
    
    input:
        tuple val(id), path(r1), path(r2), val(label), path(id_list)

    publishDir "${params.PROJECT_PATH}/result/${id}/${label}", pattern: "*.fq.gz"

    output:
        tuple val(id), path("r1.${label}.fq.gz"), path("r2.${label}.fq.gz")

    script:
        """
        fq-filter-reads --in-id-list ${id_list} --in-fastq ${r1} --inverse | gzip > r1.${label}.fq.gz && \
        fq-filter-reads --in-id-list ${id_list} --in-fastq ${r2} --inverse | gzip > r2.${label}.fq.gz 
        """
}
