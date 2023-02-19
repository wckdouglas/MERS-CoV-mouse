nextflow.enable.dsl=2

params.PROJECT_PATH="/home/wu58/projects/rna-seq/mouse"
params.N_CPU = 2
params.CONDA_ENV = "${params.PROJECT_PATH}/code/conda.yml"

DATA_PATH="${params.PROJECT_PATH}/data"
REFERENCE_PATH="${params.PROJECT_PATH}/ref"
INDEX_PATH="$REFERENCE_PATH/mm10.transcripts.kallisto.idx"
GENE_GTF="$REFERENCE_PATH/genes.gtf"
GENOME_INDEX_PATH="$REFERENCE_PATH/mm10.masked.fa"
VIRAL_GENOME="$REFERENCE_PATH/viral.fa"
RRNA="$REFERENCE_PATH/rRNA.fa"
METADATA_TABLE = "$DATA_PATH/2023-01-27_metadata.csv"
FASTQ_PATH = "/data/NGS/2023_01/mRNA_seq_Qiwen_Roy"


include {map_to_db as map_rRNA; map_to_db as map_viral} from './modules/filter_reads'
include {filter_aligned as filter_rRNA; filter_aligned as filter_viral} from './modules/filter_reads'

process cutadapt {
    conda params.CONDA_ENV
    cpus params.N_CPU
    
    input: 
        tuple val(id), path(r1), path(r2)

    publishDir "${params.PROJECT_PATH}/result/${id}/trimmed", pattern: "*"

    output:
        tuple val(id), path("r1.trimmed.fq.gz"), path("r2.trimmed.fq.gz")

    script:
        """
        cutadapt \
        --cores ${task.cpus} \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-a "A{100}" \
	-A "A{100}" \
	-a "T{100}" \
	-A "T{100}" \
        -o r1.trimmed.fq.gz \
        -p r2.trimmed.fq.gz \
        --quality-cutoff 30 \
        --minimum-length 30 \
        --report full \
        $r1 $r2 &> trim.log
        """
}


process run_kallisto {
    conda params.CONDA_ENV
    cpus params.N_CPU

    input:
        tuple val(id), path(r1), path(r2)

    publishDir "${params.PROJECT_PATH}/result/${id}", pattern: "kallisto", mode: "copy"

    output:
        tuple val(id), path("kallisto/pseudoalignments.bam")

    script:
        """
        kallisto quant \
            --genomebam \
            --gtf ${REFERENCE_PATH}/genes.gtf \
            -i $INDEX_PATH \
            -t ${task.cpus} \
            -o kallisto ${r1} ${r2}
        """
}

process index_kallisto {
    conda params.CONDA_ENV
    cpus params.N_CPU

    input:
        tuple val(id), path(bam)

    publishDir "${params.PROJECT_PATH}/result/${id}/kallisto", pattern: "*i", mode: "copy"
    
    output:
        path("pseudoalignments.bam.csi")
    
    script:
        """
        samtools index -c -@ ${task.cpus} ${bam}
        """
}

process run_hisat2 {
    conda params.CONDA_ENV
    cpus params.N_CPU

    input:
        tuple val(id), path(r1), path(r2)

    publishDir "${params.PROJECT_PATH}/result/${id}/hisat2", pattern: "*", mode: "copy"

    output:
        tuple val(id), path("aligned.bam"), path("aligned.bam.bai"), path("summary.txt")

    script:
        """
        hisat2 \
        --no-discordant --no-mixed \
        --summary-file summary.txt --mm \
        --new-summary --threads ${task.cpus} \
        --known-splicesite-infile $REFERENCE_PATH/mm10.splice_site.txt\
        -x ${GENOME_INDEX_PATH} \
        -1 ${r1} -2 ${r2} \
        | samtools view -b \
	| samtools sort \
	> aligned.bam \
	&& samtools index aligned.bam
        """    
}


workflow {
    design_file = Channel.fromPath(METADATA_TABLE) 
    fq_ch = design_file \
        .splitCsv(header: true)  \
        .filter{  it.species == "Mouse" } \
        .filter{ it.sample_id != "test"}\
        map{ row -> [row.sample_id, "$FASTQ_PATH/$row.sample_id*R1_001.fastq.gz", "$FASTQ_PATH/$row.sample_id*R2_001.fastq.gz"] } 

    trimmed_fq = cutadapt(fq_ch)
    viral_channel = Channel.of(["viral", "$VIRAL_GENOME"])
    viral_id_list = map_viral(trimmed_fq.combine(viral_channel))
    filtered_viral_fq = filter_viral(viral_id_list.map {row -> [row[0], row[1], row[2], row[3], row[4]]})
    rrna_channel = Channel.of(["rRNA", "$RRNA"])
    rRNA_id_list = map_rRNA(filtered_viral_fq.combine(rrna_channel))
    filtered_rRNA_fq = filter_rRNA(rRNA_id_list.map {row -> [row[0], row[1], row[2], row[3], row[4]]})
    kallisto_ch = run_kallisto(filtered_rRNA_fq) 
    kallisto_index_ch = index_kallisto(kallisto_ch)
    hisat2_ch = run_hisat2(filtered_rRNA_fq)
}

