PROJECT_PATH=/home/wu58/projects/rna-seq/mouse
REF_PATH=${PROJECT_PATH}/ref
GENOME_NAME=mm10
REF_FASTA_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
REF_GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ensGene.gtf.gz"
MM_RRNA_ID=BK000964.3

all: ${REF_PATH}/${GENOME_NAME}.transcripts.kallisto.idx  ${REF_PATH}/viral.fa.1.bt2 ${REF_PATH}/tx2gene.tsv ${REF_PATH}/${GENOME_NAME}.masked.fa.1.ht2 ${REF_PATH}/rRNA.fa.1.bt2 ${REF_PATH}/mm10.splice_site.txt

#HUMAN
${REF_PATH}:
	mkdir -p ${REF_PATH}

${REF_PATH}/genes.gtf: ${REF_PATH}
	curl -o ${REF_PATH}/genes.gtf.gz ${REF_GTF_URL} 
	zcat ${REF_PATH}/genes.gtf.gz > ${REF_PATH}/genes.gtf

${REF_PATH}/tx2gene.tsv: ${REF_PATH}/genes.gtf
	cat ${REF_PATH}/genes.gtf | grep -v "#" | awk '$$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print $$4"\t"$$2}' |   sed 's/\"//g' | sort |uniq > ${REF_PATH}/tx2gene.tsv

${REF_PATH}/${GENOME_NAME}.masked.fa:  ${REF_PATH}
	curl -o ${REF_PATH}/${GENOME_NAME}.masked.fa.gz ${REF_FASTA_URL}
	zcat ${REF_PATH}/${GENOME_NAME}.masked.fa.gz > ${REF_PATH}/${GENOME_NAME}.masked.fa

${REF_PATH}/${GENOME_NAME}.transcripts.fa: ${REF_PATH}/genes.gtf ${REF_PATH}/${GENOME_NAME}.masked.fa
	gffread -g ${REF_PATH}/${GENOME_NAME}.masked.fa ${REF_PATH}/genes.gtf -w ${REF_PATH}/${GENOME_NAME}.transcripts.fa 

${REF_PATH}/${GENOME_NAME}.transcripts.kallisto.idx: ${REF_PATH}/${GENOME_NAME}.transcripts.fa
	kallisto index -i ${REF_PATH}/${GENOME_NAME}.transcripts.kallisto.idx ${REF_PATH}/${GENOME_NAME}.transcripts.fa

${REF_PATH}/viral.fa: ${REF_PATH}
	python dna_to_fasta.py --data-path ${PROJECT_PATH}/data > ${REF_PATH}/viral.fa

${REF_PATH}/viral.fa.1.bt2: ${REF_PATH}/viral.fa
	bowtie2-build ${REF_PATH}/viral.fa ${REF_PATH}/viral.fa

${REF_PATH}/${GENOME_NAME}.masked.fa.1.ht2: ${REF_PATH}/${GENOME_NAME}.masked.fa
	hisat2-build ${REF_PATH}/${GENOME_NAME}.masked.fa ${REF_PATH}/${GENOME_NAME}.masked.fa

${REF_PATH}/rRNA.fa: ${REF_PATH}
	python get_ncbi_fa.py BK000964.3 > ${REF_PATH}/rRNA.fa

${REF_PATH}/rRNA.fa.1.bt2: ${REF_PATH}/rRNA.fa
	bowtie2-build ${REF_PATH}/rRNA.fa  ${REF_PATH}/rRNA.fa

${REF_PATH}/mm10.splice_site.txt: ${REF_PATH}/genes.gtf
	hisat2_extract_splice_sites.py  ${REF_PATH}/genes.gtf > ${REF_PATH}/mm10.splice_site.txt
