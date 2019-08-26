#!/bin/bash
#SBATCH --job-name=RNASeq_pipe2
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=50G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o ../log_files/%x_%A.out
#SBATCH -e ../log_files/%x_%A.err


d1="raw_data"
d2="raw_data_combined"
d3="fastqc"
d4="quality_control"
d5="sickle" #can have different folder names inside d4 (quality_control)
d6="after_sickle"
d7="mapping"
d8="hisat_sickle" #can have different folder names inside d7 (mapping)
d9="remove_duplicates"
d10="counts"

R1="_R1.fastq"

index_path="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/HISAT2/Rattus_norvegicus"
splice_site_path="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/splice_sites"
GFF_File="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gff3"
GTF_File="/isg/shared/databases/alignerIndex/animal/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gtf"

echo `hostname`
echo ${1}

##################################################################
## 			Unzip Data				##
##################################################################
echo "start : `date`"

for f in ../raw_data/${1}/*.gz; do
	gunzip -d $f
done

#gunzip -d ../raw_data/${1}/*.gz

echo "end : `date`"
echo "=========== Unzip Done ================="

##################################################################
##                      Raw Data Combine                        ##
##################################################################
if [ ! -d ../${d2} ]; then
        mkdir -p ../${d2}
fi

cd ../${d2}
#create sub directories
if [ ! -d ${1} ]; then
        mkdir -p ${1}
fi

echo "start ${d2} : `date`"

cat `ls ../${d1}/${1}/*_R1_* | sort -` >> ./${1}/${1}${R1}

echo "end ${d2} : `date`"
echo "=========== Combine Done =============="

##################################################################
## Quality Check of the Reads                                   ##
##################################################################
module load fastqc/0.11.5

if [ ! -d ../${d3} ]; then
       mkdir -p ../${d3}
fi

cd ../${d3}/

if [ ! -d ${d1} ]; then
        mkdir -p ${d1}
fi

echo "=========== fastqc before Starts :: $1 ================="
echo "start ${d3} : `date`"

fastqc --outdir ./${d1}/ ../${d2}/${1}/${1}${R1}
#fastqc --outdir ./${d1}/ ../${d2}/${1}/${1}${R2}
module unload fastqc/0.11.5
echo "end ${d3} : `date`"
echo "=========== fastqc before Done :: $1 ================="


##################################################################
## Trimming of the Reads                                        ##
##################################################################
if [ ! -d ../${d4} ]; then
        mkdir -p ../${d4}
fi

cd ../${d4}/

if [ ! -d ${d5} ]; then
        mkdir -p ${d5}
fi

echo "=========== sickle starting :: $1 ================="
module load sickle/1.33

echo "start ${d4} : `date`"

sickle se \
	-f ../${d2}/${1}/${1}${R1} \
	-t sanger \
	-o ./${d5}/trimmed_${1}${R1} \
	-q 30 \
	-l 45 -n

echo "end ${d4} : `date`"
module unload sickle/1.33
echo "=========== sickle Done :: $1 ================="


##################################################################
## FASTQC after Trimming                                        ##
##################################################################
module load fastqc/0.11.5
cd ../fastqc/

if [ ! -d ${d6} ]; then
        mkdir -p ${d6}
fi

echo "start ${d6} FASTQC : `date`"
fastqc --outdir ./${d6}/ ../${d4}/${d5}/trimmed_${1}${R1}
#fastqc --outdir ./${d6}/ ../${d4}/${d5}/trimmed_${1}${R2}

echo "end ${d6} FASTQC : `date`"
module unload fastqc/0.11.5
echo "=========== After Trimming fastqc Done :: $1 ================="



##################################################################
## Mapping                                                      ##
##################################################################
echo "=========== mapping starting :: $1 ================="
module load hisat2/2.1.0
if [ ! -d ../${d7}/${d8} ]; then
        mkdir -p ../${d7}/${d8}
fi

cd ../${d7}
echo ${1} >> ../log_files/${d8}_mapping_%x_%A.err

#index_path="../../index/HISAT2/Oryctolagus_cuniculus"
echo "start ${d7} : `date`"
hisat2 -p 8 -x ${index_path} \
	--known-splicesite-infile ${splice_site_path} \
        -U ../${d4}/${d5}/trimmed_${1}${R1} \
        -S ./${d8}/${1}_mapped.sam

echo "end ${d7} : `date`"
module unload hisat2/2.0.5
echo "\n\n" >> ../log_files/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err
cat ../log_files/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.err >> ${d8}_${1}_mapping_rate.log
echo "=========== mapping Done :: $1 =================\n\n"


##################################################################
## SAM to BAM
##################################################################
echo "=========== sam to bam :: $1 ================="
module load samtools/1.7
echo "start SAMtools : `date`"
samtools view -@ 8 -bhS ./${d8}/${1}_mapped.sam -o ./${d8}/${1}_mapped.bam
samtools sort -@ 8 ./${d8}/${1}_mapped.bam -o ./${d8}/${1}_mapped_sort.bam
echo "end SAMTOOLS : `date`"
module unload samtools/1.7
echo "=========== sam to bam Done :: $1 =================\n\n"


##################################################################
## Remove duplicates
##################################################################
echo "=========== Removing Duplicates :: $1 ================="
echo "start ${d9} : `date`"
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

if [ ! -d ../${d9}/${d8} ]; then
        mkdir -p ../${d9}/${d8}
fi

cd ../${d9}/${d8}

java -jar $PICARD MarkDuplicates \
	INPUT=../../${d7}/${d8}/${1}_mapped_sort.bam \
        OUTPUT=./${1}_dup_removed.bam \
	METRICS_FILE=./${1}_dup_removed_metrics \
	REMOVE_DUPLICATES= true \
	CREATE_INDEX= true

module unload picard/2.2.1
echo "end ${d9} : `date`"
echo "=========== Removing Duplicates Done =================\n\n"


##################################################################
## Sorting the BAM reads                                        ##
##################################################################
module load samtools/1.7
echo "=========== bam sort  :: $1 ================="
echo "start : `date`"
samtools sort -@ 8 ./${1}_dup_removed.bam -o ./${1}_dup_removed_sort.bam

module unload samtools/1.7
echo "end : `date`"
echo "=========== bam sort Done :: $1 =================\n\n"

cd ../

##################################################################
## Generate Counts                                              ##
##################################################################
echo "=========== htseq counts :: $1 ================="
echo "start : `date`"

module load htseq/0.9.1
if [ ! -d ../${d10}/${d8} ]; then
        mkdir -p ../${d10}/${d8}
fi

cd ../${d10}/${d8}
echo ${1}
htseq-count -s no -r pos -t exon -i gene_id --additional-attr=gene_name -f bam ../../${d9}/${d8}/${1}_dup_removed_sort.bam ${GTF_File} > ${1}.counts

echo "end : `date`"
module unload htseq/0.9.1
echo "=========== Htseq Count Done :: "$1" ================="


