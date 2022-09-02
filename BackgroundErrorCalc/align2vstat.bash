#!/bin/bash
#SBATCH --job-name=fq2vstat
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH -n 30 # number of cores
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 99:59:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --workdir=/groups/wyattgrp/log
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log

fr=$(free -hm);
echo "total memory:";
echo "$fr";
datenow=$(date +"%D %H:%M");
printf "starttime: $datenow\n"

## need for cp pipe
read1="$1"
read2="$2"
basen="$3"
outputdir="$4"

ref="/groups/wyattgrp/reference/grch38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa";

#targetbed="/groups/wyattgrp/reference/target/72gene_panel/grch38genes4intersect.bed";
#gcbed="/groups/wyattgrp/reference/grch38/grch38.gc50Base.bed";
#mindepth=40
#cores=

# printf "Arguments:\n";
# printf "read1 tumor = $read1_t\n";
# printf "read2 tumor = $read2_t\n";
# printf "read1 normal = $read1_n\n";
# printf "read2 normal = $read2_n\n";
# printf "basen = $basen\n";
# printf "outputdir = $outputdir\n";
# printf "reference = $ref\n";
# printf "background snv = $bgsnv\n";
# printf "background indel = $bgindel\n";

#source activate pysamstats;
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh";
source ${conda_profile_path};
conda activate pysamstats;

shopt -s extglob;
mkdir -p $outputdir/bam;
mkdir -p $outputdir/tmp;
mkdir -p $outputdir/bam/tmp;
###########################
#tumor alignment to vstat #
###########################

#run bwa samtools etc...
printf "\nstart running bwa on $basen tumor at $(date +"%D %H:%M")\n";

bwa mem -M -t 20 $ref $read1 $read2 | samtools view -u -q 10 | samtools sort -@ 10 > $outputdir/bam/$basen.bam;
samtools index $outputdir/bam/$basen.bam;

java -jar -Xmx40g /home/eritch/anaconda3/pkgs/picard-2.18.0-py36_0/share/picard-2.18.0-0/picard.jar MarkDuplicates INPUT=$outputdir/bam/$basen.bam OUTPUT=$outputdir/bam/$basen.md.bam METRICS_FILE=$outputdir/bam/$basen.output.dup_metrics REMOVE_DUPLICATES=TRUE TMP_DIR=$outputdir/tmp VALIDATION_STRINGENCY=LENIENT;
samtools index $outputdir/bam/$basen.md.bam;

java -XX:ParallelGCThreads=30 -jar /home/eritch/anaconda3/pkgs/picard-2.18.0-py36_0/share/picard-2.18.0-0/picard.jar AddOrReplaceReadGroups INPUT=$outputdir/bam/$basen.md.bam OUTPUT=$outputdir/bam/$basen.md.ARRG.bam RGLB=$basen RGPL=illumina RGPU=$basen RGSM=$basen TMP_DIR=$outputdir/tmp VALIDATION_STRINGENCY=LENIENT;
samtools index $outputdir/bam/$basen.md.ARRG.bam;

#run pysamstats

pysamstats --type variation -f $ref --no-dup --min-baseq=20  $outputdir/bam/$basen.md.ARRG.bam | grep -v ^KI270 | grep -v ^GL0002 > $outputdir/bam/${basen}.tmp.vstat;

#format pysamstats
echo "start adding sample name and linux formatting on $basen tumor at $(date +"%D %H:%M")";
paste <(cat $outputdir/bam/${basen}.tmp.vstat | tr -d '\r' | head -n1) <(echo sample) > $outputdir/bam/${basen}.vstat;
cat $outputdir/bam/${basen}.tmp.vstat | tr -d '\r' | tail -n +2 | awk -F $'\t' -v basen=$basen '{print $0,basen}' FS="\t" OFS="\t" >> $outputdir/bam/${basen}.vstat;

echo "finished adding sample name and linux formatting on $basen tumor at $(date +"%D %H:%M")";
