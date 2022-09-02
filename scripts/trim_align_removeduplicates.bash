#!/bin/bash

# =============================================================================
# Set paths and variable for allingment and primary processing
# =============================================================================

printf "total memory: $(free -hm)\n";
printf "starttime: $(date +"%D %H:%M")\n";

read1="$1"
read2="$2"
basen="$3"
outputdir="$4"
configfile="$5"

shopt -s extglob;
source ${configfile};
source ${conda_profile_path};

mkdir -p ${outputdir};
tmp_dir=$(mktemp -d -p ${outputdir});
#tmp_dir should be the full path, but its failed before, so just incase get full path
tmpdir=$(readlink -ve ${tmp_dir});

printf "Arguments:\n";
printf "read1 = $read1\n";
printf "read2 = $read2\n";
printf "base name = $basen\n";
printf "outputdir = $outputdir\n";
printf "reference = $ref\n";
printf "number of threads = $threads\n";

# =============================================================================
# Trimming and masking
# For consistency with other wyatt lab pipelines we are using the fasta rust package
# but trimming can be done with any other tool that uses a trim by quality algorithm such as:
# https://resources.qiagenbioinformatics.com/manuals/clcassemblycell/current/index.php?manual=Quality_trimming.html
# =============================================================================
conda activate $rust_env_name;

printf "\nstart trimming reads from $basen at $(date +"%D %H:%M")\n";
fasta trim by quality ${read1} 30 > ${tmpdir}/${basen}_1.trim.fq &
fasta trim by quality ${read2} 30 > ${tmpdir}/${basen}_2.trim.fq &
wait

printf "\nstart masking reads from $basen at $(date +"%D %H:%M")\n";
fasta mask by quality ${tmpdir}/${basen}_1.trim.fq 20 > ${tmpdir}/${basen}_1.trim.mask.fq &
fasta mask by quality ${tmpdir}/${basen}_2.trim.fq 20 > ${tmpdir}/${basen}_2.trim.mask.fq &
wait

printf "\nfinished trimming and masking reads from $basen at $(date +"%D %H:%M")\n";
conda deactivate;
# =============================================================================
# Align and sort
# =============================================================================

subthreads=$(echo "$threads / 3" | bc);
printf "\nstart bwa mem alignment of reads from $basen at $(date +"%D %H:%M")\n";
printf "Command:\n";
printf "${bwa_path} mem -M -t $subthreads $ref ${tmpdir}/${basen}_1.trim.mask.fq ${tmpdir}/${basen}_2.trim.mask.fq | ${samtools_path} view -u -q 20 -@ $subthreads | ${samtools_path} sort -@ $subthreads > ${tmpdir}/${basen}.bam\n";

${bwa_path} mem -M -t $subthreads $ref ${tmpdir}/${basen}_1.trim.mask.fq ${tmpdir}/${basen}_2.trim.mask.fq | ${samtools_path} view -u -q 20 -@ $subthreads | ${samtools_path} sort -@ $subthreads > ${tmpdir}/${basen}.bam;
${samtools_path} index ${tmpdir}/${basen}.bam;

printf "finished bwa mem alignment of reads from $basen at$(date +"%D %H:%M")\n";

# =============================================================================
# run picard markduplicate and reindex
# through trial and error Xmx40g was found to be good for marking and removing duplicates
# this lets us have 5 jobs per node
# =============================================================================

printf "\nstarting picard mark duplicates at $(date +"%D %H:%M")\n";
printf "Command:\n";
printf "java -jar -Xmx40g ${picard_path} MarkDuplicates INPUT=${tmpdir}/${basen}.bam OUTPUT=${tmpdir}/${basen}.md.bam METRICS_FILE=${tmpdir}/${basen}.dup_metrics REMOVE_DUPLICATES=TRUE TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT;\n";

java -jar -Xmx40g ${picard_path} MarkDuplicates INPUT=${tmpdir}/${basen}.bam OUTPUT=${tmpdir}/${basen}.md.bam METRICS_FILE=${tmpdir}/${basen}.dup_metrics REMOVE_DUPLICATES=TRUE TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT;
${samtools_path} index ${tmpdir}/${basen}.md.bam;

printf "finished picard mark duplicates at $(date +"%D %H:%M")\n";

# =============================================================================
# run picard AddOrReplaceReadGroups and reindex
# This is necessary for some bam file formatting for other downstream tools we use
# that check whether read groups are assigned such as any gatk tools
# =============================================================================

printf "\nstarting picard AddOrReplaceReadGroups at $(date +"%D %H:%M")\n";
printf "Command:\n";
printf "java -XX:ParallelGCThreads=${threads} -jar ${picard_path} AddOrReplaceReadGroups INPUT=${tmpdir}/${basen}.md.bam OUTPUT=${tmpdir}/${basen}.md.ARRG.bam RGLB=${basen} RGPL=illumina RGPU=${basen} RGSM=${basen} TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT\n";

java -XX:ParallelGCThreads=${threads} -jar ${picard_path} AddOrReplaceReadGroups INPUT=${tmpdir}/${basen}.md.bam OUTPUT=${tmpdir}/${basen}.md.ARRG.bam RGLB=${basen} RGPL=illumina RGPU=${basen} RGSM=${basen} TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT;
${samtools_path} index ${tmpdir}/${basen}.md.ARRG.bam;

printf "finished picard AddOrReplaceReadGroups at $(date +"%D %H:%M")\n";

# =============================================================================
# move relevant files from tmp to outputdir
# =============================================================================

mv ${tmpdir}/${basen}.dup_metrics ${outputdir}/${basen}.dup_metrics;
mv ${tmpdir}/${basen}.md.ARRG.bam ${outputdir}/${basen}.md.ARRG.bam;
mv ${tmpdir}/${basen}.md.ARRG.bam.bai ${outputdir}/${basen}.md.ARRG.bam.bai;
rm -r ${tmpdir};

printf "finished allignment processing at $(date +"%D %H:%M")\n";
printf "files created:\n";
readlink -ve ${outputdir}/${basen}.dup_metrics;
readlink -ve ${outputdir}/${basen}.md.ARRG.bam;
readlink -ve  ${outputdir}/${basen}.md.ARRG.bam.bai;
