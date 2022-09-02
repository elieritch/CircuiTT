#!/bin/bash

## this script takes as input tumor and normal bam files produced by the previous script and creates what we call a tnvstat file.
## a tnvstat file which is a pile up like file. Each row is a base and columns are base level statistics
## the name is derived from tumor-normal variant statistics file and it is pronounced "Tee-eN-Vee-stat"
## this script depends on two environments one called pysamstats (python 2.7) which has pysamstats installed and one called samspider (python 3.6) which has pandas and numpy installed.

# =============================================================================
# Set paths and variable
# =============================================================================

printf "total memory: $(free -hm)\n";
printf "starttime: $(date +"%D %H:%M")\n";

bamtumor="$1";
bamnormal="$2";
sample_tumor="$3";
sample_normal="$4";
outputdir="$5";
configfile="$6"

shopt -s extglob;
source ${configfile};
source ${conda_profile_path};
# chrom is an environment variable on all of our machines but redefined here for others
# we dont use chr before our chromosome names because that takes more file space without adding any extra information
chrom=$(printf "$(seq 1 22)\nX\nY\nMT\n");

#redefining these paths because if not present parallel will mask std error
bam_normal=$(readlink -ve $bamnormal);
bam_tumor=$(readlink -ve $bamtumor);
reference=$(readlink -ve $ref);
pysamstats_env=$(echo $pysamstats_env_name)
mkdir -p ${outputdir};
tmp_dir=$(mktemp -d -p ${outputdir});
cd ${tmp_dir};

# =============================================================================
# Make tumor vstat file
# =============================================================================

echo "start running pysamstats on $bam_tumor tumor file at $(date +"%D %H:%M")";

#first do this for each chrom individually and then concatenate the together
echo "$chrom" | parallel --env conda_profile_path --env pysamstats_env --env reference --env bam_tumor --env tmp_dir -j $threads "source $conda_profile_path; conda activate $pysamstats_env; pysamstats --type variation -f $reference --no-dup -c {} $bam_tumor | grep -v ^K | grep -v ^G | grep -v ^M > $tmp_dir/{}.vstat";

#concatenate the chroms while removing bases 
head -n1 $(ls *.vstat | head -n1) > tumor.vstat;
for i in $(echo $chrom); do
echo "concatenating chromosome $i at $(date +"%D %H:%M")";
tail -n +2 $i.vstat | awk -v mind=$min_depth '$4 > mind' >> tumor.vstat;
done;

echo "start adding sample name and linux formatting on $sample_tumor tumor at $(date +"%D %H:%M")";
paste <(head -n1 tumor.vstat | tr -d '\r') <(echo sample) > tumor_tmp.vstat;
cat tumor.vstat | tr -d '\r' | tail -n +2 | awk -F $'\t' -v sample_tumor=$sample_tumor '{print $0,sample_tumor}' FS="\t" OFS="\t" >> tumor_tmp.vstat;

echo "$chrom" | parallel "rm ./{}.vstat";
# =============================================================================
# Make normal vstat file
# =============================================================================

echo "start running pysamstats on $bam_normal normal file at $(date +"%D %H:%M")";

echo "$chrom" | parallel --env conda_profile_path --env pysamstats_env --env reference --env bam_normal --env tmp_dir -j $threads "source $conda_profile_path; conda activate $pysamstats_env; pysamstats --type variation -f $reference --no-dup -c {} $bam_normal | grep -v ^K | grep -v ^G | grep -v ^M > $tmp_dir/{}.vstat";

head -n1 $(ls *.vstat | head -n1) > normal.vstat;
for i in $(echo $chrom); do
echo "concatenating chromosome $i at $(date +"%D %H:%M")";
tail -n +2 $i.vstat | awk -v mind=$min_depth '$4 > mind' >> normal.vstat;
done;

echo "start adding sample name and linux formatting on $sample_normal normal at $(date +"%D %H:%M")";
paste <(cat normal.vstat | tr -d '\r' | head -n1) <(echo sample) > normal_tmp.vstat;
cat normal.vstat | tr -d '\r' | tail -n +2 | awk -F $'\t' -v sample_normal=$sample_normal '{print $0,sample_normal}' FS="\t" OFS="\t" >> normal_tmp.vstat;

echo "$chrom" | parallel "rm ./{}.vstat";
# =============================================================================
# Marry vstat files.
# =============================================================================

conda activate conda_pysam_env;
printf "Marrying tumor and normal at $(date +"%D %H:%M")\n";
tumor=${tmp_dir}/tumor_tmp.vstat;
normal=${tmp_dir}/normal_tmp.vstat;
python ${script_dir}/marry_tumor_normal_vstats.py -t ${tumor} -n ${normal} -o ${tmp_dir}/${sample_tumor}_tmp_tumnormal.tsv;
conda deactivate;

####################################
#add gc and targets annotation
####################################

printf "Finished marrying tumor and normal at $(date +"%D %H:%M")\n";
printf "Adding GC and gene names at $(date +"%D %H:%M")\n";
#printf "chrom\tpos\tref\treads_all_n\tmatches_n\tmismatches_n\tdeletions_n\tinsertions_n\tA_n\tC_n\tT_n\tG_n\tN_n\treads_all_t\tmatches_t\tmismatches_t\tdeletions_t\tinsertions_t\tA_t\tC_t\tT_t\tG_t\tN_t\tTAF_t\tCAF_t\tGAF_t\tAAF_t\tTAF_n\tCAF_n\tGAF_n\tAAF_n\tnormal_max_value\tbase_normal\ttumor_max_value\tbase_tumor\tgc\ttarget\n" > ${tmp_dir}/${sample_tumor}.tnvstat;
printf "chrom\tpos\tref\treads_all_n\tmatches_n\tmismatches_n\tdeletions_n\tinsertions_n\tA_n\tC_n\tT_n\tG_n\tN_n\tsample_n\treads_all_t\tmatches_t\tmismatches_t\tdeletions_t\tinsertions_t\tA_t\tC_t\tT_t\tG_t\tN_t\tsample_t\tTAF_t\tCAF_t\tGAF_t\tAAF_t\tTAF_n\tCAF_n\tGAF_n\tAAF_n\tnormal_max_value\tbase_normal\ttumor_max_value\tbase_tumor\tgc\ttarget\n" > ${tmp_dir}/${sample_tumor}.tnvstat;

paste <(tail -n +2 ${tmp_dir}/${sample_tumor}_tmp_tumnormal.tsv | cut -f 1,2) <(tail -n +2 ${tmp_dir}/${sample_tumor}_tmp_tumnormal.tsv | cut -f 2) <(tail -n +2 ${tmp_dir}/${sample_tumor}_tmp_tumnormal.tsv | tr '\t' ';' ) > ${tmp_dir}/${sample_tumor}_tmp_tumnormal.bed;

${bedtools_path} intersect -a ${tmp_dir}/${sample_tumor}_tmp_tumnormal.bed -b ${gcbed} -wa -wb > ${tmp_dir}/${sample_tumor}.txt;

${bedtools_path} intersect -a ${tmp_dir}/${sample_tumor}.txt -b ${targetbed} -wa -wb -loj | cut -f 4,8,12 | tr ';' '\t' | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' | sort -k1,1V -k2,2n >> ${tmp_dir}/${sample_tumor}.tnvstat;

mv ${tmp_dir}/${sample_tumor}.tnvstat ${outputdir}/${sample_tumor}.tnvstat;
cd ${outputdir};
rm -rf ${tmp_dir};
printf "Finished making tumor normal variant statistic (TNVSTAT) file at $(date +"%D %H:%M")\n";
