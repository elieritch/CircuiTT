#!/bin/bash
#SBATCH --job-name=pipe73gp
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH --cpus-per-task=62
#SBATCH --mem 100000 # memory pool for all cores
#SBATCH -t 05:59:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --workdir=/tmp
#SBATCH --output=/tmp/%j.log
#SBATCH --error=/tmp/%j.log

# =============================================================================
## INPUTS
# =============================================================================

printf "SLURM_JOB_ID=$SLURM_JOB_ID\n";

read1_t="$1"
read2_t="$2"
read1_n="$3"
read2_n="$4"
name_tumor="$5"
name_normal="$6"
outputdir="$7"
configfile="$8"

shopt -s extglob;
source ${configfile};
source ${conda_profile_path};

printf "Arguments:\n";
printf "read1 tumor = $read1_t\n";
printf "read2 tumor = $read2_t\n";
printf "read1 normal = $read1_n\n";
printf "read2 normal = $read2_n\n";
printf "tumor sample = $name_tumor\n";
printf "normal sample = $name_normal\n";
printf "outputdir = $outputdir\n";

printf "Config file contents:\n";
cat ${configfile};

#make sure it exists
mkdir -p ${outputdir};
cd ${outputdir};

# =============================================================================
## MAKE BAM FILES
# =============================================================================

printf "bash ${script_dir}/trim_align_removeduplicates.bash ${read1_t} ${read2_t} ${name_tumor} ${outputdir} ${configfile}\n";
bash ${script_dir}/trim_align_removeduplicates.bash ${read1_t} ${read2_t} ${name_tumor} ${outputdir} ${configfile};
printf "bash ${script_dir}/trim_align_removeduplicates.bash ${read1_n} ${read2_n} ${name_normal} ${outputdir} ${configfile}\n";
bash ${script_dir}/trim_align_removeduplicates.bash ${read1_n} ${read2_n} ${name_normal} ${outputdir} ${configfile};

# =============================================================================
## MAKE TNVSTAT FILES
# =============================================================================
tbam=$(readlink -ve ${outputdir}/${name_tumor}.md.ARRG.bam)
nbam=$(readlink -ve ${outputdir}/${name_normal}.md.ARRG.bam)
printf "bash ${script_dir}/make_tnvstat_file.bash ${tbam} ${nbam} ${name_tumor} ${name_normal} ${outputdir} ${configfile}\n";
bash ${script_dir}/make_tnvstat_file.bash ${tbam} ${nbam} ${name_tumor} ${name_normal} ${outputdir} ${configfile};

# =============================================================================
## COPY NUMBER
# =============================================================================

conda activate ${conda_pysam_env};
tnvstat=$(readlink -ve ${outputdir}/${name_tumor}.tnvstat);
printf "python ${script_dir}/gene_panel_copynumber.py -tnv ${tnvstat} -o ${outputdir}/${name_tumor}\n";
python ${script_dir}/gene_panel_copynumber.py -tnv ${tnvstat} -o ${outputdir}/${name_tumor};

# =============================================================================
## CALL VARIANTS
# =============================================================================

mkdir -p ${outputdir}/variant_intermediate_files;
cd ${outputdir}/variant_intermediate_files;

tnvstat=$(readlink -ve ${outputdir}/${name_tumor}.tnvstat);
error_rate=$(readlink -ve $backgrounderror);
anovar=$(readlink -ve $annovar_path);
output="${outputdir}/variant_intermediate_files/${name_tumor}";
tbam=$(readlink -ve ${outputdir}/${name_tumor}.md.ARRG.bam);

#SNV calling;
python ${script_dir}/call_snv.py -tnv ${tnvstat} -bg ${error_rate} -o ${output} -annovar $anovar;

FILE=${output}_somatic_anno_snv.tsv;
if [[ -f "$FILE" ]]; then python ${script_dir}/lowcomplex_filter_snv.py -snv ${output}_somatic_anno_snv.tsv -ref $ref -o ${output}_snv_lcf.tsv;fi;

FILE=${output}_snv_lcf.tsv;
if [[ -f "$FILE" ]]; then python ${script_dir}/read_end_filter_snv.py -snv ${output}_snv_lcf.tsv -bam $tbam -o ${output}_snv_lcf_re.tsv; fi;

#if there are no snvs then just make an snv header for the later maf file
FILE=${output}_snv_lcf_re.tsv;
if [ ! -f "$FILE" ]; then printf "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.knownGene\tGene.knownGene\tGeneDetail.knownGene\tExonicFunc.knownGene\tAAChange.knownGene\tavsnp147\tExAC_ALL\tcosmic70\tCLINSIG\tCLNDBN\tKaviar_AF\tCADD_raw\tCADD_raw_rankscore\tchrom\tpos\tpos_end\tref\tbase_var\treads_all_n\treads_all_t\tA_n\tA_t\tC_n\tC_t\tT_n\tT_t\tG_n\tG_t\tsample_n\tsample_t\tAAF_n\tAAF_t\tCAF_n\tCAF_t\tTAF_n\tTAF_t\tGAF_n\tGAF_t\tbase_normal\tbase_tumor\tmatches_n\tmismatches_n\tmatches_t\tmismatches_t\tdeletions_n\tinsertions_n\tdeletions_t\tinsertions_t\tn_zygosity\tmean_errorAtoC\tmean_errorAtoG\tmean_errorAtoT\tmean_errorCtoA\tmean_errorCtoG\tmean_errorCtoT\tmean_errorGtoA\tmean_errorGtoC\tmean_errorGtoT\tmean_errorTtoA\tmean_errorTtoC\tmean_errorTtoG\tLSEQ\tRSEQ\tMaskReg10\tlow_com\tmapflag\n" > ${output}_snv_lcf_re.tsv; fi;

#indel calling
python ${script_dir}/call_indel.py -tnv ${tnvstat} -bg ${error_rate} -bam $tbam -o ${output} -annovar $anovar;

FILE=${output}_somatic_anno_indel.tsv;
if [[ -f "$FILE" ]]; then python ${script_dir}/lowcomplex_filter_indel.py -indel ${output}_somatic_anno_indel.tsv -ref $ref -o ${output}_indel_lcf.tsv; fi;

#if there are no snvs then just make an snv header for the later maf file
FILE=${output}_indel_lcf.tsv;
if [ ! -f "$FILE" ]; then printf "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.knownGene\tGene.knownGene\tGeneDetail.knownGene\tExonicFunc.knownGene\tAAChange.knownGene\tavsnp147\tExAC_ALL\tcosmic70\tCLINSIG\tCLNDBN\tKaviar_AF\tCADD_raw\tCADD_raw_rankscore\tchrom\tpos\tpos_end\tref\ttumor_allele\tbase_var\treads_all_t\treads_all_n\tsample_t\tsample_n\ttarget\tmatches_t\tmismatches_t\tA_t\tT_t\tG_t\tC_t\tdeletions_t\tinsertions_t\tmatches_n\tmismatches_n\tA_n\tT_n\tG_n\tC_n\tdeletions_n\tinsertions_n\tDAF_t\tIAF_t\tDAF_n\tIAF_n\tdel_n_zygosity\tins_n_zygosity\tmean_errordel\tmean_errorins\tLSEQ\tRSEQ\tMaskReg10\tlow_com\tmapflag\n" > ${output}_indel_lcf.tsv; fi

# =============================================================================
## MAKE MAF FILE
# =============================================================================
snv=$(readlink -ve ${outputdir}/variant_intermediate_files/${name_tumor}_snv_lcf_re.tsv)
indel=$(readlink -ve ${outputdir}/variant_intermediate_files/${name_tumor}_indel_lcf.tsv)

python ${script_dir}/SNVandINDELtsvs_to_maf.py -snv ${snv} -indel ${indel} -o ${outputdir}/${name_tumor}.maf;

# =============================================================================
## RUN VARDICT FOR GERMLINE AND SOMATIC
# =============================================================================
#NOTE: this vardict-java has xmx changed in its enviro script!!

printf "\nstarting vardict for somatic calling at $datenow\n";
mkdir -p ${outputdir}/vardict;
cd ${outputdir}/vardict;

tbam=$(readlink -ve ${outputdir}/${name_tumor}.md.ARRG.bam);
nbam=$(readlink -ve ${outputdir}/${name_normal}.md.ARRG.bam);

#this vardict-java has xmx changed in its enviro script!!
#take note var2vcf_somatic.pl disappeared in new version!!
$vardict_java_path -G $ref -N $name_tumor -b "$tbam|$nbam" -C -c 1 -S 2 -E 3 -g 4 -th 19 -r 10 -f 0.01 $targetbed | $vardict_script_dir/testsomatic.R | $vardict_script_dir/var2vcf_somatic.pl -C -N "$tbam|$nbam" -f 0.01 > ${outputdir}/vardict/$name_tumor.vcf;
#only the strong survive
cat ${outputdir}/vardict/$name_tumor.vcf | grep -i -e "^#" -e strong > ${outputdir}/vardict/$name_tumor.strong.vcf;
perl $annovar_path ${outputdir}/vardict/$name_tumor.strong.vcf $annovar_humandb -buildver hg38 -out ${outputdir}/vardict/$name_tumor -remove -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel -operation g,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput;
echo "finished somatic vardict on $name_tumor";

printf "\nstarting vardict for germline calling at $(date +"%D %H:%M")\n";
#now germline calling:
vardict-java -G $ref -N $name_normal -b $nbam -C -c 1 -S 2 -E 3 -g 4 -th 15 -r 10 -f 0.1 $targetbed | $vardict_script_dir/teststrandbias.R | $vardict_script_dir/var2vcf_valid.pl -C -N $nbam > ${outputdir}/vardict/$name_normal.germ.vcf;
#only the pass move on
cat ${outputdir}/vardict/$name_normal.germ.vcf | grep -i -e ^# -e PASS > ${outputdir}/vardict/$name_normal.pass.germ.vcf;
echo "vardict finished, now running annovar for $name_normal";
perl $annovar_path ${outputdir}/vardict/$name_normal.pass.germ.vcf $annovar_humandb -buildver hg38 -out ${outputdir}/vardict/${name_normal}_germ -remove -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel -operation g,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput;
echo "finished germ vardict on $name_normal";

printf "SAMPLE\tTYPE\tDP\tEND\tVD\tAF\tBIAS\tREFBIAS\tVARBIAS\tPMEAN\tPSTD\tQUAL\tQSTD\tSBF\tODDRATIO\tMQ\tSN\tHIAF\tADJAF\tSHIFT3\tMSI\tMSILEN\tNM\tHICNT\tHICOV\tLSEQ\tRSEQ\n" > ${outputdir}/vardict/${name_normal}_germ.tmp;
cat ${outputdir}/vardict/${name_normal}_germ.hg38_multianno.txt | cut -f 129 | sed "s/$/#/g" | tr ';' '\n' | tr '=' '\t' | cut -f 2 | tr '\n' '\t' | tr '#' '\n' | sed 's/^\t//g' >> ${outputdir}/vardict/${name_normal}_germ.tmp;
paste <(cat ${outputdir}/vardict/${name_normal}_germ.hg38_multianno.txt | tr -s ' ' '_' | cut -f 1-119) <(cat ${outputdir}/vardict/${name_normal}_germ.tmp) | tr -s ' ' '_' | grep -e exonic -e ^Chr | awk '$6 == "exonic" || $6 == "Func.refGene"' | awk '($140 < 3) || ($140 == "MSI")' | awk '(($125 > 0.3) && ($125 < 0.7)) || ($125 == "AF")' | grep -v "BIVM-ERCC5;ERCC5" | awk '($17 < 0.05) || ($17 == ".") || ($17 == "ExAC_ALL")' > ${outputdir}/vardict/${name_normal}_germline.tsv;

cp ${outputdir}/vardict/${name_normal}_germline.tsv ${outputdir}/${name_normal}_germline.tsv;

printf "\n\n";
printf "copying log file at $(date +"%D %H:%M")\n";
mv /tmp/$SLURM_JOB_ID.log ${outputdir}/$SLURM_JOB_ID.log;
printf "pipeline finished at $(date +"%D %H:%M")\n";

