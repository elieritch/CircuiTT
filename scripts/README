# =============================================================================
# CTDNA PROSTATE PANEL PIPELINE README
# =============================================================================

This is the ctdna pipeline that I use for the prostate cancer 73 gene panel.
It depends on several tools and environments. Each must be defined in the file config.txt. The config.txt file is sourced in all of the bash scripts for the variables associated with paths to files, binaries and environments.
The three environments that this pipeline depends on are "conda_pysam_env" which includes python 3.7 and pysam, "pysamstats_env" which includes pysamstats (https://github.com/alimanfoo/pysamstats) and python 3.6 (these envs could not be concatenated into a single env as they do have conflicts). The rust environment only includes rust and mattis fasta tools for trimming and masking. If these tools are replaced, or you have a compiled binary, then the rust env can be removed. 

# =============================================================================
# PIPELINE STEPS
# =============================================================================

The pipeline is written for slurm, and is composed of 5 primary steps.

1) Primary processing:
Tumor and normal fastqs are trimmed, masked, aligned, duplicates are removed and read groups are added. This is accomplished with the script "trim_align_removeduplicates.bash". For consistency with other wyatt lab pipelines we are using the fasta rust package for trimming and masking, but trimming can be done with any other tool that uses a trim by quality algorithm such as https://resources.qiagenbioinformatics.com/manuals/clcassemblycell/current/index.php?manual=Quality_trimming.html
For masking, all bases with a base qulaity less than than 20 are replcaed with N.
trim_align_removeduplicates.bash takes as input $read1 $read2 the basename of sample (tumor name or normal name) an outputdir and the config.txt. Note tumor and normal samples cannot have the same name. This step also produces the file *.dup_metrics which is the duplication metrics from picard mark duplicates.

2) Making a tnvstat file:
We make a pileup file which contain statistics of base distributions att each genomic loci. This is done with the script "make_tnvstat_file.bash". This script takes as input tumor and normal bam files produced by the previous script and creates what we call a tnvstat file. A tnvstat file which is a pile up like file. Each row is a base and columns are base level statistics. The name is comes from "tumor-normal variant statistics" and it is pronounced "Tee-eN-Vee-stat" and it has the suffixc ".tnvstat". This script depends on the two environments "conda_pysam_env" and  "pysamstats_env". The tnvstat file is the input for the copynumber and variant calling scripts.

3) Copynumber calling:
Copynumber calling is accomplished using the python script "gene_panel_copynumber.py". It produces 2 plots plots and 2 tsvs. 
The file ${tumor_sample}_gc_correction.png explains the distribution of gc% to log depth ratio before and after lowess correction.
The file ${tumor_sample}_CN_violinplot.pdf is a violin plot of the depth ratio distribution of each gene. This plot is extremely important for determining whether a gene is amplified or deleted. 
The file ${tumor_sample}.cn contains the mean normalized Log2 depth ratio of each gene in the sample under the column header mean.gene.ndr. This is the mean of each violin in the above plot.
The file ${tumor_sample}_b_allele_table.tsv contains each heterozygous loci in the normal sample and the corresponding b allele frequency in the tumor sample.

4) Somatic variant calling:
First a subdirectory called "variant_intermediate_files" is created. It will be filled with several useful files which are often needed to assess pathogenicity of variants or other things.
Somatic snvs are called and annotated using call_snv.py to produce the file ${tumor_sample}_somatic_anno_snv.tsv. Then a filter for exac values, lowcomplexity regions and other things is applied making the file ${tumor_sample}_snv_lcf.tsv. Then a filter is applied to determine if the variant was primarily found in the ends of the reads produceing ${tumor_sample}_snv_lcf_re.tsv which is a common source of false positives. If you need cosmic values, any score or any annotation including the values in the tnvstat, they are all in ${tumor_sample}_snv_lcf_re.tsv.
A similar process occurs with the  the indels with the exception of the read end filter. 
Then the indels and snvs are concatenated and turned into a "${tumor_sample}.maf" file compatible with https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html. The maf file is moved to the primary output directory.

5) Germline varinat calling:
Any germline variant caller will do. We use vardict-java and annotate with annovar. a subdirectory called "vardict" contains intermediate files. I also do tumor-normal somatic variant calling with vardict which stay in the subdirectory and can be used for comparison. But ther germline variniant calls after filtering are moved to the main output directory and are called ${normal_sample}_germline.tsv.


# =============================================================================
# INSTALLATION
# =============================================================================

there are two environment file "pysam_env" and "pysamstats_env". 
these environments need to be created 
conda create -n <NAME> --file pysam_env
conda create -n <NAME> --file pysamstats_env
then put the names you have chosen into the appropriate values in the config.txt

The conda command  will also need to be available within subshells that dont necessarily have users environment variables (ie within parallel), this is why the conda profile profile.d/conda.sh needs to be sourced as well. This script comes with every installation of conda and just exports the conda variables and appends them to the $PATH.

This entire pipleine can be run from the submission node of your cluster using the commands:

script="/PATH/TO/gene_panel_slurm_pipeline.bash";
outputdir="/PATH/TO/DIRECTORY/${name_tumor}";
mkdir -p "/PATH/TO/DIRECTORY/${name_tumor}";
read1_t="/PATH/TO/${name_tumor}_READ_1.fq.gz";
read2_t="/PATH/TO/${name_tumor}_READ_2.fq.gz";
read1_n="/PATH/TO/${name_normal}_READ_1.fq.gz";
read2_n="/PATH/TO/${name_normal}_READ_2.fq.gz";
name_tumor="${name_tumor}";
name_normal="${name_normal}";
configfile="/PATH/TO/config.txt";
sbatch ${script} ${read1_t} ${read2_t} ${read1_n} ${read2_n} ${name_tumor} ${name_normal} ${outputdir} ${configfile};

Mostt often another script is used to create the above code for a given cohort into a file called launch.bash and the that file is executed to submit each slurm job. 

The paths to each of the other things in the config.txt also have to be filled out.
