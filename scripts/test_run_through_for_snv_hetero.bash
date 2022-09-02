script="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/gene_panel_slurm_pipeline.bash";
mkdir -p /groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero;
read1_t="/groups/wyattgrp/data/fq/GU_73gp/GU-17-150-07Nov2017-cfDNA_1.fq.gz";
read2_t="/groups/wyattgrp/data/fq/GU_73gp/GU-17-150-07Nov2017-cfDNA_2.fq.gz";
read1_n="/groups/wyattgrp/data/fq/GU_73gp/GU-17-150-WBC_1.fq.gz";
read2_n="/groups/wyattgrp/data/fq/GU_73gp/GU-17-150-WBC_2.fq.gz";
name_tumor="GU-17-150-07Nov2017-cfDNA";
name_normal="GU-17-150-WBC";
outputdir="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero/${name_tumor}";
configfile="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/config2.txt";
sbatch ${script} ${read1_t} ${read2_t} ${read1_n} ${read2_n} ${name_tumor} ${name_normal} ${outputdir} ${configfile};

#sample_7 --> 17-154
script="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/gene_panel_slurm_pipeline.bash";
mkdir -p /groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero;
read1_t="/groups/wyattgrp/data/fq/GU_73gp/GU-17-154-cfDNA-UMI-2018May02_1.fq.gz";
read2_t="/groups/wyattgrp/data/fq/GU_73gp/GU-17-154-cfDNA-UMI-2018May02_2.fq.gz";
read1_n="/groups/wyattgrp/data/fq/GU_73gp/GU-17-154-WBC_1.fq.gz";
read2_n="/groups/wyattgrp/data/fq/GU_73gp/GU-17-154-WBC_2.fq.gz";
name_tumor="GU-17-154-cfDNA-UMI-2018May02";
name_normal="GU-17-154-WBC";
outputdir="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero/${name_tumor}";
configfile="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/config2.txt";
sbatch ${script} ${read1_t} ${read2_t} ${read1_n} ${read2_n} ${name_tumor} ${name_normal} ${outputdir} ${configfile};

#sample_23 --> 17-170
script="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/gene_panel_slurm_pipeline.bash";
mkdir -p /groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero;
read1_t="/groups/wyattgrp/data/fq/M1_73gp/GU-17-170-cfDNA-2017Jul12_1.fq.gz";
read2_t="/groups/wyattgrp/data/fq/M1_73gp/GU-17-170-cfDNA-2017Jul12_2.fq.gz";
read1_n="/groups/wyattgrp/data/fq/M1_73gp/17-170-WBC_1.fq.gz";
read2_n="/groups/wyattgrp/data/fq/M1_73gp/17-170-WBC_1.fq.gz";
name_tumor="GU-17-170-cfDNA-2017Jul12";
name_normal="GU-17-170-WBC";
outputdir="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero/${name_tumor}";
configfile="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/config2.txt";
sbatch ${script} ${read1_t} ${read2_t} ${read1_n} ${read2_n} ${name_tumor} ${name_normal} ${outputdir} ${configfile};

#Wyatt_gDNA_sample_38 --> 20-043	
script="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/gene_panel_slurm_pipeline.bash";
mkdir -p /groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero;
read1_t="/groups/wyattgrp/data/fq/GU_73gp/GU-20-043-cfDNA-Baseline-2020Jan14_1.fq.gz";
read2_t="/groups/wyattgrp/data/fq/GU_73gp/GU-20-043-cfDNA-Baseline-2020Jan14_2.fq.gz";
read1_n="/groups/wyattgrp/data/fq/GU_73gp/GU-20-043-WBC-2020Jan14_1.fq.gz";
read2_n="/groups/wyattgrp/data/fq/GU_73gp/GU-20-043-WBC-2020Jan14_2.fq.gz";
name_tumor="GU-20-043-cfDNA-Baseline-2020Jan14";
name_normal="GU-20-043-WBC-2020Jan14";
outputdir="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/test_data_output_forhetero/${name_tumor}";
configfile="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/scripts/config2.txt";
sbatch ${script} ${read1_t} ${read2_t} ${read1_n} ${read2_n} ${name_tumor} ${name_normal} ${outputdir} ${configfile};
