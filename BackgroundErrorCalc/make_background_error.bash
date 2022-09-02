#!/bin/bash

#SBATCH --job-name=ErrDist
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH -n 20 # number of cores
#SBATCH --mem 720000 # memory pool for all cores
#SBATCH -t 23:59:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --workdir=/groups/wyattgrp/log
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log
conda_profile_path="/home/$(whoami)/anaconda3/etc/profile.d/conda.sh";
source ${conda_profile_path};
conda activate samspider;

#source activate pysam;
shopt -s extglob;

#exome should recieve unlimited time and max ram;
#to do: 
# make compute mean parallel;
# change the grep chromosome to something faster;
# delete the directories after done, saving just $outputdir/bg_dir/chr${chromosome}_errorrates.tsv
# right now deleting toi make sure split doesnt make binary files
# add custom logging


vstat_dir="$1"
outputdir="$2"
chromosome="$3"

mkdir -p $outputdir;
mkdir -p $outputdir/chrom_dir;
mkdir -p $outputdir/bg_dir;
mkdir -p $outputdir/chrom_split;

fr=$(free -hm);
echo "total memory:";
echo "$fr";
datenow=$(date +"%D %H:%M");
printf "starttime: $datenow\n";

########################
# subset by chromosome #
########################

cd $vstat_dir;
#tried doin this with parallel and some rows always end up getting fucked up;
for i in $(ls *.vstat | grep -v tmp)
do
	datenow=$(date +"%D %H:%M");
	printf "$datenow\n";
	printf "grep ^$chromosome[[:space:]] $i | cut -f 1-4,6,8,10,12,14,16,18,20,22,24 | awk '(\$4 >= 1)' >> $outputdir/chrom_dir/${chromosome}.tsv\n";
	grep ^$chromosome[[:space:]] $i | cut -f 1-4,6,8,10,12,14,16,18,20,22,24 | awk '($4 >= 1)' >> $outputdir/chrom_dir/${chromosome}.tsv;
done;

##################
# sort and split #
##################

datenow=$(date +"%D %H:%M");
printf "sorting: $datenow\n";
printf "sort --temporary-directory=$vstat_dir -k2,2n $outputdir/chrom_dir/${chromosome}.tsv > $outputdir/chrom_dir/${chromosome}_sorted.tsv\n";
sort --temporary-directory=$vstat_dir -k2,2n $outputdir/chrom_dir/${chromosome}.tsv > $outputdir/chrom_dir/${chromosome}_sorted.tsv;

datenow=$(date +"%D %H:%M");
printf "splitting $datenow\n";
#printf "split -d -l 1000000 --verbose $outputdir/chrom_dir/${chromosome}_sorted.tsv /groups/wyattgrp/eritch/projects/make_poisson_wxs/chrsplitnew/$prefix\n";
#printf "split -d 1000 -l 10000000 --verbose $outputdir/chrom_dir/${chromosome}_sorted.tsv $outputdir/chrom_split/${chromosome}_sort\n";
#split -d 1000 -l 10000000 --verbose $outputdir/chrom_dir/${chromosome}_sorted.tsv $outputdir/chrom_split/${chromosome}_sort;

#When i try split to 10000000 files or bigger, non deterministic creation of small binary files
printf "split --suffix-length=4 --numeric-suffixes=1000 --lines=1000000 --verbose $outputdir/chrom_dir/${chromosome}_sorted.tsv $outputdir/chrom_split/${chromosome}_sort\n";
split --suffix-length=4 --numeric-suffixes=1000 --lines=1000000 --verbose $outputdir/chrom_dir/${chromosome}_sorted.tsv $outputdir/chrom_split/${chromosome}_sort;

datenow=$(date +"%D %H:%M");
printf "finished split and sort $datenow\n";

######################################################
# Make sure that positions arent split between files #
######################################################

cd $outputdir/chrom_split/;

for j in $(ls | grep ^${chromosome}_sort | sort -n)
do
	number=$(echo $j | sed "s/${chromosome}_sort//g");
	number_plus1=$(echo "$number + 1" | bc -l);
	onenottwo=$(echo ${#number_plus1});
	
	if [ "$onenottwo" -eq "1" ]; then
		fileplusone=${chromosome}_sort0${number_plus1};
	else
		fileplusone=${chromosome}_sort${number_plus1};
	fi;
	
	if [ -f $fileplusone ]; then
		datenow=$(date +"%D %H:%M");
		printf "\nchromosome=$chromosome \nj=$j \nnumber=$number \nnumber_plus1=$number_plus1 \nonenottwo=$onenottwo \nfileplusone=$fileplusone\n";
		value_needs_adding=$(tail -n1 $j | cut -f 2);
		printf "value_needs_adding = $value_needs_adding\n";
		printf "head -n 100 $fileplusone | grep ^$chromosome[[:space:]]$value_needs_adding >> $j\n";
		head -n 100 $fileplusone | grep ^$chromosome[[:space:]]$value_needs_adding >> $j;
		#printf "sed -i.bak '/^$chromosome[[:space:]]$value_needs_adding[[:space:]]/d' $fileplusone\n";
		#sed -i.bak '/^$chromosome[[:space:]]$value_needs_adding[[:space:]]/d' $fileplusone;
		printf "grep -v ^$chromosome[[:space:]]$value_needs_adding[[:space:]] $fileplusone > ${fileplusone}_fix\n";
		grep -v ^$chromosome[[:space:]]$value_needs_adding[[:space:]] $fileplusone > ${fileplusone}_fix;
		rm $fileplusone;
		mv ${fileplusone}_fix $fileplusone;
	else
		printf "$fileplusone does not exist; $j is the last file for $chromosome\n";
	fi;
done;

datenow=$(date +"%D %H:%M");
printf "finished repair after split $datenow\n";

##################
# compute lambda #
##################

for i in $(ls | grep ^${chromosome}_sort)
do
	input=$(readlink -f $i);
	output="$outputdir/bg_dir/$i";
	python /groups/wyattgrp/scripts/make_background_error2/make_snv_error.py -in $input -out $output;
done;

datenow=$(date +"%D %H:%M");
printf "\nfinished computing means $datenow\n";

cd $outputdir/bg_dir;

###################
# remove all zero #
###################
#faster than doing in python

printf "chrom\tpos\tref\tmean_errorAtoC\tmean_errorAtoT\tmean_errorAtoG\tmean_errorTtoC\tmean_errorTtoG\tmean_errorTtoA\tmean_errorGtoC\tmean_errorGtoT\tmean_errorGtoA\tmean_errorCtoG\tmean_errorCtoT\tmean_errorCtoA\tmean_errordel\tmean_errorins\n" > $outputdir/bg_dir/chr${chromosome}_errorrates.tsv;
for i in $(ls *_chromosome_mean_errors.tsv | grep ^${chromosome}_sort | sort -n);do
grep -v "[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0" $i >> $outputdir/bg_dir/chr${chromosome}_errorrates.tsv;
done;


#########################################################################################################################
# 
# call with something like this:
# sbatch /groups/wyattgrp/eritch/projects/make_poisson_wxs/scripts2/chrom_level_error4.bash
# /groups/wyattgrp/eritch/projects/make_poisson_wxs/bam /groups/wyattgrp/eritch/projects/make_poisson_wxs/newtest 1
#
#########################################################################################################################


# mkdir -p /groups/wyattgrp/eritch/projects/make_poisson_73gp/new_chrom_level4;
# for i in $(echo "$chrom"); do 
# vstat_dir="/groups/wyattgrp/eritch/projects/make_poisson_73gp/bam";
# outputdir="/groups/wyattgrp/eritch/projects/make_poisson_73gp/new_chrom_level4";
# chromosome=$i;
# printf "sbatch /groups/wyattgrp/eritch/projects/make_poisson_wxs/scripts2/chrom_level_error4.bash $vstat_dir $outputdir $chromosome\n";
# sbatch /groups/wyattgrp/eritch/projects/make_poisson_wxs/scripts2/chrom_level_error4.bash $vstat_dir $outputdir $chromosome;
# done;


# mkdir -p /groups/wyattgrp/eritch/projects/make_poisson_wxs/new_chrom_level4;
# for i in $(echo "$chrom"); do 
# vstat_dir="/groups/wyattgrp/eritch/projects/make_poisson_wxs/bam";
# outputdir="/groups/wyattgrp/eritch/projects/make_poisson_wxs/new_chrom_level4";
# chromosome=$i;
# printf "sbatch /groups/wyattgrp/eritch/projects/make_poisson_wxs/scripts2/chrom_level_error4_wxs.bash $vstat_dir $outputdir $chromosome\n";
# sbatch /groups/wyattgrp/eritch/projects/make_poisson_wxs/scripts2/chrom_level_error4_wxs.bash $vstat_dir $outputdir $chromosome;
# done;

#########################################################
#	i did this part on interactive session:		#
#########################################################

# mkdir /groups/wyattgrp/reference/bg/73gp_snv_18may18;
# mkdir /groups/wyattgrp/reference/bg/wxs_snv_18may18;
# mkdir /groups/wyattgrp/reference/bg/73gp_indel_18may18;
# mkdir /groups/wyattgrp/reference/bg/wxs_indel_18may18;

# for i in $(echo "$chrom"); do 
# sourcedir="/groups/wyattgrp/eritch/projects/make_poisson_73gp/new_chrom_level4/bg_dir";
# foo=$(readlink -ve $sourcedir/chr${i}_errorrates.tsv);
# targetdir="/groups/wyattgrp/reference/bg/73gp_snv_18may18";
# bar=$(echo $targetdir/chr${i}_errorrates.tsv);
# printf "cp $foo $bar\n";
# cp -v $foo $bar;
# done;

# for i in $(echo "$chrom"); do 
# sourcedir="/groups/wyattgrp/eritch/projects/make_poisson_wxs/new_chrom_level4/bg_dir";
# foo=$(readlink -ve $sourcedir/chr${i}_errorrates.tsv);
# targetdir="/groups/wyattgrp/reference/bg/wxs_snv_18may18";
# bar=$(echo $targetdir/chr${i}_errorrates.tsv);
# printf "cp $foo $bar\n";
# cp -v $foo $bar;
# done;

# cd /groups/wyattgrp/reference/bg/73gp_snv_18may18;
# head -n1 $(ls *errorrates.tsv | head -n1) > allchr_errorrates.txt;
# for i in $(echo "$chrom");do 
# printf "tail -n +2 chr${i}_errorrates.tsv >> allchr_errorrates.txt\n";
# tail -n +2 chr${i}_errorrates.tsv >> allchr_errorrates.txt;
# done;

# cd /groups/wyattgrp/reference/bg/wxs_snv_18may18;
# head -n1 $(ls *errorrates.tsv | head -n1) > allchr_errorrates.txt;
# for i in $(echo "$chrom");do 
# printf "tail -n +2 chr${i}_errorrates.tsv >> allchr_errorrates.txt\n";
# tail -n +2 chr${i}_errorrates.tsv >> allchr_errorrates.txt;
# done;

# mkdir /groups/wyattgrp/reference/bg/73gp_snv_18may18/chr;
# cd /groups/wyattgrp/reference/bg/73gp_snv_18may18;
# ls | grep ^chr | prl "mv {} /groups/wyattgrp/reference/bg/73gp_snv_18may18/chr/{}";

# cat /groups/wyattgrp/reference/bg/73gp_snv_18may18/allchr_errorrates.txt | cut -f 1-3,16,17 | grep -v "[[:space:]]0\.0[[:space:]]0\.0" > /groups/wyattgrp/reference/bg/73gp_indel_18may18/allchr_errorrates.txt;

# cd /groups/wyattgrp/reference/bg/wxs_snv_18may18/chr;
# mkdir chr;
# ls | grep ^chr[1-9XYM] | prl "mv {} /groups/wyattgrp/reference/bg/wxs_snv_18may18/chr/{}";
# cat /groups/wyattgrp/reference/bg/wxs_snv_18may18/allchr_errorrates.txt | cut -f 1-3,16,17 | grep -v "[[:space:]]0\.0[[:space:]]0\.0" > /groups/wyattgrp/reference/bg/wxs_indel_18may18/allchr_errorrates.txt;

# # /groups/wyattgrp/reference/bg/73gp_snv_18may18/allchr_errorrates.txt
# # /groups/wyattgrp/reference/bg/73gp_indel_18may18/allchr_errorrates.txt
