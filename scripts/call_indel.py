# -*- coding: utf-8 -*-
"""
@author: Elie
"""

from collections import OrderedDict
import os
import sys
import argparse
import datetime
import numpy as np
import pandas as pd
import math
import subprocess
from subprocess import call
import gc
import pysam
pd.options.mode.chained_assignment = None

def main():

	##########################
	# command line arguments #
	##########################
	
	dir = os.path.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument("-tnv", "--tnvstat", help="Path to tnvstat file created with merge script", dest="tnvstat", type=str, required=True)
	parser.add_argument("-bg", "--background", help="Path to background error file", dest="be", type=str, required=True)
	parser.add_argument("-bam", help="Path to bam file", dest="bam", type=str, required=True)
	parser.add_argument("-o", "--output", help="Output file basename. Can include path to output eg /Path/to/output/basename_of_sample" ,dest="output", type=str, required=True)
	parser.add_argument("-bm", "--background_error_rate_filter", dest="bm", type=float, default=15.0, help="Background error rate multiplier. Remove Any variants that have an allele frequency below this float multiplied by the background error rate at the given position. Default=15.0")
	parser.add_argument("-nf", "--normal_filter", dest="nf", type=float, default=3.0, help="Tumor has to be this many times the allele frequency of the base found in the normal sample. Default=3.0")
	parser.add_argument("-maf", "--minimum_allele_frequency", dest="maf", type=float, default=0.01, help="Minimum allele frequency in the tumor sample to consider a base to be a somatic variant. Default=0.01")
	parser.add_argument("-mrs", "--minimum_read_support", dest="mrs", type=int, default=10, help="Minimum read support in the tumor sample to consider a base to be a somatic variant. Default=10")
	parser.add_argument("-mgf", "--minimum_germ_frequency", dest="mgf", type=float, default=0.30, help="Minimum allele frequency in the normal sample to consider a base to be a germline variant. Default=0.30")
	parser.add_argument("-mgrs", "--minimum_germ_read_support", dest="mgrs", type=int, default=20, help="Minimum read support in the normal sample to consider a base to be a germline variant. Default=20")
	parser.add_argument("-annovar", "--annovar_path", help="Path to annovar" ,dest="annovar_path", type=str, default="$annovar_path")
	parser.add_argument("-humandb", "--annovar_humandbd_path", help="Path to humandb" ,dest="humandb_path", type=str, default="/groups/wyattgrp/software/annovar/annovar/humandb")
	args = parser.parse_args()
	args.tnvstat = os.path.expanduser(args.tnvstat)
	args.be = os.path.expanduser(args.be)
	args.output = os.path.expanduser(args.output)
	args.annovar_path = os.path.expanduser(args.annovar_path)
	args.humandb_path = os.path.expanduser(args.humandb_path)
	base = os.path.basename(args.output)
	basename = os.path.splitext(base)[0]




	#######################################################
	# Loading the files for tnvstata and background error #
	#######################################################

	def load_tnvstat(path):
		#tnvstat columns needed.
		cols = ["chrom", "pos", "ref", "reads_all_n", "matches_n", "mismatches_n", "deletions_n", "insertions_n", "A_n", "C_n", "T_n", "G_n", "N_n", "sample_n", "reads_all_t", "matches_t", "mismatches_t", "deletions_t", "insertions_t", "A_t", "C_t", "T_t", "G_t", "N_t", "sample_t", "TAF_t", "CAF_t", "GAF_t", "AAF_t", "TAF_n", "CAF_n", "GAF_n", "AAF_n", "normal_max_value", "base_normal", "tumor_max_value", "base_tumor", "gc", "target"]
		typing = {"chrom": "category", "pos": np.int64, "ref": "category", "reads_all_n": np.int32, "matches_n": np.int32, "mismatches_n": np.int32, "deletions_n": np.int32, "insertions_n": np.int32, "A_n": np.int32, "C_n": np.int32, "T_n": np.int32, "G_n": np.int32, "N_n": np.int32, "sample_n": "category", "reads_all_t": np.int32, "matches_t": np.int32, "mismatches_t": np.int32, "deletions_t": np.int32, "insertions_t": np.int32, "A_t": np.int32, "C_t": np.int32, "T_t": np.int32, "G_t": np.int32, "N_t": np.int32, "sample_t": "category", "TAF_t": np.float16, "CAF_t": np.float16, "GAF_t": np.float16, "AAF_t": np.float16, "TAF_n": np.float16, "CAF_n": np.float16, "GAF_n": np.float16, "AAF_n": np.float16, "normal_max_value": np.int32, "base_normal": "category", "tumor_max_value": np.int32, "base_tumor": "category", "gc": np.float16, "target": "category"}
		df = pd.read_csv(path, sep='\t', usecols=cols, dtype=typing, low_memory=True, engine='c').query('(deletions_t > 0) or (insertions_t > 0)').reset_index(drop=True)
		gc.collect()
		return df

	def load_background_error(path):
		cols = ["chrom", "pos", "ref", "mean_errordel", "mean_errorins"]
		typing = {"chrom": "category", "pos": np.int64, "ref": "category", "mean_errordel":np.float32, "mean_errorins":np.float32}
		df = pd.read_csv(path, sep='\t', usecols=cols, dtype=typing, low_memory=True, engine='c').query('(mean_errordel > 0) or (mean_errorins > 0)').reset_index(drop=True)
		gc.collect()
		return df

	print('Loading data from sample '+basename+' '+str(datetime.datetime.now()))
	tnv = load_tnvstat(args.tnvstat)
	print('Finished loading data from sample '+basename+' '+str(datetime.datetime.now()))

	print('Loading background error '+str(datetime.datetime.now()))
	bg = load_background_error(args.be)
	print('Finished loading background error '+str(datetime.datetime.now()))
	
	######################################
	# Merge tnvstat and background error #
	######################################

	all_data = pd.merge(tnv, bg, how='left', on=['chrom', 'pos', 'ref'])
	all_data[["mean_errordel", "mean_errorins"]] = all_data[["mean_errordel", "mean_errorins"]].fillna(value=0.0)
	del tnv
	del bg
	gc.collect()
	
	
	#######################################################
	# Only call variants at strongly homozygous positions #
	#######################################################

	#call tumor normal differences
	#print("computing normal indel zygosity")
	def assign_zygosity(df):
		df['del_n_zygosity'] = np.where( ((df['DAF_n'] > 0.20) & (df['DAF_n'] < 0.95)), 'het', 'hom' )
		df['ins_n_zygosity'] = np.where( ((df['IAF_n'] > 0.20) & (df['IAF_n'] < 0.95)), 'het', 'hom' )
		return df
		
	all_data['DAF_t'] = all_data['deletions_t'] / all_data['reads_all_t']
	all_data['IAF_t'] = all_data['insertions_t'] / all_data['reads_all_t']
	all_data['DAF_n'] = all_data['deletions_n'] / all_data['reads_all_n']
	all_data['IAF_n'] = all_data['insertions_n'] / all_data['reads_all_n']
	print('Assigning zygosity from sample '+basename+' '+str(datetime.datetime.now()))
	all_data = assign_zygosity(all_data)
	print('Finished assigning zygosity from sample '+basename+' '+str(datetime.datetime.now()))
	
	
	###################################################
	# background error calculation and normal filters #
	###################################################

	#for exomes and smaller panels lambda is faster than basic math. 
	# But for whole genomes this was all converted to numba
	def error_calculation(df, background):
		df['delerror'] = df['mean_errordel'].apply(lambda x: x*background).astype(float)
		df['inserror'] = df['mean_errorins'].apply(lambda x: x*background).astype(float)
		return df

	def normal_filter_values(df, normaltimes):
		df['DFN10'] = normaltimes * df['DAF_n']
		df['IFN10'] = normaltimes * df['IAF_n']
		return df

	all_data = error_calculation(all_data, args.bm)
	all_data = normal_filter_values(all_data, args.nf)
	
	##########################
	# doing deletion calling #
	##########################
	#for deletions just need to string together pilups row where index i and i+1 meet criteria for indel
		#then use reference in tnvstat to string together deleted sequence
	def query_somatic_dels(df, min_af, min_reads):
		min_af = float(min_af)
		min_reads = int(min_reads)
	
		deletions = df.query('(del_n_zygosity != "het") & (DAF_t > DFN10) & (deletions_t > @min_reads) & (DAF_t > @min_af) & (DAF_t > delerror)')
		deletions['base_var'] = 'DEL'

		if len(deletions) > 0:
			#code modified from:
			#https://stackoverflow.com/questions/22804067/how-to-get-start-and-end-of-ranges-in-pandas
			#this is the deletion section, 
			#first make start and end if deletion positions are consecutive
			dresult = deletions.sort_values(['chrom'], ascending=True).reset_index()
			dr1 = dresult.groupby('chrom').apply(lambda x: x.sort_values(['chrom', 'pos'])) #sort by position
			ds = dr1['pos']
			start = ds[ds.diff(1) != 1].reset_index(drop=True)
			end = ds[ds.diff(-1) != -1].reset_index(drop=True)
			dr1['start'] = dr1['pos']
			dfull = pd.DataFrame({'start': start, 'end': end}, columns=['start', 'end'])
			df1 = pd.merge(dr1, dfull, on=['start'])

			#now put deleted sequence in from reference
			d = []
			for i, r in dfull.iterrows():
				dleft = dfull.at[i, 'start']
				dright = dfull.at[i, 'end']
				ddf = dr1[dr1['pos'].between(dleft, dright, inclusive=True)]
				dsequence_list = ddf['ref'].tolist()
				dsequence = ''.join(dsequence_list)
				d.append({'start': dleft, 'end': dright, 'full_ref': dsequence})
			dnewf1 = pd.DataFrame(d)
			df2 = pd.merge(df1, dnewf1, on=['start', 'end'])
			df2['tumor_allele'] = '-'
			df2.drop(['ref'], axis=1, inplace=True)
			df2.rename(columns={'full_ref': 'ref'}, inplace=True)
			df2["length"] = (df2["end"] - df2["start"]) + 1
			df2.drop_duplicates(inplace=True)
			return df2
		else:
			return deletions
	
	somaticdel = query_somatic_dels(all_data, args.maf, args.mrs)
 
	###########################
	# doing insertion calling #
	###########################	

	def query_somatic_inss(df, min_af, min_reads):
		min_af = float(min_af)
		min_reads = int(min_reads)
		#insertions
		#start and end in insertion
		#reference = "-" in tnvstat
		insertions = df.query('(ins_n_zygosity != "het") and (IAF_t > IFN10) and (insertions_t > @min_reads) and (IAF_t > @min_af) & (IAF_t > inserror)')
		insertions['base_var'] = 'INS'
		insertions['ref'] = '-'
		#insertions['tumor_allele'] = 'AT' #going to have to fix this later. Very hackey

		if len(insertions) > 0:
			#https://stackoverflow.com/questions/22804067/how-to-get-start-and-end-of-ranges-in-pandas
			#this is the deletion section, 
			#first make start and end if deletion positions are consecutive
			inresult = insertions.sort_values(['chrom'], ascending=True).reset_index()
			inr1 = inresult.groupby('chrom').apply(lambda x: x.sort_values(['chrom', 'pos']))
			ins = inr1['pos']
			start = ins[ins.diff(1) != 1].reset_index(drop=True)
			end = ins[ins.diff(-1) != -1].reset_index(drop=True)
			inr1['start'] = inr1['pos']
			infull = pd.DataFrame({'start': start, 'end': end}, columns=['start', 'end'])
			inf1 = pd.merge(inr1, infull, on=['start'])

			#now put sequence deleted in reference
			inslist = []
			for i, r in infull.iterrows():
				inleft = infull.at[i, 'start']
				inright = infull.at[i, 'end']
				#indf = inr1[inr1['pos'].between(inleft, inright, inclusive=True)]
				#insequence_list = indf['ref'].tolist()
				#insequence = ''.join(insequence_list)
				inslist.append({'start': inleft, 'end': inright})
			innewf1 = pd.DataFrame(inslist)
			inf2 = pd.merge(inf1, innewf1, on=['start', 'end'])
			inf2['start'] = inf2['pos']
			#now given the range of the insertion, need to go into the bam file, and find the sequence of the insertion
			samfile = pysam.AlignmentFile(args.bam, "rb")
			insertiondata = pd.DataFrame([])
			for index, row in inf2.iterrows():
				insertion_options = []
				position_start=row["pos"]
				chromosome = row["chrom"]
				for pileupcolumn in samfile.pileup(chromosome, position_start - 1, position_start + 1): #needs min of 3bp range
					for read in pileupcolumn.pileups:
						if (read.indel > 0): # +1 is insertion, -1 deletion, 0 is match
							#position_end = position_start + read.indel - 1 #-1 for difference in 0 based or 1 based
							start_on_read = read.query_position_or_next
							end_on_read = read.query_position_or_next + read.indel
							length = read.indel
							sequence = read.alignment.query_sequence[start_on_read+1:end_on_read+1]
							insertion_options.append({'start': position_start, 'length': length, 'tumor_allele': sequence})
				ins_option = pd.DataFrame(insertion_options)
				ins_option['count'] = ins_option.groupby(['length','start','tumor_allele'])['length'].transform('count')
				ins_option = ins_option[['start','length','tumor_allele','count']]
				ins_option = ins_option.drop_duplicates().sort_values(by=['count'], ascending=False).reset_index()
				first_row = ins_option.iloc[[0], ins_option.columns.get_indexer(['start','length','tumor_allele','count'])]
				insertiondata = insertiondata.append(first_row)
				#start_insseq = start_insseq.append(first_row, ignore_index=True)
			inf2['start'] = inf2['pos']
			inf3 = pd.merge(inf2, insertiondata, on=['start'])
			inf3.drop_duplicates(inplace=True)
			return inf3
		else:
			return insertions
			
	somaticins = query_somatic_inss(all_data, args.maf, args.mrs)
		
	#########################################
	# concatenate deletion and/or insertion #
	#########################################

	if (len(somaticins) > 0) and (len(somaticdel) > 0):
		somaticframes = [somaticdel, somaticins]
		somatic_indel = pd.concat(somaticframes)
	elif (len(somaticins) > 0) and (len(somaticdel) == 0):
		somatic_indel = somaticins
	elif (len(somaticins) == 0) and (len(somaticdel) > 0):
		somatic_indel = somaticdel
	else:
		print("no indels present")
		
	###############################
	# subset columns and annotate #
	###############################
	
	if (len(somaticdel) > 0) or (len(somaticdel) > 0):
		indel_base_info = somatic_indel[['chrom', 'start', 'end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins']]

		print("writing tables")
		indel_base_info.to_csv(args.output+"_somatic_indel.tsv", sep='\t', index=False, header=False)

		cmd_somatic = "perl %s %s %s -buildver hg38 -out %s -remove -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel -operation g,g,f,f,f,f,f,f,f,f,f,f,f -nastring ." %(args.annovar_path, args.output+"_somatic_indel.tsv", args.humandb_path, args.output+"_annotated_somatic_indel.tsv")
		subprocess.call(cmd_somatic, shell=True)

		sanno = pd.read_csv(args.output+"_annotated_somatic_indel.tsv.hg38_multianno.txt", sep='\t', low_memory=False)
		sindel = pd.read_csv(args.output+"_somatic_indel.tsv", sep='\t', low_memory=False, header=None)
		sindel.columns = ['chrom', 'start', 'end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins']

		annoindel = pd.merge(sanno, sindel, how='left', right_on=['chrom', 'start', 'end', 'ref'], left_on=['Chr', 'Start', 'End', 'Ref'])
		annoindelu = annoindel.drop_duplicates()
		annoindelu.to_csv(args.output+"_somatic_anno_indel.tsv", sep='\t', index=False)
		suffix1 = "_annotated_somatic_indel.tsv.hg38_multianno.txt"
		suffix2 = "_somatic_indel.tsv"
		inter1 = str(args.output)+suffix1
		inter2 = str(args.output)+suffix2
		cmd = f"rm {inter1}"
		subprocess.call(cmd, shell=True)
		cmd = f"rm {inter2}"
		subprocess.call(cmd, shell=True)
	else:
		pass
if __name__ == "__main__":
	main()
