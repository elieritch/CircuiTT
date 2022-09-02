# -*- coding: utf-8 -*-
"""
@author: Elie
"""

from collections import OrderedDict
import os
import sys
import argparse
import numpy as np
import datetime
import pandas as pd
import math
import subprocess
from subprocess import call
import multiprocessing as mp
import gc
pd.options.mode.chained_assignment = None

def main():

	##########################
	# command line arguments #
	##########################

	dir = os.path.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument("-tnv", "--tnvstat", help="Path to tnvstat file created with merge script", dest="tnvstat", type=str, required=True)
	parser.add_argument("-bg", "--background", help="Path to background error file", dest="be", type=str, required=True)
	#intersect with gc outside for parallel process with copynumber
	#parser.add_argument("-gc",help="Path to gc bed file" ,dest="be", type=str, required=True)
	parser.add_argument("-o", "--output", help="Output file basename. Can include path to output eg /Path/to/output/basename_of_sample" ,dest="output", type=str, required=True)
	parser.add_argument("-bm", "--background_error_rate_filter", dest="bm", type=float, default=15.0, help="Background error rate multiplier. Remove Any variants that have an allele frequency below this float multiplied by the background error rate at the given position. Default=15.0")
	parser.add_argument("-nf", "--normal_filter", dest="nf", type=float, default=3.0, help="Tumor has to be this many times the allele frequency of the base found in the normal sample. Default=3.0")
	parser.add_argument("-maf", "--minimum_allele_frequency", dest="maf", type=float, default=0.01, help="Minimum allele frequency in the tumor sample to consider a base to be a somatic variant. Default=0.01")
	parser.add_argument("-mrs", "--minimum_read_support", dest="mrs", type=int, default=10, help="Minimum read support in the tumor sample to consider a base to be a somatic variant. Default=10")
	parser.add_argument("-mgf", "--minimum_germ_frequency", dest="mgf", type=float, default=0.30, help="Minimum allele frequency in the normal sample to consider a base to be a germline variant. Default=0.30")
	parser.add_argument("-mgrs", "--minimum_germ_read_support", dest="mgrs", type=int, default=20, help="Minimum read support in the normal sample to consider a base to be a germline variant. Default=20")
	
	parser.add_argument("-het_low", "--heterozygous_lower_bound", dest="hlb", type=float, default=0.30, help="Lower bound of germline allele frequency in the normal sample to consider a base to be heterozygous. Default=0.30")
	parser.add_argument("-het_high", "--heterozygous_upper_bound", dest="hub", type=float, default=0.90, help="Upper bound of germline allele frequency in the normal sample to consider a base to be heterozygous. Default=0.90")
	
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
		#tnvstat columns needed
		cols = ["chrom", "pos", "ref", "reads_all_n", "matches_n", "mismatches_n", "deletions_n", "insertions_n", "A_n", "C_n", "T_n", "G_n", "N_n", "sample_n", "reads_all_t", "matches_t", "mismatches_t", "deletions_t", "insertions_t", "A_t", "C_t", "T_t", "G_t", "N_t", "sample_t", "TAF_t", "CAF_t", "GAF_t", "AAF_t", "TAF_n", "CAF_n", "GAF_n", "AAF_n", "normal_max_value", "base_normal", "tumor_max_value", "base_tumor", "gc", "target"]

		typing = {"chrom": "category", "pos": np.int64, "ref": "category", "reads_all_n": np.int32, "matches_n": np.int32, "mismatches_n": np.int32, "deletions_n": np.int32, "insertions_n": np.int32, "A_n": np.int32, "C_n": np.int32, "T_n": np.int32, "G_n": np.int32, "N_n": np.int32, "sample_n": "category", "reads_all_t": np.int32, "matches_t": np.int32, "mismatches_t": np.int32, "deletions_t": np.int32, "insertions_t": np.int32, "A_t": np.int32, "C_t": np.int32, "T_t": np.int32, "G_t": np.int32, "N_t": np.int32, "sample_t": "category", "TAF_t": np.float16, "CAF_t": np.float16, "GAF_t": np.float16, "AAF_t": np.float16, "TAF_n": np.float16, "CAF_n": np.float16, "GAF_n": np.float16, "AAF_n": np.float16, "normal_max_value": np.int32, "base_normal": "category", "tumor_max_value": np.int32, "base_tumor": "category", "gc": np.float16, "target": "category"}
		
		df = pd.read_csv(path, sep='\t', usecols=cols, dtype=typing, low_memory=True, engine='c').query('(mismatches_t > 1)').reset_index(drop=True)
		gc.collect()
		return df
		
	def load_background_error(path):
		cols = ["chrom", "pos", "ref", "mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]

		typing = {"chrom": "category", "pos": np.int64, "ref": "category", "mean_errorAtoC":np.float32, "mean_errorAtoT":np.float32, "mean_errorAtoG":np.float32, "mean_errorTtoC":np.float32, "mean_errorTtoG":np.float32, "mean_errorTtoA":np.float32, "mean_errorGtoC":np.float32, "mean_errorGtoT":np.float32, "mean_errorGtoA":np.float32, "mean_errorCtoG":np.float32, "mean_errorCtoT":np.float32, "mean_errorCtoA":np.float32, "mean_errordel":np.float32, "mean_errorins":np.float32}
		
		df = pd.read_csv(path, sep='\t', usecols=cols, dtype=typing, low_memory=True, engine='c')
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
	all_data[["mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]] = all_data[["mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]].fillna(value=0.0)
	del tnv
	del bg
	gc.collect()
	
	#######################################################
	# Only call variants at strongly homozygous positions #
	#######################################################
	
	def assign_zygosity(df, hub, hlb):
		df['n_zygosity'] = 'hom'
		hub = float(hub) #upper bound default 0.90
		hlb = float(hlb) #lower bound default 0.30
		#if true then n_zygosity = het, if false than hom, 
		#if any base has allele frequency between 0.3 and 0.9 then its het, otherwise its hom.
		#np.where faster than .loc methods and becomes a numba function in wgs
		df['n_zygosity'] = np.where( ((df['TAF_n'] > hlb) & (df['TAF_n'] < hub)) | ((df['CAF_n'] > hlb) & (df['CAF_n'] < hub)) | ((df['GAF_n'] > hlb) & (df['GAF_n'] < hub)) | ((df['AAF_n'] > hlb) & (df['AAF_n'] < hub)), 'het', 'hom')
		#have to add a line if tumor base is not >0.3 and <0.7 and normal base is het then tomur zygosity is loh. 
		#for some reasons, the np.where converted some categories to strings, so put back to category type and collect garbage
		df[df.select_dtypes(['object']).columns] = df.select_dtypes(['object']).apply(lambda x: x.astype('category'))
		gc.collect()
		return df
	
	print('Assigning zygosity from sample '+basename+' '+str(datetime.datetime.now()))
	all_data = assign_zygosity(all_data, args.hub, args.hlb)
	print('Finished assigning zygosity from sample '+basename+' '+str(datetime.datetime.now()))
	
	###################################################
	# background error calculation and normal filters #
	###################################################

	#for exomes and smaller panels lambda is faster than basic math. 
	# But for whole genomes this was all converted to numba
	def error_calculation(df, background):
		df['errorAtoC'] = df['mean_errorAtoC'].apply(lambda x: x*background).astype(float)
		df['errorAtoG'] = df['mean_errorAtoG'].apply(lambda x: x*background).astype(float)
		df['errorAtoT'] = df['mean_errorAtoT'].apply(lambda x: x*background).astype(float)

		df['errorTtoC'] = df['mean_errorTtoC'].apply(lambda x: x*background).astype(float)
		df['errorTtoG'] = df['mean_errorTtoG'].apply(lambda x: x*background).astype(float)
		df['errorTtoA'] = df['mean_errorTtoA'].apply(lambda x: x*background).astype(float)

		df['errorCtoG'] = df['mean_errorCtoG'].apply(lambda x: x*background).astype(float)
		df['errorCtoA'] = df['mean_errorCtoA'].apply(lambda x: x*background).astype(float)
		df['errorCtoT'] = df['mean_errorCtoT'].apply(lambda x: x*background).astype(float)

		df['errorGtoC'] = df['mean_errorGtoC'].apply(lambda x: x*background).astype(float)
		df['errorGtoA'] = df['mean_errorGtoA'].apply(lambda x: x*background).astype(float)
		df['errorGtoT'] = df['mean_errorGtoT'].apply(lambda x: x*background).astype(float)
		return df

	def normal_filter_values(df, normaltimes):
		df['cAFN10'] = normaltimes * df['CAF_n']
		df['gAFN10'] = normaltimes * df['GAF_n']
		df['tAFN10'] = normaltimes * df['TAF_n']
		df['aAFN10'] = normaltimes * df['AAF_n']
		return df

	all_data = error_calculation(all_data, args.bm)
	all_data = normal_filter_values(all_data, args.nf)

	##########################
	# query somatic variants #
	##########################

	def query_somatic_variants(df, min_af, min_reads):
		min_af = float(min_af)
		min_reads = int(min_reads)
		somaticAtoC = df.query('(base_normal == "A") & (n_zygosity != "het") & (CAF_t > cAFN10) & (C_t > @min_reads) & (CAF_t > @min_af) & (CAF_t > errorAtoC)')
		somaticAtoC['base_var'] = 'C'
		somaticAtoG = df.query('(base_normal == "A") & (n_zygosity != "het") & (GAF_t > gAFN10) & (G_t > @min_reads) & (GAF_t > @min_af) & (GAF_t > errorAtoG)')
		somaticAtoG['base_var'] = 'G'
		somaticAtoT = df.query('(base_normal == "A") & (n_zygosity != "het") & (TAF_t > tAFN10) & (T_t > @min_reads) & (TAF_t > @min_af) & (TAF_t > errorAtoT)')
		somaticAtoT['base_var'] = 'T'
		somaticGtoC = df.query('(base_normal == "G") & (n_zygosity != "het") & (CAF_t > cAFN10) & (C_t > @min_reads) & (CAF_t > @min_af) & (CAF_t > errorGtoC)')
		somaticGtoC['base_var'] = 'C'
		somaticGtoA = df.query('(base_normal == "G") & (n_zygosity != "het") & (AAF_t > aAFN10) & (A_t > @min_reads) & (AAF_t > @min_af) & (AAF_t > errorGtoA)')
		somaticGtoA['base_var'] = 'A'
		somaticGtoT = df.query('(base_normal == "G") & (n_zygosity != "het") & (TAF_t > tAFN10) & (T_t > @min_reads) & (TAF_t > @min_af) & (TAF_t > errorGtoT)')
		somaticGtoT['base_var'] = 'T'
		somaticTtoC = df.query('(base_normal == "T") & (n_zygosity != "het") & (CAF_t > cAFN10) & (C_t > @min_reads) & (CAF_t > @min_af) & (CAF_t > errorTtoC)')
		somaticTtoC['base_var'] = 'C'
		somaticTtoG = df.query('(base_normal == "T") & (n_zygosity != "het") & (GAF_t > gAFN10) & (G_t > @min_reads) & (GAF_t > @min_af) & (GAF_t > errorTtoG)')
		somaticTtoG['base_var'] = 'G'
		somaticTtoA = df.query('(base_normal == "T") & (n_zygosity != "het") & (AAF_t > aAFN10) & (A_t > @min_reads) & (AAF_t > @min_af) & (AAF_t > errorTtoA)')
		somaticTtoA['base_var'] = 'A'
		somaticCtoT = df.query('(base_normal == "C") & (n_zygosity != "het") & (TAF_t > tAFN10) & (T_t > @min_reads) & (TAF_t > @min_af) & (TAF_t > errorCtoT)')
		somaticCtoT['base_var'] = 'T'
		somaticCtoG = df.query('(base_normal == "C") & (n_zygosity != "het") & (GAF_t > gAFN10) & (G_t > @min_reads) & (GAF_t > @min_af) & (GAF_t > errorCtoG)')
		somaticCtoG['base_var'] = 'G'
		somaticCtoA = df.query('(base_normal == "C") & (n_zygosity != "het") & (AAF_t > aAFN10) & (A_t > @min_reads) & (AAF_t > @min_af) & (AAF_t > errorCtoA)')
		somaticCtoA['base_var'] = 'A'
		somaticframes = [somaticAtoC, somaticAtoG, somaticAtoT, somaticGtoC, somaticGtoA, somaticGtoT, somaticTtoC, somaticTtoG, somaticTtoA, somaticCtoT, somaticCtoG, somaticCtoA]
		
		somaticSNV = pd.concat(somaticframes)
		somaticSNV = somaticSNV.drop(['cAFN10', 'gAFN10', 'tAFN10', 'aAFN10'], axis=1).reset_index(drop=True)
		return somaticSNV

	somaticSNV = query_somatic_variants(all_data, args.maf, args.mrs)	
	
	###########################
	# query germline variants #
	###########################

	def query_germline_variants(df, min_af, min_reads):
		min_af = float(min_af)
		min_reads = int(min_reads)
		germAtoC = df.query('(ref == "A") and (n_zygosity != "het") and (CAF_n > @min_af) and (C_n > @min_reads) and (CAF_n > errorAtoC)')
		germAtoC['base_var'] = 'C'
		germAtoG = df.query('(ref == "A") and (n_zygosity != "het") and (GAF_n > @min_af) and (G_n > @min_reads) and (GAF_n > errorAtoG)')
		germAtoG['base_var'] = 'G'
		germAtoT = df.query('(ref == "A") and (n_zygosity != "het") and (TAF_n > @min_af) and (T_n > @min_reads) and (TAF_n > errorAtoT)')
		germAtoT['base_var'] = 'T'																									
		germGtoC = df.query('(ref == "G") and (n_zygosity != "het") and (CAF_n > @min_af) and (C_n > @min_reads) and (CAF_n > errorGtoC)')
		germGtoC['base_var'] = 'C'
		germGtoA = df.query('(ref == "G") and (n_zygosity != "het") and (AAF_n > @min_af) and (A_n > @min_reads) and (AAF_n > errorGtoA)')
		germGtoA['base_var'] = 'A'
		germGtoT = df.query('(ref == "G") and (n_zygosity != "het") and (TAF_n > @min_af) and (T_n > @min_reads) and (TAF_n > errorGtoT)')
		germGtoT['base_var'] = 'T'
		germTtoC = df.query('(ref == "T") and (n_zygosity != "het") and (CAF_n > @min_af) and (C_n > @min_reads) and (CAF_n > errorTtoC)')
		germTtoC['base_var'] = 'C'
		germTtoG = df.query('(ref == "T") and (n_zygosity != "het") and (GAF_n > @min_af) and (G_n > @min_reads) and (GAF_n > errorTtoG)')
		germTtoG['base_var'] = 'G'
		germTtoA = df.query('(ref == "T") and (n_zygosity != "het") and (AAF_n > @min_af) and (A_n > @min_reads) and (AAF_n > errorTtoA)')
		germTtoA['base_var'] = 'A'
		germCtoT = df.query('(ref == "C") and (n_zygosity != "het") and (TAF_n > @min_af) and (T_n > @min_reads) and (TAF_n > errorCtoT)')
		germCtoT['base_var'] = 'T'
		germCtoG = df.query('(ref == "C") and (n_zygosity != "het") and (GAF_n > @min_af) and (G_n > @min_reads) and (GAF_n > errorTtoG)')
		germCtoG['base_var'] = 'G'
		germCtoA = df.query('(ref == "C") and (n_zygosity != "het") and (AAF_n > @min_af) and (A_n > @min_reads) and (AAF_n > errorTtoA)')
		germCtoA['base_var'] = 'A'
		
		germframes = [germAtoC, germAtoG, germAtoT, germGtoC, germGtoA, germGtoT, germTtoC, germTtoG, germTtoA, germCtoT, germCtoG, germCtoA]
		germSNV = pd.concat(germframes)
		return germSNV
		
	germSNV = query_germline_variants(all_data, args.mgf, args.mgrs)	
	###########################################
	# annotate these data frames with annovar #
	###########################################
	
	somaticSNV['pos_end'] = somaticSNV['pos']
	germSNV['pos_end'] = germSNV['pos']
	
	somaticSNV = somaticSNV[['chrom', 'pos', 'pos_end', 'ref', 'base_var', 'reads_all_n', 'reads_all_t', 'A_n', 'A_t', 'C_n', 'C_t', 'T_n', 'T_t', 'G_n', 'G_t', 'sample_n', 'sample_t', 'TAF_t', 'CAF_t', 'GAF_t', 'AAF_t', 'TAF_n', 'CAF_n', 'GAF_n', 'AAF_n', 'base_normal', 'base_tumor', 'matches_n', 'mismatches_n', 'deletions_n', 'insertions_n', 'matches_t', 'mismatches_t', 'deletions_t', 'insertions_t', 'n_zygosity', 'mean_errorAtoC', 'mean_errorAtoG', 'mean_errorAtoT', 'mean_errorCtoA', 'mean_errorCtoG', 'mean_errorCtoT', 'mean_errorGtoA', 'mean_errorGtoC', 'mean_errorGtoT', 'mean_errorTtoA', 'mean_errorTtoC', 'mean_errorTtoG']]

	germSNV = germSNV[['chrom', 'pos', 'pos_end', 'ref', 'base_var', 'reads_all_n', 'reads_all_t', 'A_n', 'A_t', 'C_n', 'C_t', 'T_n', 'T_t', 'G_n', 'G_t', 'sample_n', 'sample_t', 'TAF_t', 'CAF_t', 'GAF_t', 'AAF_t', 'TAF_n', 'CAF_n', 'GAF_n', 'AAF_n', 'base_normal', 'base_tumor', 'matches_n', 'mismatches_n', 'deletions_n',  'insertions_n', 'matches_t', 'mismatches_t', 'deletions_t', 'insertions_t', 'n_zygosity', 'mean_errorAtoC', 'mean_errorAtoG', 'mean_errorAtoT', 'mean_errorCtoA', 'mean_errorCtoG', 'mean_errorCtoT', 'mean_errorGtoA', 'mean_errorGtoC', 'mean_errorGtoT', 'mean_errorTtoA', 'mean_errorTtoC', 'mean_errorTtoG']]

	print("writing tables")
	
	if len(somaticSNV) > 0:
		somaticSNV.to_csv(args.output+"_somatic_snv.tsv", sep='\t', index=False, header=False)
		cmd_somatic = "/usr/bin/perl %s %s %s -buildver hg38 -out %s -remove -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel -operation g,g,f,f,f,f,f,f,f,f,f,f,f -nastring ." %(args.annovar_path, args.output+"_somatic_snv.tsv", args.humandb_path, args.output+"_annotated_somatic_snv.tsv")
		subprocess.call(cmd_somatic, shell=True)
		sanno = pd.read_csv(args.output+"_annotated_somatic_snv.tsv.hg38_multianno.txt", sep='\t', low_memory=False)
		ssnv = pd.read_csv(args.output+"_somatic_snv.tsv", sep='\t', low_memory=False, header=None)
		ssnv.columns = ['chrom', 'pos', 'pos_end', 'ref', 'base_var', 'reads_all_n', 'reads_all_t', 'A_n', 'A_t', 'C_n', 'C_t', 'T_n', 'T_t', 'G_n', 'G_t', 'sample_n', 'sample_t', 'TAF_t', 'CAF_t', 'GAF_t', 'AAF_t', 'TAF_n', 'CAF_n', 'GAF_n', 'AAF_n', 'base_normal', 'base_tumor', 'matches_n', 'mismatches_n', 'deletions_n', 'insertions_n', 'matches_t', 'mismatches_t', 'deletions_t', 'insertions_t', 'n_zygosity', 'mean_errorAtoC', 'mean_errorAtoG', 'mean_errorAtoT', 'mean_errorCtoA', 'mean_errorCtoG', 'mean_errorCtoT', 'mean_errorGtoA', 'mean_errorGtoC', 'mean_errorGtoT', 'mean_errorTtoA', 'mean_errorTtoC', 'mean_errorTtoG']
		annossnv = pd.merge(sanno, ssnv, how='left', right_on=['chrom', 'pos', 'ref', 'base_var'], left_on=['Chr', 'Start', 'Ref', 'Alt'])
		annossnv = annossnv.drop_duplicates(subset=['Chr', 'Start', 'End'])
		annossnv.to_csv(args.output+"_somatic_anno_snv.tsv", sep='\t', index=False)
		suffix1 = "_annotated_somatic_snv.tsv.hg38_multianno.txt"
		suffix2 = "_somatic_snv.tsv"
		inter1 = str(args.output)+suffix1
		inter2 = str(args.output)+suffix2
		cmd = f"rm {inter1}"
		subprocess.call(cmd, shell=True)
		cmd = f"rm {inter2}"
		subprocess.call(cmd, shell=True)
	else:
		pass
	
	if len(germSNV) > 0:
		germSNV.to_csv(args.output+"_germline_snv.tsv", sep='\t', index=False, header=False)
		cmd_germ = "/usr/bin/perl %s %s %s -buildver hg38 -out %s -remove -protocol refGene,knownGene,avsnp147,exac03,cosmic70,clinvar_20170130,kaviar_20150923,gnomad_exome,dbnsfp33a,dbscsnv11,hrcr1,mcap,revel -operation g,g,f,f,f,f,f,f,f,f,f,f,f -nastring ." %(args.annovar_path, args.output+"_germline_snv.tsv", args.humandb_path, args.output+"_annotated_germline_snv.tsv")
		subprocess.call(cmd_germ, shell=True)
		ganno = pd.read_csv(args.output+"_annotated_germline_snv.tsv.hg38_multianno.txt", sep='\t', low_memory=False)
		gsnv = pd.read_csv(args.output+"_germline_snv.tsv", sep='\t', low_memory=False, header=None)
		gsnv.columns = ['chrom', 'pos', 'pos_end', 'ref', 'base_var', 'reads_all_n', 'reads_all_t', 'A_n', 'A_t', 'C_n', 'C_t', 'T_n', 'T_t', 'G_n', 'G_t', 'sample_n', 'sample_t', 'TAF_t', 'CAF_t', 'GAF_t', 'AAF_t', 'TAF_n', 'CAF_n', 'GAF_n', 'AAF_n', 'base_normal', 'base_tumor', 'matches_n', 'mismatches_n', 'deletions_n', 'insertions_n', 'matches_t', 'mismatches_t', 'deletions_t', 'insertions_t', 'n_zygosity', 'mean_errorAtoC', 'mean_errorAtoG', 'mean_errorAtoT', 'mean_errorCtoA', 'mean_errorCtoG', 'mean_errorCtoT', 'mean_errorGtoA', 'mean_errorGtoC', 'mean_errorGtoT', 'mean_errorTtoA', 'mean_errorTtoC', 'mean_errorTtoG']
		annogsnv = pd.merge(ganno, gsnv, how='left', right_on=['chrom', 'pos', 'ref', 'base_var'], left_on=['Chr', 'Start', 'Ref', 'Alt'])
		annogsnv = annogsnv.drop_duplicates(subset=['Chr', 'Start', 'End'])
		annogsnv.to_csv(args.output+"_germline_anno_snv.tsv", sep='\t', index=False)
		suffix1 = "_annotated_germline_snv.tsv.hg38_multianno.txt"
		suffix2 = "_germline_snv.tsv"
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

