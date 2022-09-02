# -*- coding: utf-8 -*-
"""
@author: Elie
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import gc
pd.options.mode.chained_assignment = None

def main():
	dir = os.path.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument("-t", help="Path to tumor vstat file", dest="tumor", type=str, required=True)
	parser.add_argument("-n", help="Path to normal vstat file", dest="normal", type=str, required=True)
	#intersect with gc outside for parallel process with copynumber
	#parser.add_argument("-gc",help="Path to gc bed file" ,dest="be", type=str, required=True)
	parser.add_argument("-o",help="Output file path" ,dest="output", type=str, required=True)
	parser.add_argument("-mnd", help="minimum normal depth", dest="mindepth", type=int, default=10, required=False)
	args = parser.parse_args()
	args.tumor = os.path.expanduser(args.tumor)
	args.normal = os.path.expanduser(args.normal)
	args.output = os.path.expanduser(args.output)
	minimumdepth = int(args.mindepth)

	column_names = ["chrom", "pos", "ref", "reads_all", "reads_pp", "matches", "matches_pp", "mismatches", "mismatches_pp", "deletions", "deletions_pp", "insertions", "insertions_pp", "A", "A_pp", "C", "C_pp", "T", "T_pp", "G", "G_pp", "N", "N_pp", "sample"]
	column_types = {"chrom":str, "pos":np.int32, "ref":str, "reads_all":np.int32, "reads_pp":np.int32, "matches":np.int32, "matches_pp":np.int32, "mismatches":np.int32, "mismatches_pp":np.int32, "deletions":np.int32, "deletions_pp":np.int32, "insertions":np.int32, "insertions_pp":np.int32, "A":np.int32, "A_pp":np.int32, "C":np.int32, "C_pp":np.int32, "T":np.int32, "T_pp":np.int32, "G":np.int32, "G_pp":np.int32, "N":np.int32, "N_pp":np.int32, "sample":str}

	def load_data(tumor, normal, column_names, column_types):
		print("loading "+tumor)
		dft = pd.read_csv(tumor, sep='\t', low_memory=True, usecols=column_names, dtype=column_types, engine='c')
		print("loading "+normal)
		dfn = pd.read_csv(normal, sep='\t', low_memory=True, usecols=column_names, dtype=column_types, engine='c')
		return dft, dfn

	def merge_tumor_normal(dft, dfn, column_names):
		print("renaming headers")
		old_names = column_names
		normal_names = ["chrom", "pos", "ref", "reads_all_n", "reads_pp_n", "matches_n", "matches_pp_n", "mismatches_n", "mismatches_pp_n", "deletions_n", "deletions_pp_n", "insertions_n", "insertions_pp_n", "A_n", "A_pp_n", "C_n", "C_pp_n", "T_n", "T_pp_n", "G_n", "G_pp_n", "N_n", "N_pp_n", "sample_n"]
		tumor_names = ["chrom", "pos", "ref", "reads_all_t", "reads_pp_t", "matches_t", "matches_pp_t", "mismatches_t", "mismatches_pp_t", "deletions_t", "deletions_pp_t", "insertions_t", "insertions_pp_t", "A_t", "A_pp_t", "C_t", "C_pp_t", "T_t", "T_pp_t", "G_t", "G_pp_t", "N_t", "N_pp_t", "sample_t"]
		dfn.rename(columns=dict(zip(old_names, normal_names)), inplace=True)
		dft.rename(columns=dict(zip(old_names, tumor_names)), inplace=True)

		#merge frames by genomic position
		print("merging tumor and normal")
		df = pd.merge(dfn, dft, how='inner', on=['chrom','pos','ref'])
		#remove all statistics for paired reads 
		df = df.drop(columns=["A_pp_t", "C_pp_t", "T_pp_t", "G_pp_t", "N_pp_t", "A_pp_n", "C_pp_n", "T_pp_n", "G_pp_n", "N_pp_n", "reads_pp_n", "matches_pp_n", "mismatches_pp_n", "deletions_pp_n", "insertions_pp_n", "reads_pp_t", "matches_pp_t", "mismatches_pp_t", "deletions_pp_t", "insertions_pp_t"]).reset_index(drop=True)

		#this bit is for large memory savings when dfs are very large for exome/wgs, not really necessary for small panels though
		df[df.select_dtypes(['object']).columns] = df.select_dtypes(['object']).apply(lambda x: x.astype('category'))
		del dfn
		del dft
		gc.collect()
		return df

	def add_additional_columns(df):
		print("allele fractions")
		df['TAF_t'] = df['T_t'] / df['reads_all_t']
		df['CAF_t'] = df['C_t'] / df['reads_all_t']
		df['GAF_t'] = df['G_t'] / df['reads_all_t']
		df['AAF_t'] = df['A_t'] / df['reads_all_t']
		df['TAF_n'] = df['T_n'] / df['reads_all_n']
		df['CAF_n'] = df['C_n'] / df['reads_all_n']
		df['GAF_n'] = df['G_n'] / df['reads_all_n']
		df['AAF_n'] = df['A_n'] / df['reads_all_n']

		print("normal base")
		df['normal_max_value'] = df[['A_n', 'C_n', 'G_n', 'T_n']].max(axis=1)
		df['base_normal'] = 'A'
		df.loc[df['normal_max_value'] == df['T_n'], 'base_normal'] = 'T'
		df.loc[df['normal_max_value'] == df['C_n'], 'base_normal'] = 'C'
		df.loc[df['normal_max_value'] == df['G_n'], 'base_normal'] = 'G'
		print("tumor base")
		df['tumor_max_value'] = df[['A_t', 'C_t', 'G_t', 'T_t']].max(axis=1)
		df['base_tumor'] = 'A'
		df.loc[df['tumor_max_value'] == df['T_t'], 'base_tumor'] = 'T'
		df.loc[df['tumor_max_value'] == df['C_t'], 'base_tumor'] = 'C'
		df.loc[df['tumor_max_value'] == df['G_t'], 'base_tumor'] = 'G'
		return df

	dft, dfn = load_data(args.tumor, args.normal, column_names, column_types)
	df_tn = merge_tumor_normal(dft, dfn, column_names)
	df_tn = add_additional_columns(df_tn)
	#to reduce file size incase mindepth wasnt set earlier, 10 bases is an absolute minimimum
	df_tn = df_tn.query('(A_n > @minimumdepth) or (T_n > @minimumdepth) or (C_n > @minimumdepth) or (G_n > @minimumdepth)')
	print("writing tables to "+args.output)
	df_tn.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
	main()
