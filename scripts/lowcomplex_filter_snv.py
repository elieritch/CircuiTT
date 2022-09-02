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
import pysam
pd.options.mode.chained_assignment = None
import re
import collections

def main():

	##########################
	# command line arguments #
	##########################

	dir = os.path.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument("-snv", help="Path to somatic anno snv file", dest="snv", type=str, required=True)
	parser.add_argument("-ref", help="Path to reference fasta that has been run through repeatmasker to lower the case of bases in repeats", dest="ref", type=str, required=True)
	#intersect with gc outside for parallel process with copynumber
	#parser.add_argument("-gc",help="Path to gc bed file" ,dest="be", type=str, required=True)
	parser.add_argument("-o",help="Path to output file" ,dest="output", type=str, required=True)
	parser.add_argument("-bm", "--background_error_rate_filter", dest="bm", type=float, default=15.0, help="Background error rate multiplier. Remove Any variants that have an allele frequency below this float multiplied by the background error rate at the given position. Default=15.0")
	parser.add_argument("-am", "--additional_multiplir_filter", dest="am", type=float, default=2.0, help="Additional background error rate multiplier. If the variant is flagged as suspicious due to low complexity regions, low mappability regions or next to repeats then the variant will need to have an allele frequency > (the background error)*(the above multiplier)*(this value). Default=2.0")
	parser.add_argument("-maxexac", "--maximum_exac_value", dest="maxexac", type=float, default=0.05, help="All variant must have no exac value ar be less than this value to be retained. Default=0.05")
	parser.add_argument("-nf", "--normal_filter", dest="nf", type=int, default=5, help="Variant flagged for filtering must have a higher allele frequency in the tumor than this value times the allele frequency in the normal. Default=5")
	parser.add_argument("-maf", "--minimum_allele_frequency", dest="maf", type=float, default=0.05, help="Minimum allele frequency in the tumor sample to consider a base to be a somatic variant for flagged variants. Default=0.05")
	parser.add_argument("-mrs", "--minimum_read_support", dest="mrs", type=int, default=20, help="Minimum read support in the tumor sample to consider a base to be a somatic variant for flagged variants. Default=20")
	args = parser.parse_args()
	args.snv = os.path.expanduser(args.snv)
	args.ref = os.path.expanduser(args.ref)
	args.output = os.path.expanduser(args.output)
	base = os.path.basename(args.output)
	basename = os.path.splitext(base)[0]
	
	df = pd.read_csv(args.snv, sep='\t', low_memory=False)
	reference = pysam.Fastafile(args.ref)
	df.chrom = df.chrom.astype(str)

	#mapping errors
	df = df[~df['GeneDetail.refGene'].str.match('dist=')]
	df = df[~df['GeneDetail.knownGene'].str.match('dist=')]	
	df.replace('\s+', '_',regex=True,inplace=True)

	df['ExAC_ALL'] = df['ExAC_ALL'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_AMR'] = df['ExAC_AMR'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_AFR'] = df['ExAC_AFR'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_EAS'] = df['ExAC_EAS'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_FIN'] = df['ExAC_FIN'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_NFE'] = df['ExAC_NFE'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_OTH'] = df['ExAC_OTH'].replace(to_replace=".", value=0).astype(float)
	df['ExAC_SAS'] = df['ExAC_SAS'].replace(to_replace=".", value=0).astype(float)
	df['Kaviar_AF'] = df['Kaviar_AF'].replace(to_replace=".", value=0).astype(float)
	maxexac = float(args.maxexac)
	df = df.query('((ExAC_ALL < @maxexac) or (ExAC_AFR < @maxexac) or (ExAC_AMR < @maxexac) or (ExAC_EAS < @maxexac) or (ExAC_FIN < @maxexac) or (ExAC_NFE < @maxexac) or (ExAC_OTH < @maxexac) or (ExAC_SAS < @maxexac) or (Kaviar_AF < @maxexac)) or (cosmic70 != ".")')

	# turn this into a flag
	df['mapflag'] = np.where( ((df['Gene.refGene'] != df['Gene.knownGene']) | (df['ExonicFunc.refGene'] != df['ExonicFunc.knownGene']) | (df['Func.refGene'] != df['Func.knownGene'])), 'MAPFLAG', 'NA' )
	
	seq_columns = []
	for r, c in df.iterrows():
		chromosome = str(df.at[r, "chrom"])
		position = df.at[r, "pos"]
		left_seq_start = position - 21
		left_seq_end = position - 1
		right_seq_start = position + 1
		right_seq_end = position + 21
		leftseq = reference.fetch(reference=chromosome, start=left_seq_start, end=left_seq_end)
		rightseq = reference.fetch(reference=chromosome, start=right_seq_start, end=right_seq_end)
		seq_columns.append({"chrom": chromosome, "pos": position, "LSEQ": leftseq, "RSEQ": rightseq})
	left_right_seqs = pd.DataFrame(seq_columns)
	df1 = pd.merge(df, left_right_seqs, on=["chrom", "pos"])

	masked_reg = []
	for r, c in df.iterrows():
		chromosome = str(df.at[r, "chrom"])
		position = df.at[r, "pos"]
		seq_start = position - 10
		seq_end = position + 10
		sequence = reference.fetch(reference=chromosome, start=seq_start, end=seq_end)
		MR10 = "MR10"
		notapp = "NA"
		if re.search(r"[acgt]", sequence):
			masked_reg.append({"chrom": chromosome, "pos": position, "MaskReg10": MR10})
		else:
			masked_reg.append({"chrom": chromosome, "pos": position, "MaskReg10": notapp})
	maskedregions = pd.DataFrame(masked_reg)
	df2 = pd.merge(df1, maskedregions, on=["chrom", "pos"])

	low_complex_list = []
	for r, c in df.iterrows():
		chromosome = str(df.at[r, "chrom"])
		position = df.at[r, "pos"]
		seq_start = position - 15
		seq_end = position + 15
		sequence = reference.fetch(reference=chromosome, start=seq_start, end=seq_end)
		sequence_list = list(sequence)
		total_len = len(sequence_list)
		counter = collections.Counter(sequence_list)
		maxbase = max(counter.values())
		fraction = maxbase / total_len
		low_comp_indicator = "lowcom"
		notapp = "NA"
		if fraction > 0.8:
			low_complex_list.append({"chrom": chromosome, "pos": position, "low_com": low_comp_indicator})
		else:
			low_complex_list.append({"chrom": chromosome, "pos": position, "low_com": notapp})
	low_complex_df = pd.DataFrame(low_complex_list)
	df3 = pd.merge(df2, low_complex_df, on=["chrom", "pos"])

	#cluster50bp
	# c50_list = []
	# for r, c in df.iterrows():
		# chromosome = str(df.at[r, "chrom"])
		# position = df.at[r, "pos"]
		# position_left = df.at[r, "Start"]
		# position_right = df.at[r, "End"]
		# left_window = position_left - 25
		# right_window = position_right + 25
		# notapp = "NA"
		# C50bp_indicator = "C50bp"
		# values_plus_minus_50 = df3.query('(chrom == @chromosome) and (pos > @left_window) and (pos < @right_window)')
		# if len(values_plus_minus_50) > 1:
			# c50_list.append({"chrom": chromosome, "pos": position, "C50bp": C50bp_indicator})
		# else:
			# c50_list.append({"chrom": chromosome, "pos": position, "C50bp": notapp})
	# c50_df = pd.DataFrame(c50_list)
	# df4 = pd.merge(df3, c50_df, on=["chrom", "pos"])

	df3['flag'] = np.where( ((df3['MaskReg10'] != "NA") | (df3['low_com'] != "NA") | (df3['Func.refGene'] != "exonic") | (df3['mapflag'] != "NA")), 'FLAG', 'PASS' )
	
	# df4 = df3[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'base_var', 'reads_all_n', 'reads_all_t', 'A_n', 'A_t', 'C_n', 'C_t', 'T_n', 'T_t', 'G_n', 'G_t', 'sample_n', 'sample_t', 'AAF_n', 'AAF_t', 'CAF_n', 'CAF_t', 'TAF_n', 'TAF_t', 'GAF_n', 'GAF_t', 'base_normal', 'base_tumor', 'matches_n', 'mismatches_n', 'matches_t', 'mismatches_t', 'deletions_n', 'insertions_n', 'deletions_t', 'insertions_t', 'n_zygosity', 'mean_errorAtoC', 'mean_errorAtoG', 'mean_errorAtoT', 'mean_errorCtoA', 'mean_errorCtoG', 'mean_errorCtoT', 'mean_errorGtoA', 'mean_errorGtoC', 'mean_errorGtoT', 'mean_errorTtoA', 'mean_errorTtoC', 'mean_errorTtoG', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']]
	# df4.to_csv(args.output+"_unfiltered.tsv", sep='\t', index=False, na_rep='NAN', header=True)


	flagdf = df3.query('(flag == "FLAG")')
	notflagdf = df3.query('(flag == "PASS")')
	notflagdf["filter"] = "PASS"

	flagdf['errorAtoC'] = flagdf['mean_errorAtoC'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorAtoG'] = flagdf['mean_errorAtoG'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorAtoT'] = flagdf['mean_errorAtoT'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorTtoC'] = flagdf['mean_errorTtoC'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorTtoG'] = flagdf['mean_errorTtoG'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorTtoA'] = flagdf['mean_errorTtoA'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorCtoG'] = flagdf['mean_errorCtoG'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorCtoA'] = flagdf['mean_errorCtoA'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorCtoT'] = flagdf['mean_errorCtoT'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorGtoC'] = flagdf['mean_errorGtoC'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorGtoA'] = flagdf['mean_errorGtoA'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorGtoT'] = flagdf['mean_errorGtoT'].apply(lambda x: args.am*x*args.bm).astype(float)
		
	flagdf['cAFN10'] = args.nf * flagdf['CAF_n']
	flagdf['gAFN10'] = args.nf * flagdf['GAF_n']
	flagdf['tAFN10'] = args.nf * flagdf['TAF_n']
	flagdf['aAFN10'] = args.nf * flagdf['AAF_n']

	##########################
	# query somatic variants #
	##########################
	
	mrs = int(args.mrs)
	maf = float(args.maf)
	somaticAtoC = flagdf.query('(base_normal == "A") & (n_zygosity != "het") & (CAF_t > cAFN10) & (C_t > @mrs) & (CAF_t > @maf) & (CAF_t > errorAtoC)')
	somaticAtoC['base_var'] = 'C'
	somaticAtoG = flagdf.query('(base_normal == "A") & (n_zygosity != "het") & (GAF_t > gAFN10) & (G_t > @mrs) & (GAF_t > @maf) & (GAF_t > errorAtoG)')
	somaticAtoG['base_var'] = 'G'
	somaticAtoT = flagdf.query('(base_normal == "A") & (n_zygosity != "het") & (TAF_t > tAFN10) & (T_t > @mrs) & (TAF_t > @maf) & (TAF_t > errorAtoT)')
	somaticAtoT['base_var'] = 'T'
	somaticGtoC = flagdf.query('(base_normal == "G") & (n_zygosity != "het") & (CAF_t > cAFN10) & (C_t > @mrs) & (CAF_t > @maf) & (CAF_t > errorGtoC)')
	somaticGtoC['base_var'] = 'C'
	somaticGtoA = flagdf.query('(base_normal == "G") & (n_zygosity != "het") & (AAF_t > aAFN10) & (A_t > @mrs) & (AAF_t > @maf) & (AAF_t > errorGtoA)')
	somaticGtoA['base_var'] = 'A'
	somaticGtoT = flagdf.query('(base_normal == "G") & (n_zygosity != "het") & (TAF_t > tAFN10) & (T_t > @mrs) & (TAF_t > @maf) & (TAF_t > errorGtoT)')
	somaticGtoT['base_var'] = 'T'
	somaticTtoC = flagdf.query('(base_normal == "T") & (n_zygosity != "het") & (CAF_t > cAFN10) & (C_t > @mrs) & (CAF_t > @maf) & (CAF_t > errorTtoC)')
	somaticTtoC['base_var'] = 'C'
	somaticTtoG = flagdf.query('(base_normal == "T") & (n_zygosity != "het") & (GAF_t > gAFN10) & (G_t > @mrs) & (GAF_t > @maf) & (GAF_t > errorTtoG)')
	somaticTtoG['base_var'] = 'G'
	somaticTtoA = flagdf.query('(base_normal == "T") & (n_zygosity != "het") & (AAF_t > aAFN10) & (A_t > @mrs) & (AAF_t > @maf) & (AAF_t > errorTtoA)')
	somaticTtoA['base_var'] = 'A'
	somaticCtoT = flagdf.query('(base_normal == "C") & (n_zygosity != "het") & (TAF_t > tAFN10) & (T_t > @mrs) & (TAF_t > @maf) & (TAF_t > errorCtoT)')
	somaticCtoT['base_var'] = 'T'
	somaticCtoG = flagdf.query('(base_normal == "C") & (n_zygosity != "het") & (GAF_t > gAFN10) & (G_t > @mrs) & (GAF_t > @maf) & (GAF_t > errorCtoG)')
	somaticCtoG['base_var'] = 'G'
	somaticCtoA = flagdf.query('(base_normal == "C") & (n_zygosity != "het") & (AAF_t > aAFN10) & (A_t > @mrs) & (AAF_t > @maf) & (AAF_t > errorCtoA)')
	somaticCtoA['base_var'] = 'A'
	somaticframes = [somaticAtoC, somaticAtoG, somaticAtoT, somaticGtoC, somaticGtoA, somaticGtoT, somaticTtoC, somaticTtoG, somaticTtoA, somaticCtoT, somaticCtoG, somaticCtoA]

	flaggedsomaticSNV = pd.concat(somaticframes)
	#somaticSNV.drop(['cAFN10', 'gAFN10', 'tAFN10', 'aAFN10'], axis=1)

	cols=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'base_var', 'reads_all_n', 'reads_all_t', 'A_n', 'A_t', 'C_n', 'C_t', 'T_n', 'T_t', 'G_n', 'G_t', 'sample_n', 'sample_t', 'AAF_n', 'AAF_t', 'CAF_n', 'CAF_t', 'TAF_n', 'TAF_t', 'GAF_n', 'GAF_t', 'base_normal', 'base_tumor', 'matches_n', 'mismatches_n', 'matches_t', 'mismatches_t', 'deletions_n', 'insertions_n', 'deletions_t', 'insertions_t', 'n_zygosity', 'mean_errorAtoC', 'mean_errorAtoG', 'mean_errorAtoT', 'mean_errorCtoA', 'mean_errorCtoG', 'mean_errorCtoT', 'mean_errorGtoA', 'mean_errorGtoC', 'mean_errorGtoT', 'mean_errorTtoA', 'mean_errorTtoC', 'mean_errorTtoG', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']
	flaggedsomaticSNV = flaggedsomaticSNV[cols]
	notflagdf = notflagdf[cols]
	dfs = [notflagdf, flaggedsomaticSNV]
	somaticsnv = pd.concat(dfs)

	somaticsnv.to_csv(args.output, sep='\t', index=False, na_rep='NAN', header=True)

if __name__ == "__main__":
	main()
