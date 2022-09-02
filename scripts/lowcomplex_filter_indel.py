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
	parser.add_argument("-indel", help="Path to somatic anno indel file", dest="indel", type=str, required=True)
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
	args.indel = os.path.expanduser(args.indel)
	args.ref = os.path.expanduser(args.ref)
	args.output = os.path.expanduser(args.output)
	#args.bm = args.bm.astype(float)
	#they provide path with basename
	#args.output = os.path.expanduser(args.output)
	base = os.path.basename(args.output)
	basename = os.path.splitext(base)[0]
	
	df = pd.read_csv(args.indel, sep='\t', low_memory=False)
	reference = pysam.Fastafile(args.ref)
	df.chrom = df.chrom.astype(str)

	#mapping errors
	df = df[~df['GeneDetail.refGene'].str.match('dist=')]
	df = df[~df['GeneDetail.knownGene'].str.match('dist=')]	
	df.replace('\s+', '_',regex=True,inplace=True)
	
	df = df.rename(columns={ 'start' : 'pos', 'end' : 'pos_end'})

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
	
	# df4 = df3[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']]
	# df4.to_csv(args.output+"_unfiltered.tsv", sep='\t', index=False, na_rep='NAN', header=True)


	flagdf = df3.query('(flag == "FLAG")')
	notflagdf = df3.query('(flag == "PASS")')
	notflagdf["filter"] = "PASS"

	flagdf['errordeletion'] = flagdf['mean_errordel'].apply(lambda x: args.am*x*args.bm).astype(float)
	flagdf['errorinsertion'] = flagdf['mean_errorins'].apply(lambda x: args.am*x*args.bm).astype(float)
	
	flagdf['DAFN10'] = args.nf * flagdf['DAF_n']
	flagdf['IAFN10'] = args.nf * flagdf['IAF_n']

	##########################
	# query somatic variants #
	##########################
	
	mrs = int(args.mrs)
	maf = float(args.maf)
	
	flagdel = flagdf.query('(del_n_zygosity != "het") & (DAF_t > DAFN10) & (deletions_t > @mrs) & (DAF_t > @maf) & (DAF_t > errordeletion)')
	flagins = flagdf.query('(ins_n_zygosity != "het") & (IAF_t > IAFN10) & (insertions_t > @mrs) & (IAF_t > @maf) & (IAF_t > errordeletion)')
	
	flagdel = flagdel[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']]
	flagins = flagins[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']]
	notflagdf = notflagdf[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']]

	columnnames = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'ExonicFunc.knownGene', 'AAChange.knownGene', 'avsnp147', 'ExAC_ALL', 'cosmic70', 'CLINSIG', 'CLNDBN', 'Kaviar_AF', 'CADD_raw', 'CADD_raw_rankscore', 'chrom', 'pos', 'pos_end', 'ref', 'tumor_allele', 'base_var', 'reads_all_t', 'reads_all_n', 'sample_t', 'sample_n', 'target', 'matches_t', 'mismatches_t', 'A_t', 'T_t', 'G_t', 'C_t', 'deletions_t', 'insertions_t', 'matches_n', 'mismatches_n', 'A_n', 'T_n', 'G_n', 'C_n', 'deletions_n', 'insertions_n', 'DAF_t', 'IAF_t', 'DAF_n', 'IAF_n', 'del_n_zygosity', 'ins_n_zygosity', 'mean_errordel', 'mean_errorins', 'LSEQ', 'RSEQ', 'MaskReg10', 'low_com', 'mapflag']

	
	somatic_indel = pd.DataFrame(columns=columnnames)

	somatic_indel = somatic_indel.append(notflagdf)
	somatic_indel = somatic_indel.append(flagdel)
	somatic_indel = somatic_indel.append(flagins)
	somatic_indel.to_csv(args.output, sep='\t', index=False, na_rep='NAN', header=True)
if __name__ == "__main__":
	main()
