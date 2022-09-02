# -*- coding: utf-8 -*-
"""
@author: Elie
"""

import argparse
import pandas as pd
import os
pd.options.mode.chained_assignment = None

##########################
# command line arguments #
##########################

def main():
	dir = os.path.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument("-snv", help="Path to filtered snv file", dest="snv", type=str, required=True)
	parser.add_argument("-indel", help="Path to filtered indel file", dest="indel", type=str, required=True)
	parser.add_argument("-o",help="Path to output maf file" ,dest="output", type=str, required=True)
	args = parser.parse_args()
	args.snv = os.path.expanduser(args.snv)
	args.indel = os.path.expanduser(args.indel)
	args.output = os.path.expanduser(args.output)
	#args.bm = args.bm.astype(float)
	#they provide path with basename
	#args.output = os.path.expanduser(args.output)
	base = os.path.basename(args.output)

	snv = pd.read_csv(args.snv, sep='\t', low_memory=False)
	indel = pd.read_csv(args.indel, sep='\t', low_memory=False)

	snv = snv.rename(columns={ 'Gene.refGene' : 'Hugo_Symbol', 'Chr' : 'Chromosome', 'Start' : 'Start_Position','End' : 'End_Position', 'Ref' : 'Reference_Allele', 'Alt' : 'Tumor_Seq_Allele2', 'ExonicFunc.refGene' : 'Variant_Classification', 'sample_t' : 'Tumor_Sample_Barcode', 'reads_all_t' : 't_depth'})

	indel = indel.rename(columns={ 'Gene.refGene' : 'Hugo_Symbol', 'Chr' : 'Chromosome', 'Start' : 'Start_Position','End' : 'End_Position', 'Ref' : 'Reference_Allele', 'Alt' : 'Tumor_Seq_Allele2', 'ExonicFunc.refGene' : 'Variant_Classification', 'base_var' : 'Variant_Type', 'sample_t' : 'Tumor_Sample_Barcode', 'reads_all_t' : 't_depth'})

	snv["Variant_Type"] = "SNP"
	snv['VAF'] = snv["AAF_t"]
	snv.loc[snv['Tumor_Seq_Allele2'] == "T", 'VAF'] = snv["TAF_t"]
	snv.loc[snv['Tumor_Seq_Allele2'] == "G", 'VAF'] = snv["GAF_t"]
	snv.loc[snv['Tumor_Seq_Allele2'] == "C", 'VAF'] = snv["CAF_t"]

	indel['VAF'] = indel["DAF_t"]
	indel.loc[indel['Variant_Type'] == "INS", 'VAF'] = indel["IAF_t"]

	snv.loc[snv['Func.refGene'] == "splicing", 'Variant_Classification'] = "Splice_Site"
	indel.loc[indel['Func.refGene'] == "splicing", 'Variant_Classification'] = "Splice_Site"

	snv = snv.replace('synonymous_SNV', 'Synonymous_Mutation')
	snv = snv.replace('nonsynonymous_SNV', 'Missense_Mutation')
	snv = snv.replace('stopgain', 'Nonsense_Mutation')

	indel = indel.replace('frameshift_deletion', 'Frame_Shift_Del')
	indel = indel.replace('nonframeshift_deletion', 'In_Frame_Del')
	indel = indel.replace('frameshift_insertion', 'Frame_Shift_Ins')
	indel = indel.replace('nonframeshift_insertion', 'In_Frame_Ins')
	indel = indel.replace('stopgain', 'Nonsense_Mutation')

	snv.loc[snv['Variant_Classification'] == ".", 'Variant_Classification'] = "Non-Coding"
	indel.loc[indel['Variant_Classification'] == ".", 'Variant_Classification'] = "Non-Coding"

	snv["Protein_Change"] = ""
	indel["Protein_Change"] = ""

	if (len(snv) > 0):
		for i, r in snv.iterrows():
			barsnv=[]
			pcsnv = snv.at[i, "AAChange.refGene"].split(',')
			for item in pcsnv:
				#print(pc)
				#print(item)
				pc1snv = str(item)
				pc2snv = pc1snv.split(':')[-1]
				#print(pc2)
				#pc3 = pc1[-1]
				barsnv.append(pc2snv)
			uniquesnv = list(set(barsnv))
			finalsnv = ",".join(uniquesnv)
			snv.at[i, "Protein_Change"] = finalsnv
		snv = snv[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "VAF", "Protein_Change", "t_depth"]]
	else:
		pass
		
	if (len(indel) > 0):
		for i, r in indel.iterrows():
			barindel=[]
			pcindel = indel.at[i, "AAChange.refGene"].split(',')
			for item in pcindel:
				#print(pc)
				#print(item)
				pc1indel = str(item)
				pc2indel = pc1indel.split(':')[-1]
				#print(pc2)
				#pc3 = pc1[-1]
				barindel.append(pc2indel)
			uniqueindel = list(set(barindel))
			finalindel = ",".join(uniqueindel)
			indel.at[i, "Protein_Change"] = finalindel
		indel = indel[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "VAF", "Protein_Change", "t_depth"]]
	else:
		pass

	if (len(snv) > 0) and (len(indel) > 0):
		somaticframes = [snv, indel]
		maf = pd.concat(somaticframes)
		maf = maf.fillna("NaN")
	elif (len(snv) > 0) and (len(indel) < 1):
		maf = snv
		maf = maf.fillna("NaN")
	elif (len(snv) < 1) and (len(indel) > 0):
		maf = indel
		maf = maf.fillna("NaN")
	else:
		print("no mutations present")
		maf = snv[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "VAF", "Protein_Change", "t_depth"]]
		
	maf.to_csv(args.output, sep='\t', index=False)
if __name__ == "__main__":
	main()
