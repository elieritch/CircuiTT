#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: eritch
"""
import pysam
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
from collections import OrderedDict

##########################
# COMMAND LINE ARGUMENTS #
##########################

dir = os.path.dirname(__file__)
parser = argparse.ArgumentParser()
parser.add_argument("-snv", help="Path to somatic anno snv file", dest="snv", type=str, required=True)
parser.add_argument("-bam", help="Path to bamfile", dest="bam", type=str, required=True)
parser.add_argument("-md", "--mean_distance_from_end", dest="md", type=int, default=10, help="Include variants with mean distance from read ends greater than this value")
parser.add_argument("-o",help="Path to output file." ,dest="output", type=str, required=True)
args = parser.parse_args()
args.snv = os.path.expanduser(args.snv)
args.bam = os.path.expanduser(args.bam)
args.output = os.path.expanduser(args.output)
distance = args.md
base = os.path.basename(args.output)
basename = os.path.splitext(base)[0]

def main():

	###########################
	# SNV END-DISTANCE FILTER #
	###########################
	print("starting readend filter")
	df = pd.read_csv(args.snv, sep='\t', low_memory=False)
	df['chrom'] = df['chrom'].astype(str)
	samfile = pysam.AlignmentFile(args.bam, "rb")

	seq_columns = []
	for r, c in df.iterrows():
		chromosome = str(df.at[r, "chrom"])
		position = df.at[r, "pos"]
		position_i = position - 1 #position is zero indexed
		region_start = position_i - 2 #region +- 3 of position to accomadate 3bp indels
		region_end = position_i + 3
		alt = df.at[r, "Alt"]
		ref = df.at[r, "Ref"]

		distance_from_start = []
		distance_from_end = []
		for read in samfile.fetch(chromosome, region_start, region_end):
			read_positions = read.get_reference_positions(full_length=False) #false only get aligned position. If true then would inlude none values for softclip or other unaligned
			#because region = 7bp only use read if positions list contains position
			if position_i in read_positions:
				qs = read.query_sequence
				sequence_list = list(qs)
				#make dict with key=position and value=base
				pos_seq_dict = OrderedDict(zip(read_positions, sequence_list))
				#if position is variant value
				if pos_seq_dict[position_i] == alt:
					readstart = read.reference_start
					readend = read.reference_end
					var_from_start = int(position_i) - int(readstart)
					var_from_end = int(readend) - int(position_i)
					distance_from_start.append(var_from_start)
					distance_from_end.append(var_from_end)
		#take means of each list, then lists to dict to then dict to df.
		mean_distance_from_start = np.mean(distance_from_start)
		mean_distance_from_end = np.mean(distance_from_end)
		seq_columns.append({"chrom": chromosome, "pos": position, "mdfs": mean_distance_from_start, "mdfe": mean_distance_from_end})
	end_stuff = pd.DataFrame(seq_columns)
	#merge with df, remove rows with mean distances greater than, and drop the two extra columns.
	df_with_ends = pd.merge(end_stuff, df, on=['chrom','pos']).query('((mdfe > @distance) and (mdfs > @distance))').drop(['mdfe', 'mdfs'], axis=1).drop_duplicates(subset=['Chr', 'Start', 'End'])
	df_with_ends = df_with_ends[df.columns] #after new columns are removed, restore the order of columns from orig df.
	#df_with_ends = df_with_ends.drop_duplicates(subset=['Chr', 'Start', 'End'])
	df_with_ends.to_csv(args.output, sep='\t', index=False)
	print("finished readend filter")
if __name__ == "__main__":
	main()
