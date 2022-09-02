from collections import OrderedDict
import os
import sys
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import math
from scipy.stats import chi2
import gc
import datetime
pd.options.mode.chained_assignment = None

dir = os.path.dirname(__file__)
parser = argparse.ArgumentParser()
parser.add_argument("-in", help="input file" ,dest="input", type=str, required=True)
parser.add_argument("-out", help="output filename" ,dest="output", type=str, required=True)
# parser.add_argument("input_file", help="seqz file with gene name added", action='store')
# parser.add_argument("input_file", help="seqz file with gene name added", action='store')
args = parser.parse_args()
args.input = os.path.expanduser(args.input)
args.output = os.path.expanduser(args.output)
basename = os.path.splitext(args.input)[0]

# print('running script to create '+basename+' frequency file')
# print('Loading data from sample '+basename)
print("loading "+args.input)
currentDT = datetime.datetime.now()
print (str(currentDT))
#tp = pd.read_csv(args.input, sep='\t',iterator=True, chunksize=10000, low_memory=False)
#df = pd.concat(tp, ignore_index=True)
#df.chromosome = df.chromosome.astype(str)
df = pd.read_csv(args.input, sep='\t', low_memory=False, names=['chrom', 'pos', 'ref', 'reads_all', 'matches', 'mismatches', 'deletions', 'insertions', 'A', 'C', 'T', 'G', 'N', 'sample'])
df.chrom = df.chrom.astype(str)
print("finished loading "+args.input)
currentDT = datetime.datetime.now()
print (str(currentDT))
print('Proceeding with caculations...')
#del tp

#df.drop(columns=["reads_pp", "matches", "matches_pp", "mismatches", "mismatches_pp", "deletions", "deletions_pp", "insertions", "insertions_pp", "A_pp", "C_pp", "T_pp", "G_pp", "N_pp"], inplace=True)
df[df.select_dtypes(['object']).columns] = df.select_dtypes(['object']).apply(lambda x: x.astype('category'))
gc.collect()
print('finding alternates in for ref=A')
#find the number of each non reference base
df['alt_AtoC'] = df['C'].where(df['ref'] == "A", 0)
df['alt_AtoT'] = df['T'].where(df['ref'] == "A", 0)
df['alt_AtoG'] = df['G'].where(df['ref'] == "A", 0)
gc.collect()
print('finding alternates in for ref=T')
df['alt_TtoC'] = df['C'].where(df['ref'] == "T", 0)
df['alt_TtoG'] = df['G'].where(df['ref'] == "T", 0)
df['alt_TtoA'] = df['A'].where(df['ref'] == "T", 0)
gc.collect()
print('finding alternates in for ref=G')
df['alt_GtoC'] = df['C'].where(df['ref'] == "G", 0)
df['alt_GtoA'] = df['A'].where(df['ref'] == "G", 0)
df['alt_GtoT'] = df['T'].where(df['ref'] == "G", 0)
gc.collect()
print('finding alternates in for ref=C')
df['alt_CtoG'] = df['G'].where(df['ref'] == "C", 0)
df['alt_CtoA'] = df['A'].where(df['ref'] == "C", 0)
df['alt_CtoT'] = df['T'].where(df['ref'] == "C", 0)
gc.collect()
print('finished finding alternates')
currentDT = datetime.datetime.now()
print (str(currentDT))

#get sums of total reads
#df['Data3'].groupby(df['Date']).transform('sum')
df['reads_total_for_pos'] = df['reads_all'].groupby(df['pos']).transform('sum')
gc.collect()
print('finished sum of reads all')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['sum_alt_AtoC'] = df['alt_AtoC'].groupby(df['pos']).transform('sum')
df['sum_alt_AtoT'] = df['alt_AtoT'].groupby(df['pos']).transform('sum')
df['sum_alt_AtoG'] = df['alt_AtoG'].groupby(df['pos']).transform('sum')
gc.collect()
print('finished sum of alt=A')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['sum_alt_TtoC'] = df['alt_TtoC'].groupby(df['pos']).transform('sum')
df['sum_alt_TtoG'] = df['alt_TtoG'].groupby(df['pos']).transform('sum')
df['sum_alt_TtoA'] = df['alt_TtoA'].groupby(df['pos']).transform('sum')
gc.collect()
print('finished sum of alt=T')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['sum_alt_GtoC'] = df['alt_GtoC'].groupby(df['pos']).transform('sum')
df['sum_alt_GtoA'] = df['alt_GtoA'].groupby(df['pos']).transform('sum')
df['sum_alt_GtoT'] = df['alt_GtoT'].groupby(df['pos']).transform('sum')
gc.collect()
print('finished sum of alt=G')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['sum_alt_CtoG'] = df['alt_CtoG'].groupby(df['pos']).transform('sum')
df['sum_alt_CtoA'] = df['alt_CtoA'].groupby(df['pos']).transform('sum')
df['sum_alt_CtoT'] = df['alt_CtoT'].groupby(df['pos']).transform('sum')
gc.collect()
print('finished sum of alt=C')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['sum_alt_del'] = df['deletions'].groupby(df['pos']).transform('sum')
df['sum_alt_ins'] = df['insertions'].groupby(df['pos']).transform('sum')
gc.collect()
print('finished sum of alt=indel')
currentDT = datetime.datetime.now()
print (str(currentDT))

#df['mean_errordel'] = df['sum_alt_del'] / df['reads_total_for_pos']
#df['mean_errorins'] = df['sum_alt_ins'] / df['reads_total_for_pos']
print('computing means')
df['mean_errorAtoC'] = df['sum_alt_AtoC'] / df['reads_total_for_pos']
df['mean_errorAtoT'] = df['sum_alt_AtoT'] / df['reads_total_for_pos']
df['mean_errorAtoG'] = df['sum_alt_AtoG'] / df['reads_total_for_pos']
gc.collect()
print('finished mean of alt=C')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['mean_errorTtoC'] = df['sum_alt_TtoC'] / df['reads_total_for_pos']
df['mean_errorTtoG'] = df['sum_alt_TtoG'] / df['reads_total_for_pos']
df['mean_errorTtoA'] = df['sum_alt_TtoA'] / df['reads_total_for_pos']
gc.collect()
print('finished mean of alt=T')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['mean_errorGtoC'] = df['sum_alt_GtoC'] / df['reads_total_for_pos']
df['mean_errorGtoA'] = df['sum_alt_GtoA'] / df['reads_total_for_pos']
df['mean_errorGtoT'] = df['sum_alt_GtoT'] / df['reads_total_for_pos']
gc.collect()
print('finished mean of alt=G')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['mean_errorCtoG'] = df['sum_alt_CtoG'] / df['reads_total_for_pos']
df['mean_errorCtoA'] = df['sum_alt_CtoA'] / df['reads_total_for_pos']
df['mean_errorCtoT'] = df['sum_alt_CtoT'] / df['reads_total_for_pos']
gc.collect()
print('finished mean of alt=C')
currentDT = datetime.datetime.now()
print (str(currentDT))

df['mean_errordel'] = df['sum_alt_del'] / df['reads_total_for_pos']
df['mean_errorins'] = df['sum_alt_ins'] / df['reads_total_for_pos']
gc.collect()
print('finished mean of alt=indel')
currentDT = datetime.datetime.now()
print (str(currentDT))

print('subsetting')
df = df[["chrom", "pos", "ref", "mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]]
error_cols = ["mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]
df[error_cols] = df[error_cols].round(decimals = 5)
print('finished subsetting')
currentDT = datetime.datetime.now()
print (str(currentDT))
gc.collect()

df.drop_duplicates(inplace=True)
gc.collect()
print('finished dropping dupplicates')
currentDT = datetime.datetime.now()
print (str(currentDT))

print('writing file '+args.output+"_chromosome_mean_errors.tsv")
df.to_csv(args.output+"_chromosome_mean_errors.tsv", sep='\t', index=False, header=False)
