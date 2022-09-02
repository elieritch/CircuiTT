# -*- coding: utf-8 -*-
"""
@author: Elie
"""

import os
import sys
import argparse
import datetime
import gc

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
import seaborn as sns

import statsmodels.api as sm
import math
pd.options.mode.chained_assignment = None

def main():

	violin_dir = os.path.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument("-tnv", help="tnvstat file created from merge_vstat script", dest="input", type=str, required=True)
	parser.add_argument("-o", help="output basename including the path to output. ie /path/to/output/basename_of_sample. Several files will be made with this prefix", dest="output", type=str, required=True)
	parser.add_argument("-min_baf", help="minimum allele frequency to call allele as heterozygous", dest="minvaf", type=float, default=0.40, required=False)
	parser.add_argument("-max_baf", help="maximum allele frequency to call allele as heterozygous", dest="maxvaf", type=float, default=0.60, required=False)
	parser.add_argument("-min_rd", help="minimum read depth to call allele as heterozygous", dest="read_depth", type=int, default=100, required=False)
	args = parser.parse_args()
	args.input = os.path.expanduser(args.input)
	args.output = os.path.expanduser(args.output)
	base = os.path.basename(args.output)
	basename = os.path.splitext(base)[0]

	# =============================================================================
	## load data
	# =============================================================================

	def load_tnvstat(path):
		#tnvstat columns needed.
		cols = ["chrom", "pos", "ref", "reads_all_n", "matches_n", "mismatches_n", "deletions_n", "insertions_n", "A_n", "C_n", "T_n", "G_n", "N_n", "sample_n", "reads_all_t", "matches_t", "mismatches_t", "deletions_t", "insertions_t", "A_t", "C_t", "T_t", "G_t", "N_t", "sample_t", "TAF_t", "CAF_t", "GAF_t", "AAF_t", "TAF_n", "CAF_n", "GAF_n", "AAF_n", "normal_max_value", "base_normal", "tumor_max_value", "base_tumor", "gc", "target"]

		typing = {"chrom": "category", "pos": np.int64, "ref": "category", "reads_all_n": np.int32, "matches_n": np.int32, "mismatches_n": np.int32, "deletions_n": np.int32, "insertions_n": np.int32, "A_n": np.int32, "C_n": np.int32, "T_n": np.int32, "G_n": np.int32, "N_n": np.int32, "sample_n": "category", "reads_all_t": np.int32, "matches_t": np.int32, "mismatches_t": np.int32, "deletions_t": np.int32, "insertions_t": np.int32, "A_t": np.int32, "C_t": np.int32, "T_t": np.int32, "G_t": np.int32, "N_t": np.int32, "sample_t": "category", "TAF_t": np.float16, "CAF_t": np.float16, "GAF_t": np.float16, "AAF_t": np.float16, "TAF_n": np.float16, "CAF_n": np.float16, "GAF_n": np.float16, "AAF_n": np.float16, "normal_max_value": np.int32, "base_normal": "category", "tumor_max_value": np.int32, "base_tumor": "category", "gc": np.float16, "target": "category"}
		
		df = pd.read_csv(path, sep='\t', usecols=cols, dtype=typing, low_memory=True, engine='c').query('(target != ".")').reset_index(drop=True)
		gc.collect()
		return df

	print('Loading data from sample '+basename+' '+str(datetime.datetime.now()))
	df = load_tnvstat(args.input)
	print('Finished loading data from sample '+basename+' '+str(datetime.datetime.now()))

	# =============================================================================
	## median normalization and log2 of depth ratio
	# =============================================================================

	def median_normalization(df):
		tumor_median = df.loc[(df["chrom"]=="2") | (df["chrom"]=="4") | (df["chrom"]=="12") | (df["chrom"]=="14") | (df["chrom"]=="15") | (df["chrom"]=="19") | (df["chrom"]=="20"), ['reads_all_t']].apply(np.median)
		norm_median = df.loc[(df["chrom"]=="2") | (df["chrom"]=="4") | (df["chrom"]=="12") | (df["chrom"]=="14") | (df["chrom"]=="15") | (df["chrom"]=="19") | (df["chrom"]=="20"), ['reads_all_n']].apply(np.median)
		df['intnorm.depth.tumor'] = df['reads_all_t'] / int(tumor_median)
		df['intnorm.depth.normal'] = df['reads_all_n'] / int(norm_median)
		df['intnorm.depth.ratio'] = (df['intnorm.depth.tumor'] / df['intnorm.depth.normal'])
		df['log2.intnorm.depth.ratio'] = np.log2(df['intnorm.depth.ratio'])
		gc.collect()
		return df

	print('Starting median normalization from sample '+basename+' '+str(datetime.datetime.now()))	
	df = median_normalization(df)
	print('Finished median normalization from sample '+basename+' '+str(datetime.datetime.now()))

	# =============================================================================
	## lowess normalization 
	## produces plot for gc correction and adds the GC and median normalized 'Log2.depth.ratio'
	# =============================================================================

	def lowess_correction(df):
		lowess = sm.nonparametric.lowess
		x = df['gc'].values
		y = df['log2.intnorm.depth.ratio'].values
		ys = lowess(y, x, return_sorted=False)
		df['lowess'] = ys
		df['Log2.depth.ratio'] = df['log2.intnorm.depth.ratio'] - df['lowess']
		df.rename(columns={ 'target' : 'gene'}, inplace=True)
		return df

	print('Starting lowess normalization of '+basename+' '+str(datetime.datetime.now()))	
	df = lowess_correction(df)
	print('Finished lowess normalization of '+basename+' '+str(datetime.datetime.now()))

	dfsub = df.sample(frac=0.001)
	#graph lowess correction, before after
	print('Graphing lowess normalization of '+basename+' '+str(datetime.datetime.now()))
	sns.set_style("whitegrid")
	fig, ax = plt.subplots(figsize=(15,7), ncols=2, sharey='all')
	fig.suptitle('Lowess GC Correction '+basename)
	sns.regplot(x = "gc", y = "log2.intnorm.depth.ratio", data = dfsub, fit_reg=True, lowess=True, scatter_kws={"s": 5}, ax=ax[0])
	ax[0].set_xlim(0,100)
	ax[0].set_ylabel('log2 depth ratio')
	ax[0].set_title('Before subtract lowess')
	sns.regplot(x = "gc", y = "Log2.depth.ratio", data = dfsub, fit_reg=False, scatter_kws={"s": 5}, ax=ax[1])
	ax[1].set_xlim(0,100)
	ax[1].set_ylabel('log2 depth ratio')
	ax[1].set_title('After subtract lowess')
	print("Rendering GC corrected plot...")
	plt.savefig(args.output+'_gc_correction.png', dpi=300)
	plt.close()
	print('Finished graphing lowess normalization of '+basename+' '+str(datetime.datetime.now()))

	# =============================================================================
	## heterozygous and b alleles
	# =============================================================================

	def make_b_allele_table(df, minfrequency, maxfrequency, readdepth):
		minfreq = float(minfrequency)
		maxfreq = float(maxfrequency)
		rd = int(readdepth)
		df['n_zygosity'] = "hom"
		#where is faster than loc or truth table here for some reason
		df['n_zygosity'] = np.where( ((df['TAF_n'] > minfreq) & (df['TAF_n'] < maxfreq)) | ((df['CAF_n'] > minfreq) & (df['CAF_n'] < maxfreq)) | ((df['GAF_n'] > minfreq) & (df['GAF_n'] < maxfreq)) | ((df['AAF_n'] > minfreq) & (df['AAF_n'] < maxfreq)), 'het', 'hom')
		df_het = df.query('(n_zygosity == "het") and (reads_all_n > @rd) and (chrom != "X") and (chrom != "Y") and (chrom != "MT")')
		
		Major_freq_list = []
		Minor_freq_list = []
		for rows in df_het.itertuples():
			list_of_values = [float(rows.TAF_n),float(rows.CAF_n),float(rows.GAF_n),float(rows.AAF_n)]
			list_of_values = sorted(list_of_values)
			Major_freq = list_of_values[3]
			Minor_freq = list_of_values[2]
			Major_freq_list.append(Major_freq)
			Minor_freq_list.append(Minor_freq)
		df_het["Major_freq_n"] = Major_freq_list
		df_het["Minor_freq_n"] = Minor_freq_list
		
		df_het.loc[df_het['Minor_freq_n'] == df_het['TAF_n'], 'Minor_Allele'] = "T"
		df_het.loc[df_het['Minor_freq_n'] == df_het['AAF_n'], 'Minor_Allele'] = "A"
		df_het.loc[df_het['Minor_freq_n'] == df_het['CAF_n'], 'Minor_Allele'] = "C"
		df_het.loc[df_het['Minor_freq_n'] == df_het['GAF_n'], 'Minor_Allele'] = "G"

		#there will still be zeros in minor freq n. These are heterozygous germline indels. Dont graph but still output
		
		df_het['Major_Allele'] = df_het["base_normal"]
		df_het.loc[df_het['Minor_Allele'] == "T", 'B_allele_freq'] = df_het['TAF_t']
		df_het.loc[df_het['Minor_Allele'] == "C", 'B_allele_freq'] = df_het['CAF_t']
		df_het.loc[df_het['Minor_Allele'] == "G", 'B_allele_freq'] = df_het['GAF_t']
		df_het.loc[df_het['Minor_Allele'] == "A", 'B_allele_freq'] = df_het['AAF_t']
		df_het['B_allele_freq'] = df_het['B_allele_freq'].astype('float')
		#df_het = df_het.query('(Minor_Allele != ref)')
		df_het = df_het[['chrom', 'pos', 'ref', 'Log2.depth.ratio', 'n_zygosity', 'Major_freq_n', 'Minor_freq_n', 'Minor_Allele', 'Major_Allele', 'B_allele_freq']].reset_index(drop=True)
		return df_het

	print('Making b allele freq table '+basename+' '+str(datetime.datetime.now()))
	df_het = make_b_allele_table(df, args.minvaf, args.maxvaf, args.read_depth)
	print('Finished making b allele freq table '+basename+' '+str(datetime.datetime.now()))
	df_het.to_csv(args.output+'_b_allele_table.tsv',sep='\t', index=False)

	# =============================================================================
	## mean normalized log2 depth ratio per gene table
	# =============================================================================	

	def add_mean_gene(df):
		df['mean.gene.ndr'] = df['Log2.depth.ratio'].mean()
		return df
	mean_gene_table = df.groupby(['gene']).apply(add_mean_gene)

	print('making copynumber cn file...')
	sub_mean_gene_table = mean_gene_table[['gene', 'mean.gene.ndr']].drop_duplicates()
	sub_mean_gene_table['sample'] = basename
	sub_mean_gene_table.to_csv(args.output+".cn", sep='\t', index=False)

	# =============================================================================
	## make the violin plot
	# =============================================================================		

	aa = mean_gene_table['mean.gene.ndr'].values
	_, idx = np.unique(aa, return_index=True)
	aabb = aa[np.sort(idx)]

	color_array = np.asarray(aabb)
	color = []
	for i in np.nditer(color_array):
		if float(i) >= 0.7:
			color.append("bright red")
		elif float(i) >= 0.4 and float(i) < 0.7:
			color.append("peachy pink")
		elif float(i) > -0.4 and float(i) < 0.4:
			color.append("black")
		elif float(i) <= -0.4 and float(i) > -0.7:
			color.append("periwinkle")
		elif float(i) <= -0.7:
			color.append("strong blue")
		else:
			color.append("green")
			
	mean_gene_table['col5'] = mean_gene_table['Log2.depth.ratio']
	mean_gene_table['upper'] = mean_gene_table.groupby('gene').col5.transform(lambda p: p.mean() + 1.5*p.std())
	mean_gene_table['lower'] = mean_gene_table.groupby('gene').col5.transform(lambda p: p.mean() - 1.5*p.std())
	main_graphing_table = mean_gene_table[(mean_gene_table["Log2.depth.ratio"] <= mean_gene_table.upper) & (mean_gene_table["Log2.depth.ratio"] >= mean_gene_table.lower)]

	#graph the final violin.
	sns.set()
	sns.set_style("whitegrid")
	fig, ax = plt.subplots()
	fig.set_size_inches(18, 5.5) 
	title = str(basename+' CNV Violin Plot')
	ax.set_title(title, fontsize=18)
	#automatically adjust ylim (the ceiling (floor) of the max (min) mean log2-ratio across all genes is used to set ylim,; if this is above (below) 3 (-3), then ylim is set between -3 and 3 as default. A buffer of 1 is added to nearly guarantee the violins won't get cut off by the graph's edge in either case).
	table_to_get_minmax = mean_gene_table.copy(deep=True)
	table_to_get_minmax['mean'] = mean_gene_table.groupby('gene').col5.transform(lambda p: p.mean())
	ymax = math.ceil(table_to_get_minmax['mean'].max())
	ymin = math.floor(table_to_get_minmax['mean'].min())
	if ymax < 3:
		ymax = 2
	if ymin > (-3):
		ymin = (-2)
	ax.set_ylim(ymin-1, ymax+1)
	#apply shading by chrom
	shade_table_1 = main_graphing_table.drop_duplicates(subset=['chrom', 'gene'])
	shade_table_2 = shade_table_1.groupby(['chrom'], sort=False).size().reset_index(name='count').query('(count > 0)').reset_index(drop=True)
	shade_table_2 = shade_table_2.apply(pd.to_numeric, errors='coerce')
	shade_table_2 = shade_table_2.fillna('X')
	shade_table_2["index_of_chrom"] = shade_table_2.index

	#this loop does the shading and annotating
	shade_table_2["chrom"] = shade_table_2["chrom"].astype(str).replace('\.0', '', regex=True)
	total_genes = [0]
	for i in range(0,len(shade_table_2.index)):
		number_of_genes = int(shade_table_2['count'][i])
		min_number = sum(total_genes)
		total_genes.append(number_of_genes)
		max_number = sum(total_genes)
		if (i % 2 == 0): #even 
			ax.axvspan(min_number+0.5, max_number+0.5, color='0.75', alpha=0.4, zorder=0)
		xlocation = min_number + ((max_number-min_number)/2)
		chromosome_name = shade_table_2['chrom'][i]
		ax.annotate(chromosome_name, xy=(xlocation+0.5, ymax), xytext=(xlocation+0.5, ymax+0.2), size=12, ha="center", va="center")
		

	plt.xticks(rotation=90, ha='center')
	sns.despine()
	col_list_palette = sns.xkcd_palette(color)
	sns.violinplot(x=main_graphing_table["gene"], y=main_graphing_table["Log2.depth.ratio"], data=main_graphing_table, scale="width", width=1, cut=0, bw=.3, linewidth=0, inner=None, split=False, ax=ax, palette=col_list_palette)
	ax.set_ylabel(r'$\log_2$ Depth Ratio')
	ax.set_xlabel("Gene")
	ax.set_xlim(0.5, len(shade_table_1.index)+0.5)
	plt.tight_layout()
	print("Rendering "+basename+"_CN_violinplot.pdf")
	plt.savefig(args.output+"_CN_violinplot.pdf")
	print("Rendering complete")
	plt.close()
	
if __name__ == "__main__":
	main()
