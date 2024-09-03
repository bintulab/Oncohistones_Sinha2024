import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('pdf') # do this because environment does not have GUI backend
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(sys.argv[1])
# locus = sys.argv[2]

antibodies = ['H3K9me3', 'H3K36me2', 'H3K36me3', 'H3K27me3', 'H3Kac', 'HA-tag', 'IgG']
# elemOrder = ['AAVS1', 'ACTB', 'KCNQ1']
condOrder = ['nodox', 'dox']

df['Cell Type'] = pd.Categorical(df['Cell Type'], categories=['PARENT', 'K36M'], ordered=True)
df.Antibody = pd.Categorical(df.Antibody, categories=antibodies, ordered=True)
# df.Element = pd.Categorical(df.Element, categories=elemOrder, ordered=True)
df.Condition = pd.Categorical(df.Condition, categories=condOrder, ordered=True)
df = df.sort_values(by=['Condition', 'Cell Type', 'Antibody', 'Element'])

#plot signal for all conditions 
pal = ['#999999', '#333333']
# g = sns.FacetGrid(data=df, col='Antibody', hue='Condition', aspect=3, height=5)
g = sns.catplot(data=df, col='Element', y='Signal per kb', hue='Condition',
				x='Antibody', row='Cell Type', kind='bar', height=4, aspect=1, palette=pal, sharey=False)
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
for ax in g.axes.flat:
	_ = ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.tight_layout()
plt.savefig('bedgraphAUC_plot_col-antibody_sigperkb.png')
plt.close()

colors = sns.color_palette(['#BBBBBB', '#E61B22'])
sns.set_theme(style='ticks')
for locus in list(set(list(df['Element']))):
	g = sns.catplot(data=df[(df['Element']==locus)], col='Antibody', col_wrap=3, hue='Cell Type', sharey=False,
					x='Condition', y='Signal per kb', kind='bar', palette=colors, aspect=0.9, height=2)
# g = sns.FacetGrid(data=df[(df['Element']==locus)], col='Antibody', col_wrap=4, hue='Cell Type')
# g.map(sns.barplot, 'Condition', 'Signal per kb')
	g.set_titles(col_template='{col_name}')
	# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cell Line')
	sns.move_legend(g, loc='center', bbox_to_anchor=(0.58, 0.2), title='Cell Line')
	g.fig.suptitle(locus)
	plt.tight_layout()
	plt.savefig('bedgraphAUC_plot_faceted_barplot_by_Ab_%s.png' % locus, dpi=300)
	plt.close()


#split dataframe by dox/no dox
df_nodox = df[df['Condition']=="nodox"]
print(df_nodox)
df_dox = df[df['Condition']=="dox"]
print(df_dox)

#calculate signal difference 
sig_diff = np.array(df_dox['Signal'].values.tolist())-np.array(df_nodox['Signal'].values.tolist())
print(sig_diff)
df_dox["Signal_Difference"]=sig_diff
#normalize signal difference by length of element 
df_dox['Length'] = df_dox['End'].astype(int) - df_dox['Start'].astype(int)
df_dox['Signal Difference per kb'] = 1000*df_dox['Signal_Difference']/df_dox['Length']

#plot signal difference between dox and no dox for each antibody 
# df_dox = df_dox[df_dox['Element']=='AAVS1']
# g = sns.catplot(data=df_dox, col='Element', y='Signal Difference per kb', x='Antibody',
# 				hue='Cell Type', kind='bar', height=4, aspect=1, palette=pal, sharey=False)
g = sns.catplot(data=df_dox, col='Antibody', y='Signal Difference per kb', x='Element',
				hue='Cell Type', kind='bar', height=4, aspect=1, palette=pal, sharey=False)
for ax in g.axes.flat:
	_ = ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.tight_layout()
plt.savefig('bedgraphAUC_plot_col-antibody_sig_diff_norm.png')




plt.close()