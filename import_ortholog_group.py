#!/usr/bin/python

import numpy as np
import random
import pandas as pd
import sys
import time
from matplotlib import pyplot as plt
from collections import Counter

#f = open('causal_gene_ortholog_input.csv', encoding="utf-8", errors='replace')
#contents =f.read()
#contents #check if any error strings as ?  

''' 
#preprocessing
df=pd.read_csv('causal_gene_ortholog_input.csv',sep=',',index_col='Index' )
df=df.replace('\r','',regex=True)
with open('causal_gene_ortholog_input_clean.csv', "w+") as outfile:
    df.to_csv('causal_gene_ortholog_input_clean.csv',sep=',', header=True)
'''
#
df_n=pd.read_csv('causal_gene_ortholog_input_clean.csv', sep=',')

Arabidopsis_thaliana=[]
Glycine_max=[]
Brassica_rapa=[]
Solanum_lycopersicum=[]
Oryza_sativa_Japonica=[]
Oryza_sativa_Indica=[]
Setaria_italica=[]
Zea_mays=[]
Sorghum_bicolor=[]
Brachypodium_distachyon=[]
Hordeum_vulgare=[]

class species:
    def __init__(self,sp_name,raw):
        self.sp_name=sp_name
        self.raw= raw
        self.total_gene=[]
        self.total_gene_flat=[]
        self.count_per_entry=[]
        self.trait_cat=[]
        self.trait_cat_flat=[]
        self.trait_cat_uni=[]
        self.cate_count=[]
        self.main_cat= ['development','abiotic_stress_response','biotic_stress_response', 'color_texture_flavor_and_fragrance','other']
        self.ortho_ent=[]
        self.no_ortho=[]
        self.known_causal_count=[]
    def get_gene(self): 
        sp_name=self.sp_name
        df=self.raw
        self.known_causal_count=len(df[df['organism']==sp_name])
        for i in range(len(df[sp_name])):
            if str(df[sp_name][i])=='nan' : 
                self.count_per_entry.append(0)
                self.total_gene.append([])
                self.trait_cat.append(df['Trait_Category'][i])
            else: 
                self.total_gene.append(df[sp_name][i].split(','))
                self.count_per_entry.append(len(df[sp_name][i].split(',')))
                self.trait_cat.append(df['Trait_Category'][i])
    def uni_gene_list(self):         
        for i in range(len(self.total_gene)): 
            for j in self.total_gene[i]: 
                self.total_gene_flat.append(j)
                self.trait_cat_flat.append(self.trait_cat[i])
        total_gene_flat_a=[] 
        for i in self.total_gene_flat:
            if len(i.rsplit('.',1))>1:
                if len(i.rsplit('.',1)[1])<7: # for some maize IDs in special format 
                    total_gene_flat_a.append(i.rsplit('.',1)[0])
                else: 
                    total_gene_flat_a.append(i)
            else: 
                total_gene_flat_a.append(i)
            total_gene_flat_uni_before_f=list(set(total_gene_flat_a))
            self.total_gene_flat_uni=list(filter(None,total_gene_flat_uni_before_f)) # unique flat gene list *
        self.cat_dict=dict(zip(total_gene_flat_a, self.trait_cat_flat)) # for trait category search
        for i in range(len(self.total_gene_flat_uni)):
            self.trait_cat_uni.append(self.cat_dict.get(self.total_gene_flat_uni[i]))
    def cate_counter(self): 
        for i in self.main_cat: 
            self.cate_count.append(Counter(self.trait_cat_uni).get(i))
    def ortho_per_ent(self): 
        for i in self.total_gene: 
            self.ortho_ent.append (len(i))
            if len(i)<1:
                self.no_ortho.append(1)
            else: 
                self.no_ortho.append(0)
        non_zero=list(filter((0).__ne__, self.ortho_ent)) # remove entry with 0 ortholog
        average_ortho_per_ent=( (np.sum(non_zero)-self.known_causal_count)/ len(non_zero) ) #average number of ortholog per entry (excluding no-ortho and known causal genes)
        percent_no_otho= sum(self.no_ortho)/ (len(self.total_gene)-self.known_causal_count)  
        self.ortho_entry_fig= [average_ortho_per_ent,percent_no_otho] #  and percentage of entries with no-ortho in the target species (excluding known causal genes)



AT= species('Arabidopsis_thaliana',df_n)
AT.get_gene()
AT.uni_gene_list()
AT.cate_counter()
AT.ortho_per_ent()

GM=species('Glycine_max',df_n)
GM.get_gene()
GM.uni_gene_list()
GM.cate_counter()
GM.ortho_per_ent()

BR=species('Brassica_rapa',df_n)
BR.get_gene()
BR.uni_gene_list()
BR.cate_counter()
BR.ortho_per_ent()

SL=species('Solanum_lycopersicum',df_n)
SL.get_gene()
SL.uni_gene_list()
SL.cate_counter()
SL.ortho_per_ent()


OSj=species('Oryza_sativa_Japonica',df_n)
OSj.get_gene()
OSj.uni_gene_list()
OSj.cate_counter()
OSj.ortho_per_ent()

OSi=species('Oryza_sativa_Indica',df_n)
OSi.get_gene()
OSi.uni_gene_list()
OSi.cate_counter()
OSi.ortho_per_ent()


SI=species('Setaria_italica',df_n)
SI.get_gene()
SI.uni_gene_list()
SI.cate_counter()
SI.ortho_per_ent()

ZM=species('Zea_mays',df_n)
ZM.get_gene()
ZM.uni_gene_list()
ZM.cate_counter()
ZM.ortho_per_ent()
#ZM.total_gene_flat_uni

SB=species('Sorghum_bicolor',df_n)
SB.get_gene()
SB.uni_gene_list()
SB.cate_counter()
SB.ortho_per_ent()
#SB.total_gene_flat_uni


BD=species('Brachypodium_distachyon',df_n)
BD.get_gene()
BD.uni_gene_list()
BD.cate_counter()
BD.ortho_per_ent()


BV=species('Hordeum_vulgare',df_n)
BV.get_gene()
BV.uni_gene_list()
BV.cate_counter()
BV.ortho_per_ent()
BV.total_gene_flat_uni

spe_abbrv={'AT':'Arabidopsis_thaliana','GM':'Glycine_max', 'BR':'Brassica_rapa','SL':'Solanum_lycopersicum',
'OSj': 'Oryza_sativa_Japonica', 'OSi':'Oryza_sativa_Indica','SI':'Setaria_italica','ZM':'Zea_mays',
'SB':'Sorghum_bicolor', 'BD':'Brachypodium_distachyon','BV':'Hordeum_vulgare' }

### summary figures
def Stack_bar_plot(spe_count): 
    plt.figure(figsize=(5,8))
    enter_species=spe_count
    cate_coun_array=np.column_stack(enter_species) # each row is a category, will be made into bar
    N=len(enter_species)
    development_bar= plt.bar(np.arange(N), list(cate_coun_array[0]),color='#ffa500')
    abiotic_bar=plt.bar(np.arange(N), list(cate_coun_array[1]),bottom= cate_coun_array[0], color= '#00ced1')
    biotic_bar =plt.bar(np.arange(N), list(cate_coun_array[2]),bottom= (cate_coun_array[0]+cate_coun_array[1]), color= '#ff0000')
    color_bar=plt.bar(np.arange(N), list(cate_coun_array[3]),bottom= (cate_coun_array[0]+cate_coun_array[1]+cate_coun_array[2]),color= '#8470ff')
    other_bar=plt.bar(np.arange(N), list(cate_coun_array[4]),bottom= (cate_coun_array[0]+cate_coun_array[1]+cate_coun_array[2]+cate_coun_array[3]),color='#556b2f' )
    plt.ylabel('Putative causal gene')
    plt.xticks(np.arange(N), ('Arabidopsis_thaliana', 'Oryza_sativa_Japonica','Setaria_italica' ),rotation=90)
    plt.legend((development_bar[0], abiotic_bar[0],biotic_bar[0],color_bar[0],other_bar[0]), ('development','abiotic_stress_response','biotic_stress_response', 'color_texture_flavor_fragrance','other'))
    plt.tight_layout()
    #plt.show()
    plt.savefig('putative_causal_genes_spe_cate.tiff', dpi=900)

Stack_bar_plot([AT.cate_count,OS.cate_count,SI.cate_count] )

def ent_otho_spes(spe): 
    name=[]
    Avg_orth=[]
    NOP=[]
    for i in spe:
        name.append(spe_abbrv.get(i))
        Avg_orth.append(globals()[i].ortho_entry_fig[0])
        NOP.append(globals()[i].ortho_entry_fig[1])
    f1=plt.figure(figsize=(3,4.5)) # average ortho plot
    bars = name
    height=Avg_orth
    y_pos = np.arange(len(bars))
    plt.bar(y_pos, height)
    plt.ylabel('orthologs per entry')
    plt.ylim([1,round(max(height))+0.5])
    plt.xticks(y_pos, bars,rotation=90)
    plt.tight_layout()
    f1.savefig('average number of ortholog per entry.tiff', dpi=900)
    plt.clf()
    f2=plt.figure(figsize=(3,4.5)) # percent of no-orthog
    bars = name
    height=NOP 
    plt.bar(y_pos, height)
    plt.ylabel('percent of entry with no ortholog')
    plt.ylim([0,max(height)+0.1])
    plt.xticks(y_pos, bars,rotation=90)
    plt.tight_layout()
    f2.savefig('percent of entry with no ortholog.tiff', dpi=900)


ent_otho_spes(['AT','GM', 'BR','SL','OSj', 'OSi','SI','ZM','SB','BD','BV'])


### update feature list with putative causal genes.

AT_feature_list=pd.read_csv('Arabidopsis_features_v3.05_noMAF_noindel.csv',sep=',' )

for i in AT.total_gene_flat_uni: 
    if i in list(AT_feature_list['ID']): 
        if (AT_feature_list['class'][AT_feature_list['ID']==i]!=1).bool(): 
            row_index=AT_feature_list[AT_feature_list['ID']==i].index.tolist()
            AT_feature_list.iloc [row_index,29]=1 #  column index of 'class'

sum(AT_feature_list['class'])

with open('Arabidopsis_features_v3.05_ortholog.csv','w') as output: 
    AT_feature_list.to_csv('Arabidopsis_features_v3.05_ortholog.csv',sep=',', index=False, header=True)


OS_feature_list=pd.read_csv('rice_features_v1.3.11_3000poly_noindel.csv',sep=',' )

for i in OSj.total_gene_flat_uni: 
    if i in list(OS_feature_list['ID']): 
        if (OS_feature_list['class'][OS_feature_list['ID']==i]!=1).bool(): 
            row_index=OS_feature_list[OS_feature_list['ID']==i].index.tolist()
            OS_feature_list.iloc [row_index,28]=1 # class column index 28

sum(OS_feature_list['class'])

with open('rice_features_v1.3.11_ortholog.csv','w') as output: 
    OS_feature_list.to_csv('rice_features_v1.3.11_ortholog.csv',sep=',', index=False, header=True)


### export setaria feature  
with open('setaria_causal_ortho_labels.csv','w') as f: 
    pd.Series(SI.total_gene_flat_uni).to_csv('setaria_causal_ortho_labels.csv', sep=',',index=False, header=False)


