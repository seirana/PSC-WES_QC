#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 11 023

@author: s.hashemi

The function reads the bed file and finds alleles and variants for a spesific postin on chromosoms, and fineds indels

INPUT: 
    ~/KIEL_Freeze_Three_OQFE.norm_PSCCON
    
OUTPUT:    
    ~/variants
"""
def QC_variants():

    import numpy as np
    import pandas as pd
    import read_write as rw
    
    file = '/home/shashemi/Desktop/0-OneDrive/Data/KIEL_Freeze_Three_OQFE.norm_PSCCON'
    bed, bim, fam = rw.read_bedbimfam(file)
    
    snps = bim.loc[:,['snp']]
    dm = bim.shape
    
    snp_lst = pd.DataFrame()
    snp_lst = pd.DataFrame(snp_lst, columns=['i','chr','position','a0','a1','multi allelic', 'indels']) 
    snp_lst[snp_lst.columns[0]] = [j for j in range(dm[0])]
    
    rep = np.zeros((dm[0],2), dtype=int)
    for i in range(dm[0]):
        s = snps.loc[i]
        s = s.to_string()
        
        #find SNP postion
        l = len(s)
        s = s[7:l]
        l = len(s)
        semi = [j for j in range(len(s)) if s.startswith(':', j)]
        snp_lst.loc[i,'chr'] = s[0:semi[0]]
        snp_lst.loc[i,'position'] = s[semi[0]+1:semi[1]] #only SNP position on a chromosome -- no chromosome number, no allels
        if len(semi) == 3:
            snp_lst.loc[i,'a0'] = s[semi[1]+1:semi[2]] 
            snp_lst.loc[i,'a1'] = s[semi[2]+1:l]
                
            if semi[2]-semi[1] != 2: #indels
                rep[i,1] = 1
                
            if l-semi[2] > 2: #indels
                rep[i,1] = 1 
        else:
            #indels
            rep[i,1] = 1
            snp_lst.loc[i,'a0'] = s[semi[1]+1:l] 
            snp_lst.loc[i,'a1'] = 'nan'
        
        if i > 0:
            if snp_lst.loc[i,'position'] == snp_lst.loc[i-1,'position']:
                if rep[i-1,0] == 0:
                    rep[i-1,0] = 1;
                rep[i,0] = rep[i-1,0]+1;
                
    snp_lst.loc[:,['multi allelic','indels']] = rep
    rw.write_csv(snp_lst,'/home/shashemi/Desktop/0-OneDrive/Output/QC-variants')