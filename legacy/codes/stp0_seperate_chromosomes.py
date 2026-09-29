#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 14 2023

@author: s.hashemi

The function reperates .bed files seperate file based on chromosomes

INPUT: 
    ~/KIEL_Freeze_Three_OQFE.norm_PSCCON
    ~/SNPs_per_chromosome
    
OUTPUT:    
    ~/bed_chr'+str(i+1)
"""

def seperate_chromosomes():

    from dask import dataframe as dd
    import read_write as rw
    
    file = '/home/shashemi/Desktop/0-OneDrive/Data/KIEL_Freeze_Three_OQFE.norm_PSCCON'
    bed, bim, fam = rw.read_bedbimfam(file)
    
    file = '/home/shashemi/Desktop/0-OneDrive/Data/SNPs_per_chromosome'
    SNPcount = rw.read_xlsx(file)
    
    s = SNPcount.loc[0:23,['start']].to_numpy()
    e = SNPcount.loc[0:23,['end']].to_numpy()
    
    dm = bed.shape
    clm_id = []
    for i in range(dm[1]):
        clm_id.append('smp'+str(i))
        
    # store bed files per chromosome in parquet format
    path = '/home/shashemi/Desktop/0-OneDrive/Output/bed_chr'
    for i in range(24):   
        name = path+str(i+1)
        bed_chr = bed[s[i,0]:e[i,0]+1,:]
        bed_chr = dd.from_dask_array(bed_chr, columns=clm_id, index=None, meta=None) #change object Array -> DataFrame (dask.)
        dd.to_parquet(df=bed_chr,path=name) #equls to -> import bed_chr.to_parquet(file)  