#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 14 2023

@author: s.hashemi

The function checkes the validity of samples based on available SNPs

INPUT: 
    ~/SNPs_per_chromosome
    ~/QC-variants
    
OUTPUT:    
    ~/variants_chr'+str(i+1)
"""

def seperate_QC_variants():

    import read_write as rw
    
    file = '/home/shashemi/Desktop/0-OneDrive/Output/SNPs_per_chromosome'
    SNPcount = rw.read_xlsx(file)
    
    s = SNPcount.loc[0:23,['start']].to_numpy()
    e = SNPcount.loc[0:23,['end']].to_numpy()
    
    s = s.astype(int)
    e = e.astype(int)
    
    file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-variants'
    variants = rw.read_csv(file)
    
    for i in range(24): #repeat it for all chromosomes
    
        variants_chr = variants[s[i,0]:e[i,0]+1]
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-variants-chr/variants_chr'+str(i+1)
        rw.write_csv(variants_chr,file)