#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 20 2023

@author: s.hashemi

The function merges results of all quality control and statistical functions

INPUT: 
    ~/variants_chr'+str(i+1)
    ~/SNPs_chr'+str(i+1)
    ~/validSNPs_chr'+str(i+1)
    ~/MAF_chr'+str(i+1)
    ~/hwe_chr'+str(i+1)
    
OUTPUT:    
    ~/SNPs_info_chr'+str(i+1)
"""

def merge_files():
    
    import read_write as rw
    import pandas as pd
    
    for i in range(24): #repeat it for all chromosomes
       
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-variants-chr/variants_chr'+str(i+1)
        variants = rw.read_csv(file)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        snp = rw.read_csv(file)
        snp.drop('snp', inplace=True, axis=1)      
     
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-MAF-chr/MAF_chr'+str(i+1)
        maf = rw.read_csv(file)
        maf.drop('PSC(0= MAF<0.001)', 'cntrl(1=0.001<MAF<0.1)', 'all(2=0.1<MAF)',inplace=True, axis=1)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-validSNPs-chr/validSNPs_chr'+str(i+1)
        validSNPs = rw.read_csv(file)   
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/SNP-hwe-chr/hwe_chr'+str(i+1)
        hwe = rw.read_csv(file)   
        hwe.drop('HWE_PSC','HWE-all',inplace=True, axis=1)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC_final_step_validSNPs-chr/validSNPs_chr'+str(i+1)
        final_valid_SNPs = rw.read_csv(file) 
            
        data = pd.concat([variants, snp, validSNPs, maf, hwe, final_valid_SNPs], axis=1)    
        file = '/home/shashemi/Desktop/0-OneDrive/Output/SNPs-info-chr/SNPs_info_chr'+str(i+1)
        rw.write_csv(data,file)