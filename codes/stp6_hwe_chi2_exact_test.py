#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 19 2023

@author: s.hashemi

The function calculates Hardy Weinberg Equilibrium (HWE) for PSC samples, controls and, all samples with snphwe function from snphwe 
chi2- square test with one degree freedom

INPUT: 
    ~/validSNPs_chr'+str(i+1)
    ~/SNPs_chr'+str(i+1)
    
OUTPUT:    
    ~/hwe_chr'+str(i+1)
"""

def QC_validSNPs():

    import numpy as np
    import read_write as rw
    import pandas as pd
    from snphwe import snphwe
    
    for i in range(24): #repeat it for all chromosomes
    
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-validSNPs-chr/validSNPs_chr'+str(i+1)
        validSNPs = rw.read_csv(file)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        snp = rw.read_csv(file)
        dm = snp.shape   
        
        hwe = np.zeros((dm[0],3), dtype=float)
        for j in range(dm[0]):  
            if validSNPs.loc[j,'valid SNPs (=1)'] == 1:                 
                
                hom1 = snp.loc[j,'PSC-a0a0'] # a0a0
                hets = snp.loc[j,'PSC-a0a1'] # a0a1
                hom2 = snp.loc[j,'PSC-a1a1'] # a1a1            
                hwe[j,0] = snphwe(hets, hom1, hom2)
                
                hom1 = snp.loc[j,'cntrl-a0a0'] # a0a0
                hets = snp.loc[j,'cntrl-a0a1'] # a0a1
                hom2 = snp.loc[j,'cntrl-a1a1'] # a1a1
                hwe[j,1] = snphwe(hets, hom1, hom2)
                
                hom1 = snp.loc[j,'PSC-a0a0']+snp.loc[j,'cntrl-a0a0'] # a0a0
                hets = snp.loc[j,'PSC-a0a1']+snp.loc[j,'cntrl-a0a1'] # a0a1
                hom2 = snp.loc[j,'PSC-a1a1']+snp.loc[j,'cntrl-a1a1'] # a1a1 
                hwe[j,2] = snphwe(hets, hom1, hom2)
                
            else:
                hwe[j] = float('nan')   
            
        HWE_test = pd.DataFrame(hwe, columns=['HWE_PSC', 'HWE-cntrl', 'HWE-all'])
        file = '/home/shashemi/Desktop/0-OneDrive/Output/SNP-hwe-chr/hwe_chr'+str(i+1)
        rw.write_csv(HWE_test,file)