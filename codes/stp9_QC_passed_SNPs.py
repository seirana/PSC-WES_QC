#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 7 2023

@author: s.hashemi
"""

#seperates the variants from position number

def QC_passed_SNPs():

    import numpy as np
    import read_write as rw  
    import pandas as pd
    
    l = 0
    for i in range(24): #repeat it for all chromosomes

        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-validSNPs-chr/validSNPs_chr'+str(i+1)
        validSNPs = rw.read_csv(file)
        indices = [j for j, x in enumerate(validSNPs[validSNPs.columns[0]]) if x == 1] #index of valid SNPs per chromosome  
        l = l + len(indices)

    thrs = 0.05 / l
    
    n_indices = []
    for i in range(24): #repeat it for all chromosomes
    
        file = '/home/shashemi/Desktop/0-OneDrive/Output/SNP-hwe-chr/hwe_chr'+str(i+1)
        hwe = rw.read_csv(file)
           
        dm = hwe.shape
        valid_snp = np.zeros((dm[0],1),dtype=int)
        
        indices = [j for j, x in enumerate(hwe['HWE-cntrl']) if x > thrs] #HWE p value is under threshold
        valid_snp[indices] = 1
        n_indices = n_indices + indices 
        
        valid_snp = pd.DataFrame(valid_snp, columns=['QC_final_step_validSNPs'])
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC_final_step_validSNPs-chr/validSNPs_chr'+str(i+1)
        rw.write_csv(valid_snp,file)

    file = '/home/shashemi/Desktop/0-OneDrive/Output/valid_SNPs.txt'
    np.savetxt(file, n_indices, fmt='%s', delimiter='\t')   