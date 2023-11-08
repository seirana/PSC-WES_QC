#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 8 2023

@author: s.hashemi

The function checkes the validity of samples based on available SNPs

INPUT: 
    ~/KIEL_Freeze_Three_OQFE.norm_PSCCON
    ~/validSNPs_chr'+str(i+1)
    ~/SNPs_chr'+str(i+1)
    
OUTPUT:    
    ~/validSamples
"""

def QC_samples():
    
    import numpy as np
    import pandas as pd
    import read_write as rw
    
    file = '/home/shashemi/Desktop/0-OneDrive/Data/KIEL_Freeze_Three_OQFE.norm_PSCCON'
    bed, bim, fam = rw.read_bedbimfam(file)
    samples = fam.loc[:,['fid']]
    dm = samples.shape
    del bed, bim, fam
    
    nans_per_smpl = np.zeros((dm[0],1), dtype=int)
    snp_sum = 0
    for i in range(24): #repeat it for all chromosomes
    
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-validSNPs-chr/validSNPs_chr'+str(i+1)
        validSNPs = rw.read_csv(file)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/bed-chr/bed_chr'+str(i+1)
        bed = rw.read_parquet(file)
        
        indices = [j for j, x in enumerate(validSNPs[validSNPs.columns[0]]) if x == 1] #index of valid SNPs per chromosome
        bed = bed.loc[indices]
        snp_sum = snp_sum + len(indices)
        
        for k in range(dm[0]):
            indices = [j for j, x in enumerate(bed.loc[:,bed.columns[k]]) if np.isnan(x)] 
            nans_per_smpl[k] = nans_per_smpl[k] + len(indices)
            
    v_smps = pd.DataFrame()
    v_smps = pd.DataFrame(v_smps, columns=['samples','valid_SNPs'])
    v_smps[v_smps.columns[0]] = samples
    v_smps[v_smps.columns[1]] = np.around((1-nans_per_smpl/snp_sum)*100,2)
    
    file = '/home/shashemi/Desktop/0-OneDrive/Output/validSamples'
    rw.write_csv(v_smps,file)