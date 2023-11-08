#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 14 2023

@author: s.hashemi

The function checkes the validity of samples based on available SNPs

INPUT: 
    ~/KIEL_Freeze_Three_OQFE.norm_PSCCON
    ~/SNPs_per_chromosome
    ~/bed_chr'+str(i+1)
    
OUTPUT:    
    ~/SNPs_chr'+str(i+1)
"""


def QC_SNPs():
        
    import numpy as np
    import pandas as pd
    import read_write as rw
    
    file = '/home/shashemi/Desktop/0-OneDrive/Data/KIEL_Freeze_Three_OQFE.norm_PSCCON'
    bed, bim, fam = rw.read_bedbimfam(file)
    dm = bed.shape
    del bed, fam
    
    # smp = ['smp4509','smp5391'] #883 PCS, from ['smp4509','smp5391']; 4509 cntlr, from ['smp0','smp4508']
    
    file = '/home/shashemi/Desktop/0-OneDrive/Output/SNPs_per_chromosome'
    SNPcount = rw.read_xlsx(file)
    
    s = SNPcount.loc[0:23,['start']].to_numpy()
    e = SNPcount.loc[0:23,['end']].to_numpy()
    
    s = s.astype(int)
    e = e.astype(int)
    
    for i in range(24): #repeat it for all chromosomes
    
        snp_lst = pd.DataFrame()
        snp_lst = pd.DataFrame(snp_lst, columns=['snp','PSC-a0a0','PSC-a0a1','PSC-a1a1','PSC-nan','cntrl-a0a0','cntrl-a0a1','cntrl-a1a1','cntrl-nan'])
        snp_lst[snp_lst.columns[0]] = bim.loc[s[i,0]:e[i,0],'snp']  ####check loci = e[i]-s[i]+1
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/bed-chr/bed_chr'+str(i+1)
        bed = rw.read_parquet(file)
        dm = bed.shape
        
        stat = np.zeros((dm[0],8), dtype=int) 
        for j in range(dm[0]): # PSC
    
            indices = [k for k, x in enumerate(bed.loc[j,'smp4509':'smp5391']) if x == 0]
            stat[j,0] = len(indices) # a0a0
            
            indices = [k for k, x in enumerate(bed.loc[j,'smp4509':'smp5391']) if x == 1]
            stat[j,1] = len(indices) # a0a1
            
            indices = [k for k, x in enumerate(bed.loc[j,'smp4509':'smp5391']) if x == 2]
            stat[j,2] = len(indices) # a1a1
            
            stat[j,3] = 883 - stat[j,0]- stat[j,1]- stat[j,2] # nan   
            
        for j in range(dm[0]): # cntrl
    
            indices = [k for k, x in enumerate(bed.loc[j,'smp0':'smp4508']) if x == 0]
            stat[j,4] = len(indices) # a0a0
            
            indices = [k for k, x in enumerate(bed.loc[j,'smp0':'smp4508']) if x == 1]
            stat[j,5] = len(indices) # a0a1
            
            indices = [k for k, x in enumerate(bed.loc[j,'smp0':'smp4508']) if x == 2]
            stat[j,6] = len(indices) # a1a1
            
            stat[j,7] = 4509 - stat[j,4]- stat[j,5]- stat[j,6]# nan
                    
        snp_lst[snp_lst.columns[1:9]] = stat[:,0:8]
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        rw.write_csv(snp_lst,file)