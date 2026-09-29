#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 11:14:00 2023

@author: s.hashemi

The function calculates MAF per SNP, for PSC cases, controls and all samples,
also checks if MAF is less and 0.001 or between 0.001 and 0.1 or is > 0.1

INPUT: 
    ~/SNPs_chr'+str(i+1)
    
OUTPUT:    
    ~/MAF_chr'+str(i+1)
"""

def QC_MAF():

    import read_write as rw
    import numpy as np
    import pandas as pd
    
    l_trs = 0.001
    u_trs = 0.1
    
    for i in range(24): #repeat it for all chromosomes
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        snp = rw.read_csv(file)
        dm = snp.shape
                
        maf = np.zeros((dm[0],6), dtype=float)  
        maf = pd.DataFrame(maf, columns=['PSC(0= MAF<0.001)', 'cntrl(1=0.001<MAF<0.1)', 'all(2=0.1<MAF)', 'MAF-PSC', 'MAF-cntrl', 'MAF-all'])
        for j in range(dm[0]):
            
            sum_psc = 2*(snp.loc[j,'PSC-a0a0'] + snp.loc[j,'PSC-a0a1'] + snp.loc[j,'PSC-a1a1'])
            if sum_psc == 0:
                print(i)
            sum_cntrl = 2*(snp.loc[j,'cntrl-a0a0'] + snp.loc[j,'cntrl-a0a1'] + snp.loc[j,'cntrl-a1a1'])
            if sum_psc == 0:
                print(i)
            sum_all = sum_psc + sum_cntrl
            
            maf_psc = (2*snp.loc[j,'PSC-a0a0'] + snp.loc[j,'PSC-a0a1'])/sum_psc        
            if  maf_psc > l_trs: 
                maf.loc[j,'PSC'] = maf_psc
                if maf_psc < u_trs:  
                    maf.loc[j,'PSC(0= MAF<0.001)'] = 1  
                else:
                    maf.loc[j,'PSC(0= MAF<0.001)'] = 2  
                    
            maf_cntrl = (2*snp.loc[j,'cntrl-a0a0'] + snp.loc[j,'cntrl-a0a1'])/sum_cntrl        
            if  maf_cntrl > l_trs:
                maf.loc[j,'cntrl'] = maf_cntrl
                if maf_cntrl < u_trs:  
                    maf.loc[j,'cntrl(1=0.001<MAF<0.1)'] = 1  
                else:
                    maf.loc[j,'cntrl(1=0.001<MAF<0.1)'] = 2  
                    
            maf_all = (maf_psc*sum_psc + maf_cntrl*sum_cntrl)/sum_all   
            if  maf_all > l_trs: 
                maf.loc[j,'all'] = maf_all
                if maf_all < u_trs:  
                    maf.loc[j,'all(2=0.1<MAF)'] = 1  
                else:
                    maf.loc[j,'all(2=0.1<MAF)'] = 2                  
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-MAF-chr/MAF_chr'+str(i+1)
        rw.write_csv(maf,file)