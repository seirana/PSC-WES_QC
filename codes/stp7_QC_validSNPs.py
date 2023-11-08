#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 18 2023

@author: s.hashemi

The function finds valid SNPs

INPUT: 
    ~/SNPs_per_chromosome
    ~/variants_chr'+str(i+1)
    ~/SNPs_chr'+str(i+1)
    
OUTPUT:    
    ~/validSNPs_chr'+str(i+1)
"""

def QC_hwe_chi2_exact_test():

    import numpy as np
    import pandas as pd
    import read_write as rw
    
    snp_lst = pd.DataFrame()
    snp_lst = pd.DataFrame(snp_lst, columns=['i','chr','position','will be deleted']) 
    
   
    omitted_SNPs = pd.DataFrame()
    omitted_SNPs = pd.DataFrame(omitted_SNPs, columns=['chr','multi allelic','nan','min 2 PSC with minor allele','all-omitted SNPs']) 
    
    for i in range(24): #repeat it for all chromosomes
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-variants-chr/variants_chr'+str(i+1)
        variants = rw.read_csv(file)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        SNPs = rw.read_csv(file)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-MAF-chr/MAF_chr'+str(i+1)
        maf = rw.read_csv(file)
            
        dm = SNPs.shape
        stat = np.ones((dm[0],1), dtype=float)
     
        indices = [j for j, x in enumerate(variants.loc[:,'multi allelic']) if x > 0] # multi allelic
        stat[indices] = 0  
        var1 = len(indices)
       
        trd = 0.05*883 #maximum 5% missing data for PSC samples is allowed
        indices = [j for j, x in enumerate(SNPs.loc[:,'PSC-nan']) if x > trd] # nan PSC
        stat[indices] = 0
        var2 = len(indices)
        
        trd = 0.05*4509 #maximum 5% missing data for control samples is allowed
        indices = [j for j, x in enumerate(SNPs.loc[:,'cntrl-nan']) if x > trd] # nan cntrls
        stat[indices] = 0
        var2 = var2 + len(indices)
           
        tmp_psc = SNPs.loc[:,'PSC-a0a0'] + SNPs.loc[:,'PSC-a0a1']
        indices_psc = [j for j, x in enumerate(tmp_psc) if x > 1] # at least two PSC samples with minor allele     
        tmp_cntrl = SNPs.loc[:,'cntrl-a0a0'] + SNPs.loc[:,'cntrl-a0a1']
        indices_cntrl = [j for j, x in enumerate(tmp_cntrl) if x > 1] # at least two cntrl samples with minor allele    
        indices = set(indices_psc).union(set(indices_cntrl))
        a = set(range(dm[0]))
        indices = list(a - indices)    
        stat[indices] = 0
        var3 = len(indices)
        
        indices = [j for j, x in enumerate(maf['cntrl(1=0.001<MAF<0.1)']) if x == 1]
        stat[indices] = 0
        var3 = var3 + len(indices)
                 
        omitted = int(dm[0] - sum(stat))
        
        # report
        omitted_SNPs.loc[i,'chr'] = i+1
        omitted_SNPs.loc[i,'multi allelic'] = str(var1)+'('+str(round(var1/dm[0]*100,2))+'%)'
        omitted_SNPs.loc[i,'nan'] = str(var2)+'('+str(round(var2/dm[0]*100,2))+'%)'
        omitted_SNPs.loc[i,'min 2 PSC with minor allele'] = str(var3)+'('+str(round(var3/dm[0]*100,2))+'%)'
        omitted_SNPs.loc[i,'all-omitted SNPs'] = str(omitted)+'('+str(round(omitted/dm[0]*100,2))+'%)'
    
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-validSNPs-chr/validSNPs_chr'+str(i+1)
        stat = pd.DataFrame(stat, columns=['valid SNPs (=1)'])
        rw.write_csv(stat,file)
        
    file = '/home/shashemi/Desktop/0-OneDrive/Output/validSNPs'
    omitted_SNPs = pd.DataFrame(omitted_SNPs)
    rw.write_csv(omitted_SNPs,file)