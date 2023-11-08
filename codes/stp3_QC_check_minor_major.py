#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 11 2023

@author: s.hashemi

The function reads the bed file and checkes if a1 is the major allel and a0 is the minor or not(for controls only), if not then it swaps the allels.

INPUT: 
    ~/SNPs_chr'+str(i+1)
    ~/variants_chr'+str(i+1)
    
OUTPUT:    
    ~/SNPs_chr'+str(i+1)
    ~/variants_chr'+str(i+1)
"""

def QC_check_minor_major():

    import read_write as rw
    
    for i in range(24): #repeat it for all chromosomes
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        snp = rw.read_csv(file)
        dm = snp.shape
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-variants-chr/variants_chr'+str(i+1)
        var = rw.read_csv(file)
        
        for j in range(dm[0]):
     
            if snp.loc[j,'cntrl-a0a0'] > snp.loc[j,'cntrl-a1a1']: #swap allel a0 and a1  
                t = var.loc[j,'a0'] 
                var.loc[j,'a0'] = var.loc[j,'a1']
                var.loc[j,'a1'] = t            
               
                t = snp.loc[j,'PSC-a0a0'] 
                snp.loc[j,'PSC-a0a0'] = snp.loc[j,'PSC-a1a1']
                snp.loc[j,'PSC-a1a1'] = t
                
                t = snp.loc[j,'cntrl-a0a0'] 
                snp.loc[j,'cntrl-a0a0'] = snp.loc[j,'cntrl-a1a1']
                snp.loc[j,'cntrl-a1a1'] = t
                
                s = snp.loc[j,'snp']
                
                #find SNP postion
                l = len(s)
                semi = [j for j in range(len(s)) if s.startswith(':', j)]
                if len(semi) == 3:
                    p = s[0:semi[1]+1]
                    a1 = s[semi[1]+1:semi[2]]
                    a0 = s[semi[2]+1:l]
                    swapped_snp = p+a0+':'+a1
                    snp.loc[j,'snp'] = swapped_snp
                            
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-variants-chr/variants_chr'+str(i+1)
        rw.write_csv(var,file)
        
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        rw.write_csv(snp,file)