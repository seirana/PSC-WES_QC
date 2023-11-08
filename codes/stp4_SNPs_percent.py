#!/usr/bin/env python3
# -*- coding: utf-8 -*-
   
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sep 20 2023

@author: s.hashemi

The function  calculates thepercentage of PSC-a0a0, PSC-a0a1, cntrl-a0a0, and cntrl-a0a1 per each SNP location
                        it counts them for based on available samples, NAs are ignored

INPUT: 
    ~/SNPs_chr'+str(i+1)
    
OUTPUT:    
    ~/SNPs_percent'+str(i+1)
"""

def SNPs_percent():
    
    import read_write as rw
    import pandas as pd
    import numpy as np
    
    for i in range(24): #repeat it for all chromosomes
     
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        snp = rw.read_csv(file)   
        dm = snp.shape     
        
        wt = np.zeros((dm[0],6), dtype=float)
        
        sm = snp.loc[:,'PSC-a0a0']+snp.loc[:,'PSC-a0a1']+snp.loc[:,'PSC-a1a1']    
        wt[:,0] = snp.loc[:,'PSC-a0a0']/sm
        wt[:,1] = snp.loc[:,'PSC-a0a1']/sm
        wt[:,2] = snp.loc[:,'PSC-a1a1']/sm
        
        sm = snp.loc[:,'cntrl-a0a0']+snp.loc[:,'cntrl-a0a1']+snp.loc[:,'cntrl-a1a1']
        wt[:,3] = snp.loc[:,'cntrl-a0a0']/sm
        wt[:,4] = snp.loc[:,'cntrl-a0a1']/sm
        wt[:,5] = snp.loc[:,'cntrl-a1a1']/sm
       
        wt =np.around(wt*100,4)
        file = '/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr/SNPs_chr'+str(i+1)
        wt = pd.DataFrame(wt, columns=['PSC-a0a0%','PSC-a0a1%','PSC-a1a1%', 'cntrl-a0a0%','cntrl-a0a1%','cntrl-a1a1%'])
        data = pd.concat([snp, wt], axis=1)   
        rw.write_csv(data,file)