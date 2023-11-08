#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Aug 25 2023

@author: s.hashemi

This functions read and write different formats of files

Reading:
    .txt
    .parquete
    .xlsx
    .csv
    .bed
    
write:
    .csv    
"""

def read_txt(x):
    import pandas as pd
    with open(x) as f: 
        lines = f.readlines()
        
    l = len(lines)
    clmns = lines[0]  
    clmns = list(clmns.split(","))

    y = pd.DataFrame()
    y = pd.DataFrame(y, columns=clmns )
    for i in range(l):
        line_i = lines[i]  
        y.loc[i,:] = list(line_i.split(","))      
    return y
        
def read_parquet(x):
    # read partitioned parquet files
    import pandas as pd
    from pathlib import Path
    data_dir = Path(x)
    y = pd.concat(
        pd.read_parquet(parquet_file)
        for parquet_file in data_dir.glob('*.parquet')
    )
    return y

def read_xlsx(x):
    import pandas as pd     
    # read by default 1st sheet of an excel file
    y = pd.read_excel(x+'.xlsx')
    return y

def read_csv(x):
    # read by default 1st sheet of an csv file
    import pandas as pd
    y = pd.read_csv(x+'.csv')
    return y

def write_csv(x,y):
    from pathlib import Path  
    filepath = Path(y+'.csv')  
    filepath.parent.mkdir(parents=True, exist_ok=True)  
    x.to_csv(filepath,index=False)     

def read_bedbimfam(x):
	# Pandas_plink 
	# .bed file 
	# The genotype values can be either 0, 1, 2, or math.nan: 
	# 0 Homozygous having the first allele (given by coordinate a0) 
	# 1 Heterozygous 
	# 2 Homozygous having the second allele (given by coordinate a1) 
    # math.nan Missing genotype 
    from os.path import join
    from pandas_plink import read_plink
    from pandas_plink import get_data_folder
    (bim, fam, bed) = read_plink(join(get_data_folder(), x+'.bed'),
                                 verbose=False)
    
    return bed, bim, fam