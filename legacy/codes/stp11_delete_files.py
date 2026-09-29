#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 8 2023

@author: s.hashemi
"""

def delete_files():

    import os, shutil
    
    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/bed-chr") #stp0
    
    os.remove("/home/shashemi/Desktop/0-OneDrive/Output/QC-variants.csv") #stp1_1
    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/QC-variants-chr") #stp1_2
    
    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/QC-SNPs-chr") #stp2
    
    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/QC-MAF-chr") #stp4
    
    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/SNP-hwe-chr") #stp5
    
    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/QC-validSNPs-chr") #stp6    

    shutil.rmtree("/home/shashemi/Desktop/0-OneDrive/Output/QC_final_step_validSNPs-chr") #stp8   