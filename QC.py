#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Aug 30 2023

@author: s.hashemi

"""
from IPython import get_ipython
get_ipython().magic('reset -sf')

import os
os.system('clear')


import stp0_seperate_chromosomes as stp0
stp0.seperate_chromosomes()

import stp1_1_QC_variants as stp1_1
stp1_1.QC_variants()

import stp1_2_seperate_QC_variants as stp1_2
stp1_2.seperate_QC_variants()
    
import stp2_QC_SNPs as stp2
stp2.QC_SNPs()

import stp3_QC_check_minor_major as stp3
stp3.QC_check_minor_major()

import stp4.SNPs_percent as stp4
stp4.SNPs_percent()

import stp5_QC_MAF as stp5
stp5.QC_MAF()

import stp6_QC_validSNPs as stp6
stp6.QC_validSNPs()

import stp7_QC_hwe_chi2_exact_test as stp7
stp7.QC_hwe_chi2_exact_test()

import stp8_QC_samples as stp8
stp8.QC_samples()

import stp9_QC_passed_SNPs as stp9
stp9.QC_passed_SNPs

import stp10_merge_files as stp10
stp10.merge_files()

import stp11_delete_files as stp11
stp11.delete_files()