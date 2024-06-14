#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 22:48:06 2019

@author: bengi
"""


def tumor_type_stage_as_DataFrame():
    import pandas as pd
    xls = pd.ExcelFile("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Data/NIHMS978596-supplement-1.xlsx")
    df1 = pd.read_excel(xls, sheet_name="TCGA-CDR")
    subgroup=['bcr_patient_barcode', 'type']
    data= df1[df1[subgroup]]
    return data
tumor_type_stage_as_DataFrame()
##################

    
def survival_data():
    import pandas as pd
    xls = pd.ExcelFile("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Data/NIHMS978596-supplement-1.xlsx")
    df1 = pd.read_excel(xls, sheet_name="TCGA-CDR")
#    col_names=list(df1.columns) 
#    col_names
    subgroup=['bcr_patient_barcode', 
          'vital_status',
          #'OS', #vital status yerine
          'death_days_to',
          'last_contact_days_to', #days to last follow up
#          'OS',
          'OS.time'
          ]

    data= df1[subgroup]
    df = data[(data.vital_status != 'Dead') | ( data.vital_status != 'Alive')]
    #Replace dead > 1, alive >0
    df.vital_status.replace(['Dead', 'Alive'], [1, 0], inplace=True)
    df.death_days_to.fillna(0, inplace=True)
    df.last_contact_days_to.fillna(0, inplace=True)
    survival_data = df.dropna(subset = [ 'OS.time']) 
    survival_data = df.dropna(subset = [ 'vital_status']) 

    survival_data = survival_data[(survival_data['last_contact_days_to']>=0)]
    return survival_data

####################
