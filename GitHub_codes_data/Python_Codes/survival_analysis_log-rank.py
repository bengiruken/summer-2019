#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 22:00:03 2019

@author: bengi
"""

#step 1: read and process the clinical data

import Python_Codes.maf_clinical_data_processing as clinical
df_survival = clinical.survival_data()

#?????patient_group_1, patient_group_2
#hangi iki grubu karşılaştıracaksak 


def kaplan_meier(df_survival, patient_group_1, patient_group_2):
    