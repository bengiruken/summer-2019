#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 22:00:03 2019

@author: bengi
"""

#step 1: read and process the clinical data

import Python_Codes.maf_clinical_data_processing as clinical
df_survival = clinical.survival_data()




  
def kaplan_meier(df_survival, patient_group_1, patient_group_2):
    """
    Perform Kaplan-Meier survival analysis for two patient groups.

    Parameters:
    - df_survival (DataFrame): DataFrame containing 'time' (time to event) and 'event' (1 if event occurred, 0 if censored) columns.
    - patient_group_1 (str): Label for patient group 1 (e.g., 'Group A').
    - patient_group_2 (str): Label for patient group 2 (e.g., 'Group B').

    Returns:
    - kmf (KaplanMeierFitter): KaplanMeierFitter object fitted with the survival data.
    """

    kmf = KaplanMeierFitter()

    # Split data into two groups
    group_1 = df_survival[df_survival['group'] == patient_group_1]
    group_2 = df_survival[df_survival['group'] == patient_group_2]

    # Fit Kaplan-Meier curves
    kmf.fit(group_1['time'], event_observed=group_1['event'], label=patient_group_1)
    ax = kmf.plot()

    kmf.fit(group_2['time'], event_observed=group_2['event'], label=patient_group_2)
    kmf.plot(ax=ax)

    plt.title('Kaplan-Meier Survival Curve')
    plt.xlabel('Time')
    plt.ylabel('Survival Probability')
    plt.show()

    return kmf


if __name__ == "__main__":
    import pandas as pd

 
    df = df_survival
    df.columns = ['time','event','group']

    # Perform Kaplan-Meier analysis for two patient groups
    kmf_result = kaplan_meier(df, 'Group A', 'Group B')

    # Print survival probabilities at specific times (optional)
    print("Survival probabilities at t=30:")
    

    
