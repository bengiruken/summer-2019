#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 22:32:48 2019

@author: bengi
"""

from Python_Codes.HGVNp_column_parsing import HGVNp_short_parse
import Python_Codes.read_maf_file_for_certain_gene as y
import Python_Codes.tumor_patient_id as z


y.create_folder("PIK3CA")
y.gene_mutants_from_pan_cancer_data("PIK3CA")
y.gene_specific_mutation_file("PIK3CA")
z.double_mutations("PIK3R1")
z.write_statistics_to_file("PIK3CA")
