#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 21:07:21 2019

@author: bengi
"""

''' "2_{}_point_mutations.txt".format(gene_name) dosyasını parse et'''
def tumor_ids(gene_name):
    filename = "/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/2_{}_point_mutations.txt".format(gene_name,gene_name)
    infile = open( "{}".format(filename) , 'r')
    tumors=[]
    #dosyadan tumor-mutated residue sözlüğü oluştur
    for line in infile:
        splitted=line.split("\t")
        #tumor id'ler 2. kolonda
        if splitted[2] not in tumors:
            tumors.append(splitted[2])
    infile.close()
    return  tumors 
x=tumor_ids('PIK3R1') 
#len(set(x))
#########33
def patient_ids(gene_name):
    filename = "/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/2_{}_point_mutations.txt".format(gene_name,gene_name)
    infile = open( "{}".format(filename) , 'r')
    patients=[]
    #dosyadan tumor-mutated residue sözlüğü oluştur
    for line in infile:
        splitted=line.split("\t")
        #tumor id'ler 2. kolonda
        if splitted[1] not in patients:
            patients.append(splitted[1])
    infile.close()
    return  patients   

#############
#y=patient_ids('PIK3R1') 
#len(set(y))
#
#print(y)
##########
'''find mutated residues for each patient'''

def patient_mutated_residues(gene_name):
    #tumor id'ler key, boş listeler value
    #dosyadan tumor-mutated residue sözlüğü oluştur
    filename = "/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/2_{}_point_mutations.txt".format(gene_name,gene_name)
    infile = open( "{}".format(filename) , 'r')
    patient_dict = {}
    for line in infile:
        splitted=line.split("\t")
        if splitted[1] not in patient_dict.keys():
            patient_dict[splitted[1]] = set()
            patient_dict[splitted[1]].add(splitted[-4])
        elif splitted[1] in patient_dict.keys():
            patient_dict[splitted[1]].add(splitted[-4])
    return patient_dict
#f=patient_mutated_residues('PIK3R1')
#patient_one_and_greater_one_mutated
def patient_one_and_greater_one_mutated(gene_name):
    one_mutated_residue = set() 
    more_than_one_mutated_residue = set()
    patient_mutated_residue = patient_mutated_residues(gene_name)
    for ids in patient_mutated_residue.keys():
        if len(patient_mutated_residue[ids]) == 1:
            one_mutated_residue.add(ids)
        else:
             more_than_one_mutated_residue.add(ids)
    return list(one_mutated_residue),  list(more_than_one_mutated_residue)

'''2 den fazla mutasyonu olan hastalarda mutasyonların 2'li kombinasyonlarını oluşturma'''
def double_mutations(gene_name):
    from itertools import combinations
    more_than_one_mutation = patient_one_and_greater_one_mutated(gene_name)[1]
    Dict = patient_mutated_residues(gene_name)
    double_mutation_set=set()
    for patient in more_than_one_mutation:
        each = [int(i) for i in  Dict[patient]]
        each=sorted(each)
        each = [str(i) for i in  each]
        double_mutation_set = double_mutation_set.union(set(combinations( each , 2 )))
    return double_mutation_set

#a=double_mutations('PIK3R1')     
#list(a)
#################
#hangi mutasyon çifti hangi hastalarda var
def double_mutation_patient_dict(gene_name):
    double_mutation_set = double_mutations(gene_name)
    patient_mutation_dict = patient_mutated_residues(gene_name)
    double_mut_pat_dictionary = {}
    for pair in double_mutation_set:
        for patient in  patient_mutation_dict.keys():
            if (pair[0] in patient_mutation_dict[patient] ) and (pair[1] in patient_mutation_dict[patient]):
                if pair not in double_mut_pat_dictionary.keys():
                    double_mut_pat_dictionary[pair]=set()
                    double_mut_pat_dictionary[pair].add(patient)
                elif pair in double_mut_pat_dictionary.keys():
                    double_mut_pat_dictionary[pair].add(patient)
    return double_mut_pat_dictionary


#write all sttaistics to file
def write_statistics_to_file(gene_name):
    import scipy
    import numpy as np
    from scipy import stats
    outfile1=open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/result_{}_equal_one.txt".format(gene_name,gene_name),"a")
    outfile2=open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/result_{}_greater_one.txt".format(gene_name,gene_name),"a")
    outfile3=open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/result_{}_not_having_both.txt".format(gene_name,gene_name),"a")
    pairs_dict=double_mutation_patient_dict(gene_name) 
    for pair in pairs_dict.keys():

        together12=set()
        only1=set()
        only2=set()
        not12=set()
        
        patients_dict=patient_mutated_residues(gene_name)  
    
        for patient in patients_dict.keys():
            if (pair[0] in patients_dict[patient]) and (pair[1] in patients_dict[patient]):
                together12.add(patient)
            elif (pair[0] in patients_dict[patient]) and (pair[1] not in patients_dict[patient]):
                only1.add(patient)
            elif (pair[0] not in patients_dict[patient]) and (pair[1] in patients_dict[patient]):
                only2.add(patient)
            elif (pair[0] not in patients_dict[patient]) and (pair[1] not in patients_dict[patient]):
                not12.add(patient)
        
        obs1 = np.array([[len(set(together12)), len(set(only2))], [len(set(only1)), len(set(not12))]])
        obs2=  np.array([[len(set(together12)), len(set(only1))], [len(set(only2)), len(set(not12))]])
        p_valu1=scipy.stats.chi2_contingency(obs1)[1]
        p_valu2=scipy.stats.chi2_contingency(obs2)[1]
        oddsratio, pvalue = stats.fisher_exact(obs1)
        if len(set(together12)) >= 2:
            print(set(together12))
            outfile2.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(pair,len(set(only1)),len(set(only2)),len(set(together12)),len(set(not12)),p_valu1, oddsratio, pvalue,together12))
            outfile3.write("{0}\t{1}\n".format(pair,set(not12)))
        if len(set(together12)) == 1:
            outfile1.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(pair,len(set(only1)),len(set(only2)),len(set(together12)),len(set(not12)),p_valu1, oddsratio, pvalue,together12))
        
#write_statistics_to_file("PIK3R1")
