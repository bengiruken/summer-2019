#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 17:53:02 2019

@author: bengi
"""

#create folder for each gene name
# Create directory
def create_folder(gene_name):
    import os
    dirName = '{}'.format(gene_name)
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

#create_folder("PIK3R1")
#create_folder("PIK3R1")

'''pan_cancer data file read mutatnts of specific genes'''
#cd '/Users/bengi/Desktop/pik3ca'
def gene_mutants_from_pan_cancer_data(gene_name):
    gene_file = open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/1_{}.txt".format(gene_name,gene_name),"w")
    with open( "/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Data/mc3.v0.2.8.PUBLIC.maf" , "r") as infile:
        for lines in infile:
            if lines.startswith("{}".format(gene_name)):
                gene_file.write(lines)
                
#gene_mutants_from_pan_cancer_data("PIK3R1")
'''header of the maf_file, the first line of the file,output=list'''
def header_of_maf_file(maf_file):
    with open(maf_file) as f:
        first_line = f.readline()
    maf_file_header = first_line.split("\t")   
    return maf_file_header

#header=header_of_maf_file("/Users/bengi/Desktop/pik3ca/GitHub_codes_data_pik3ca/Data/mc3.v0.2.8.PUBLIC.maf" )
def gene_specific_mutation_file(gene_name):
    from Python_Codes.HGVNp_column_parsing import HGVNp_short_parse
    outfile=open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/2_{}_point_mutations.txt'.format(gene_name,gene_name), 'a')
    problems=open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/0_{}_problematic_lines.txt'.format(gene_name,gene_name), 'a') 
    not_exclude_mutations=["5'UTR","Frame_Shift_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
    exclude_mutations=["3'UTR","3'Flank",".","5'Flank","In_Frame_Del","In_Frame_Ins",'Intron','RNA',"Silent","Splice_Site","Translation_Start_Site"]
    exclude_from_parsed=['.','?','_']
    with open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/1_{}.txt".format(gene_name,gene_name),'r') as infile:
        for line in infile:
            splitted = line.split("\t")
            if (splitted[15][13:15]=='01'):
                #print(splitted[15][:15])
                #splitted[36]: variant classification
                if (len(splitted)>1) and (splitted[36]!='.') and (splitted[8] not in exclude_mutations) and (splitted[36]!='p.%3D'):
                    if (splitted[8]=="Frame_Shift_Del") or (splitted[8]=="Frame_Shift_Ins"):
                        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(splitted[0],splitted[15][:12],splitted[15][:15],splitted[36],HGVNp_short_parse(splitted[36])[0],HGVNp_short_parse(splitted[36])[1],HGVNp_short_parse(splitted[36])[2].split('fs*')[0],splitted[8]))
                    elif (splitted[8]=="Nonstop_Mutation"):
                        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(splitted[0],splitted[15][:12],splitted[15][:15],splitted[36],HGVNp_short_parse(splitted[36])[0],HGVNp_short_parse(splitted[36])[1],HGVNp_short_parse(splitted[36])[2].split('ext*')[0],splitted[8]))
                    elif (splitted[8]=="Nonsense_Mutation"):
                        if '_' not in HGVNp_short_parse(splitted[36])[2]:
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(splitted[0],splitted[15][:12],splitted[15][:15],splitted[36],HGVNp_short_parse(splitted[36])[0],HGVNp_short_parse(splitted[36])[1],HGVNp_short_parse(splitted[36])[2],splitted[8]))
                    elif (splitted[8]=="Missense_Mutation"):
                        if ('_' not in HGVNp_short_parse(splitted[36])[2]) and (HGVNp_short_parse(splitted[36])[1] not in exclude_from_parsed) and (HGVNp_short_parse(splitted[36])[2] not in exclude_from_parsed) :
                            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(splitted[0],splitted[15][:12],splitted[15][:15],splitted[36],HGVNp_short_parse(splitted[36])[0],HGVNp_short_parse(splitted[36])[1],HGVNp_short_parse(splitted[36])[2],splitted[8]))
                    else: 
                        problems.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(splitted[0],splitted[15][:12],splitted[15][:15],splitted[36],HGVNp_short_parse(splitted[36])[0],HGVNp_short_parse(splitted[36])[1],HGVNp_short_parse(splitted[36])[2],splitted[8]))
    outfile.close() 

#gene_specific_mutation_file("PIK3R1")
#Add header to the existing file 
# output:   '3_{}_point_mutations_with_header.txt'.format(gene_name)
#gene_specific_mutation_file(gene_name)    
def add_header_gene_specific_mutation_file(gene_name):
    header=['Hugo_Symbol','Patient_ID','Tumor_ID','Protein_Change_HGVSp_Short','Residue#','Residue','New_Residue','Variant_Classification']
    with open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/2_{}_point_mutations.txt'.format(gene_name,gene_name), 'r') as original: data = original.read()
    with open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/{}/3_{}_point_mutations_with_header.txt'.format(gene_name,gene_name), 'w') as modified: modified.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(header[0],header[1],header[2],header[3],header[4],header[5],header[6],header[7])+data)

#add_header_gene_specific_mutation_file("PIK3R1")
