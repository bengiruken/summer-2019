#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 18:30:36 2019

@author: bengi
"""
'''aşağıdaki kısım oncogene-tumor supressor- maf dosyasını entegre edip dosyalara yazdırma'''
#''' ''' arasında kalan kısım başka dosyalara yazılabilir
'''
def create_folder():
    import os
    dirName = 'Tumor_Supp_Oncogene'
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
#cancer_gene_sensus_tumor_supressor ve oncogene olanları alam
def cancer_gene_sensus_as_DataFrame():
    import pandas as pd
    df1 = pd.read_csv("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Data/cancer_gene_census.csv", usecols=['Gene Symbol', 'Role in Cancer'])
    df = df1[df1['Role in Cancer'].isin(['TSG', 'oncogene']) ]
    return df


df=cancer_gene_sensus_as_DataFrame()
#write these genes to a file
oncogene = df[df['Role in Cancer'].isin(['oncogene'])]
oncogene['Gene Symbol'].to_csv("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/oncogene.txt", header=None, index=None, sep=' ', mode='a')

tumor_supressor = df[df['Role in Cancer'].isin(['TSG'])]
tumor_supressor['Gene Symbol'].to_csv("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/tumor_supressor.txt", header=None, index=None, sep=' ', mode='a')
tumor_supressor_list = tumor_supressor['Gene Symbol'].to_list()
   tumor_supressor_tuple = tuple(tumor_supressor['Gene Symbol'].to_list())

############


#####

#header=header_of_maf_file("/Users/bengi/Desktop/pik3ca/GitHub_codes_data_pik3ca/Data/mc3.v0.2.8.PUBLIC.maf" )
#tumor supressor'da mutasyonu olanları dosyaya yazma
def tumor_supressor_mutation_file():
    outfile=open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/maf_tumor_supressor_mutations.txt', 'a')
    #not_exclude_mutations=["5'UTR","Frame_Shift_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
    #exclude_mutations = [".",'RNA',"Silent"]
    #exclude_from_parsed = ['.','?','_']
    tumor_supressor_tuple = tuple(tumor_supressor['Gene Symbol'].to_list())
    with open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/mc3.v0.2.8.PUBLIC.maf","r") as infile:
        for line in infile:
            if line.startswith(tumor_supressor_tuple):
                outfile.write(line)
    outfile.close() 
    
tumor_supressor_mutation_file()    

####################

#oncogene'de mutasyonu olanları dosyaya yazma
def oncogene_mutation_file():
    outfile=open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/maf_oncogene_mutations.txt', 'a')
    #not_exclude_mutations=["5'UTR","Frame_Shift_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
    #exclude_mutations = [".",'RNA',"Silent"]
    #exclude_from_parsed = ['.','?','_']
    oncogene_tuple = tuple(oncogene['Gene Symbol'].to_list())
    with open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/mc3.v0.2.8.PUBLIC.maf","r") as infile:
        for line in infile:
            if line.startswith(oncogene_tuple):
                outfile.write(line)
    outfile.close() 
    
oncogene_mutation_file()  

####################

from Python_Codes.HGVNp_column_parsing import HGVNp_short_parse
import Python_Codes.read_maf_file_for_certain_gene as y
import Python_Codes.tumor_patient_id as z

#header=header_of_maf_file("/Users/bengi/Desktop/pik3ca/GitHub_codes_data_pik3ca/Data/mc3.v0.2.8.PUBLIC.maf" )
def tumor_supressor_mutation_file():
    from Python_Codes.HGVNp_column_parsing import HGVNp_short_parse

    outfile=open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/tsg_not_silent.txt", 'a')
    problems=open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/tumor_supressor_problematic_lines.txt', 'a') 
    not_exclude_mutations=["5'UTR","Frame_Shift_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
    exclude_mutations=["3'UTR","3'Flank",".","5'Flank","In_Frame_Del","In_Frame_Ins",'Intron','RNA',"Silent","Splice_Site","Translation_Start_Site"]
    exclude_from_parsed=['.','?','_']
    with open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/maf_tumor_supressor_mutations.txt','r') as infile:
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

TS = tumor_supressor_mutation_file()


##################

def oncogene_mutation_file():
    from Python_Codes.HGVNp_column_parsing import HGVNp_short_parse

    outfile=open("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/oncogene_not_silent.txt", 'a')
    problems=open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/oncogene_problematic_lines.txt', 'a') 
    not_exclude_mutations=["5'UTR","Frame_Shift_Del","In_Frame_Del","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"]
    exclude_mutations=["3'UTR","3'Flank",".","5'Flank","In_Frame_Del","In_Frame_Ins",'Intron','RNA',"Silent","Splice_Site","Translation_Start_Site"]
    exclude_from_parsed=['.','?','_']
    with open('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/maf_oncogene_mutations.txt','r') as infile:
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

ONG=oncogene_mutation_file()
'''



############
#read tsg and oncogene not silent mutations for observation of commonalities and differences
import pandas as pd
header=['Hugo_Symbol','Patient_ID','Tumor_ID','Protein_Change_HGVSp_Short','Residue#','Residue','New_Residue','Variant_Classification']

TSG_df = pd.read_csv("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/tsg_not_silent.txt", sep = "\t", header= None)
TSG_df.columns = header
TSG_patient_list = list(set(TSG_df['Patient_ID'].to_list()))


oncogene_df = pd.read_csv("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/oncogene_not_silent.txt", sep = "\t", header=None)
oncogene_df.columns = header
oncogene_patient_list = list(set(oncogene_df['Patient_ID'].to_list()))

#find differences of tsg&oncogene mutated patients
def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 

only_tumor_supressor =  Diff(TSG_patient_list, oncogene_patient_list)
only_oncogene =  Diff(oncogene_patient_list, TSG_patient_list)

##intersection of bot tsg&oncogene mutated
def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 

both_tsg_oncogene= intersection(TSG_patient_list, oncogene_patient_list )
##########
#oncogene vs tsg survival plot 
import Python_Codes.maf_clinical_data_processing as clinical
df_survival = clinical.survival_data()

def survival(cancer_type,days,patient_list_1,patient_list_2, name1="oncogene",name2="tumor_supressor",alpha=99):
    !pip install lifelines
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test as multiv_lr_test
    title = 'KM Survival oncogene vs tumor supressor {} days'.format(days)
    # Initialize KM plotter
    kmf = KaplanMeierFitter()
    #load cluster file
    #survival bilgisini uygun başlıklarla değiştrime
    subgroup1=["bcr_patient_barcode","vital_status","death_days_to","last_contact_days_to","overall_survival"]
    survival_info = df_survival
    survival_info.columns=subgroup1
    surv = survival_info[subgroup1]
    surv = surv.dropna(subset = [ 'vital_status']) 
    
    # Number of clusters
    clusters = [1,2]
    k = len(clusters)
    # Initialize KM Plot Settings
    fig = plt.figure(figsize=(10, 7)) 
    ax = plt.subplot(111)
    colors = sns.color_palette('hls', k)
    cluster_cmap = {clusters[i]:colors[i] for i in range(k)}
    
    # Plot each cluster onto KM Plot
    #only patient_group_1 
    # only_oncogene
    filter1 = surv["bcr_patient_barcode"].isin(patient_list_1) 
    clust_surv_data_1 = surv[filter1]
    clust_surv_data_1=clust_surv_data_1.dropna(subset = ["bcr_patient_barcode","vital_status","death_days_to","last_contact_days_to","overall_survival"]) 
    kmf.fit(clust_surv_data_1.overall_survival, clust_surv_data_1.vital_status, label='{}'.format(name1)+' (n=' +  str(len(clust_surv_data_1)) + ')')
    kmf.plot(ax=ax, color=cluster_cmap[1], ci_show=False)
    # Set KM plot limits to 5 years and labels
    
    # Plot each cluster onto KM Plot
    #pik3ca'deki pairları içerenler
    #only_tumor_supressor
    
    filter2 = surv["bcr_patient_barcode"].isin(patient_list_2) 
    clust_surv_data_2 = surv[filter2]
    clust_surv_data_2=clust_surv_data_2.dropna(subset = ["bcr_patient_barcode","vital_status","death_days_to","last_contact_days_to","overall_survival"]) 
    
    kmf.fit(clust_surv_data_2.overall_survival, clust_surv_data_2.vital_status, label='{} '.format(name2)+' (n=' +  str(len(clust_surv_data_2)) + ')')
    kmf.plot(ax=ax, color=cluster_cmap[2], ci_show=False)
    # Plot each cluster onto KM Plot
    #pik3ca'deki pairları içerenler
    from lifelines.statistics import logrank_test
    summary= logrank_test(clust_surv_data_1.overall_survival, clust_surv_data_2.overall_survival, clust_surv_data_1.vital_status, clust_surv_data_2.vital_status, alpha)

    plt.xlim((0,days))
    plt.xlabel('Time_Span (Days) = {0}  '.format(days), fontsize=16) #, p_value = {1:.3g} format'ın içine,summary.p_value
    plt.ylabel('Survival Probability', fontsize=16)
    plt.savefig('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/{}_oncogene_vs_tumor_supressor_{}days.png'.format(cancer_type,days),bbox_inches='tight')
    plt.show()
    
    return summary.p_value, summary.test_statistic

 #end of function   
    
    
survival(4500, only_oncogene, only_tumor_supressor)    
#print("p_value: {}, test_statistic: {}".format(summary.p_value,summary.test_statistic) )   


survival(5000, only_oncogene, only_tumor_supressor)        

###########3

#tumor grade and tumor section information

#import Python_Codes.maf_clinical_data_processing as clinical
#df_tumor_type_stage = clinical.tumor_type_stage_as_DataFrame()

def tumor_type_stage_as_DataFrame():
    import pandas as pd
    xls = pd.ExcelFile("/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Data/NIHMS978596-supplement-1.xlsx")
    df1 = pd.read_excel(xls, sheet_name="TCGA-CDR")
#    col_names=list(df1.columns) 
#    col_names
    subgroup=['bcr_patient_barcode', 
          'type'
          ]

    data= df1[subgroup]
    return data


df_tumor_type=tumor_type_stage_as_DataFrame()
#oncogene
filter_oncogene =  df_tumor_type["bcr_patient_barcode"].isin(only_oncogene) 
df_oncogene_tumor_type =  df_tumor_type[filter_oncogene]
count_oncogene=df_oncogene_tumor_type.groupby(["type"]).count()
count_oncogene.to_csv('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/tumor_type_oncogene.csv', sep=",")
#
#oncogene
filter_tsg =  df_tumor_type["bcr_patient_barcode"].isin(only_tumor_supressor) 
df_tumor_supressor_tumor_type =  df_tumor_type[filter_tsg]
count_tumor_supressor=df_tumor_supressor_tumor_type.groupby(["type"]).count()
count_tumor_supressor.to_csv('/Users/bengi/Desktop/pik3ca/GitHub_codes_data/Tumor_Supp_Oncogene/tumor_type_tumor_supressor.csv', sep=",")

################survival plots for tumor type
df_onco_type_patient=df_oncogene_tumor_type.groupby(["type"])["bcr_patient_barcode"].apply(list)
type(df_onco_type_patient["BRCA"])
df_tsg_type_patient=df_tumor_supressor_tumor_type.groupby(["type"])["bcr_patient_barcode"].apply(list)

survival("PRAD",5000,df_onco_type_patient["PRAD"],df_tsg_type_patient["PRAD"], "PRAD_oncogene" , "PRAD_tumor_supressor",95)
survival("BRCA",5000,df_onco_type_patient["BRCA"],df_tsg_type_patient["BRCA"], "BRCA_oncogene" , "BRCA_tumor_supressor",95)
survival("SARC",5000,df_onco_type_patient["SARC"],df_tsg_type_patient["SARC"], "SARC_oncogene" , "SARC_tumor_supressor",95)
survival("GBM",5000,df_onco_type_patient["GBM"],df_tsg_type_patient["GBM"], "GBM_oncogene" , "GBM_tumor_supressor",95)



####her bir tumor type için 0-1 vital status'a sahip hasta sayılarını karşılaştırma
#returns: clinical data'dan sadece o kanser tipinde hastların data framei
#vitsl statusu 0 olan hasta sayısı
#vitsl statusu 1 olan hasta sayısı

def vital_status_numbers(patient_list_1):
    import pandas as pd
    import numpy as np
    import matplotlib
    import Python_Codes.maf_clinical_data_processing as clinical
    df_survival = clinical.survival_data()
    df_survival = df_survival.dropna(subset = [ 'vital_status'])
    filter1 = df_survival["bcr_patient_barcode"].isin(patient_list_1) 
    cluster1= df_survival[filter1]
    #şu olmayınca sayılar farklı çıkıyor
    #hangi row'da Nan varsa onu atıyoruz, hem km plot, hem burası için
    cluster1=cluster1.dropna() 

    #count_cluster1 = df_survival.groupby(["vital_status"])["bcr_patient_barcode"].count()
    return cluster1, len(cluster1[cluster1["vital_status"]==0]), len(cluster1[cluster1["vital_status"]==1])
#end of function

prad_oncogene = vital_status_numbers(df_onco_type_patient["PRAD"])
prad_oncogene[1]
prad_oncogene[2]
prad_tsg= vital_status_numbers(df_tsg_type_patient["PRAD"])
prad_tsg[1]
prad_tsg[2]

brca_oncogene = vital_status_numbers(df_onco_type_patient["BRCA"])
brca_oncogene[1]
brca_oncogene[2]
brca_tsg= vital_status_numbers(df_tsg_type_patient["BRCA"])
brca_tsg[1]
brca_tsg[2]
len(brca_tsg[0])
