#!/home/adriano/anaconda3/bin/python3
import pandas as pd 
import numpy as np
import os, sys, argparse

### DEFAULT VALUES AS MODULE ##############################################################################
# About output
output_matrix_sep = '\t'
output_matrixes_encoding = 'utf-8'
output_na_rep = ''  # how to WRITE missing data
# About input
input_matrixes_sep = '\t'
input_matrixes_encoding = 'utf-8'
#color bash
g = '\033[92m'
e = '\033[0m'
w = '\033[93m'
r = '\033[91m'
###########################################################################################################

header = """
+------------------------------------------------------+
             ***  Merging genes list  ***
+------------------------------------------------------+
 Author:    Adriano De Marino
 Date:      January 2020
 Contact:   adriano.demarino@igenomix.com
+------------------------------------------------------+
""".format(g=g,e=e)

description = """ 

Description:
This program takes in input different list of genes and give in output the merging without duplicated rows

"""
print(header)

usage_example = """
 {g}-f1{e}   input file1 required
 {g}-f2{e}   input file2 required
 {g}-f3{e}   input file3 optional
 {g}-f4{e}   input file4 optional
 {g}-fn{e}   input filen optional

Examples of usage: {w} python merge_list_of_genes.py -f1 /path_to_file/filename1.list -f2 /path_to_file/filename2.list  {e}


""".format(g=g,e=e,w=w)

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ Igenomix - Bioinformatics Italy ] \n", description = description)
parser.add_argument('-f1', dest="file1", help="required file1", action="store", required=True)
parser.add_argument('-f2', dest="file2", help="required file2", action="store", required=True)
parser.add_argument("-f3", dest='file3', help="optional file3", action="store", required=False, default="")
parser.add_argument("-f4", dest='file4', help="optional file4", action="store", required=False, default="")
parser.add_argument("-f5", dest='file5', help="optional file5", action="store", required=False, default="")
parser.add_argument("-f6", dest='file6', help="optional file6", action="store", required=False, default="")
parser.add_argument("-f7", dest='file7', help="optional file7", action="store", required=False, default="")
parser.add_argument("-f8", dest='file8', help="optional file8", action="store", required=False, default="")

args = parser.parse_args()

file1 = args.file1
file2 = args.file2
file3 = args.file3
file4 = args.file4
file5 = args.file5
file6 = args.file6
file7 = args.file7
file8 = args.file8

#file1 = "gene_list.bed"
#file2 = "gene_list1.bed"
#file3 = "gene_list2.bed"
#file4 = ""
#file5 = ""
#file6 = ""
#file7 = ""
#file8 = ""

#create a list of files
list_file = [file1,file2,file3,file4,file5,file6,file7,file8]
#initialize new dataframe
df = pd.DataFrame()
#read each file and append to create a unique dataframe
for file in list_file:
	if file != "":
		df = df.append(pd.read_csv(file,sep=output_matrix_sep,header=None,encoding=output_matrixes_encoding))
#remove duplicates
df.drop_duplicates(inplace=True)
#add columns name
df.columns = ["chr","start","end","ensID","gene_type"]

#save file
df.to_csv("merged_output_list.bed",sep=output_matrix_sep,encoding=output_matrixes_encoding)

#quit
exit()