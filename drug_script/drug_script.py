import pandas as pd
import csv
import numpy as np
from pandas import DataFrame

db_miscarriage = pd.read_csv("miscarriage_database.csv")

db_miscarriage_drug=db_miscarriage[['id','type_of_miscarriage','periconceptional_drug_1','periconceptional_drug_2','periconceptional_drug_3','periconceptional_drug_4']]
#db_miscarriage_drug.to_csv('/home/gianluca/prova1.csv')

melted= pd.melt(db_miscarriage_drug, id_vars=['id','type_of_miscarriage'], value_vars =['periconceptional_drug_1','periconceptional_drug_2','periconceptional_drug_3','periconceptional_drug_4'],var_name='number_of_pcd',value_name='periconceptional_drug')

#melted.to_csv('/home/gianluca/periconceptional_drug.csv')


drug_chem_agent = csv.reader(open('drug_pa.csv', 'r'))

dictionary_chem_agent = {}

for row in drug_chem_agent:
    k, v = row
    dictionary_chem_agent[k] = v
#print(dictionary_chem_agent)


melted['chem_agent'] = melted['periconceptional_drug'].map(dictionary_chem_agent)

#melted.to_csv('/home/gianluca/periconceptional_drug.csv')

drug_smiles = csv.reader(open('drug_smile.csv', 'r'))

dictionary_smiles = {}

for row in drug_smiles:
    k, v = row
    dictionary_smiles[k]= v

melted['smiles'] = melted['periconceptional_drug'].map(dictionary_smiles)

melted_complete = melted.replace(np.nan,'NA', regex=True)

melted_complete.to_csv('/home/gianluca/draft/draft_data/periconceptional_drug.csv')



