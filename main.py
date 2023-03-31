import pandas as pd
import networkx as nx
import numpy as np

# Read dataset (CSV)
# df_interact = pd.read_csv('data/processed_drug_interactions.csv')
single_drug_adr = pd.read_csv('data/Single_drug_adr.csv',index_col=0)
ddi = pd.read_csv('data/drug_drug.csv')

# Define list of selection options and sort alphabetically
drug_list = list(np.unique(single_drug_adr['source']))
drug_list.sort()

# Implement multiselect dropdown menu for option selection (returns a list)
selected_drugs = [drug_list[0],drug_list[1]]

# Create network graph when user selects >= 1 item

df_select = single_drug_adr.loc[single_drug_adr['source'].isin(selected_drugs)]
df_select = df_select.reset_index(drop=True)

ddi_select = ddi.loc[ddi['DRUG_1_CONCEPT_NAME'].isin(selected_drugs) | \
                            ddi['DRUG_2_CONCEPT_NAME'].isin(selected_drugs)]
ddi_select = ddi_select.reset_index(drop=True)

# Create networkx graph object from pandas dataframe
G = nx.from_pandas_edgelist(df_select, 'DRUG_CONCEPT_NAME', 'EVENT_CONCEPT_NAME', 'EVENT_TYPE')