import pandas as pd
import networkx as nx
import numpy as np

# Read dataset (CSV)
# df_interact = pd.read_csv('data/processed_drug_interactions.csv')
single_drug_adr = pd.read_csv('data/Single_drug_adr.csv',index_col=0)
ddi = pd.read_csv('data/drug_drug.csv',index_col=0)
dpi = pd.read_csv('data/drug_gene.csv')
ppi = pd.read_csv('data/gene_gene.csv')

# MICROMEDEX_SEV_LEVEL

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

dpi_select = dpi.loc[dpi['drug_node_name'].isin(selected_drugs)]
dpi_select = dpi_select.reset_index(drop=True)

selected_genes = list(dpi_select['Gene Name'])
ppi_select = ppi.loc[ppi['source'].isin(selected_genes)| \
                            ppi['target'].isin(selected_genes)]
ppi_select = ppi_select.reset_index(drop=True)

# Create networkx graph object from pandas dataframe
G = nx.from_pandas_edgelist(df_select, 'source', 'target', 'rel')

dpi_select.columns = ['source','target']
L = pd.concat([dpi_select,ppi_select])
G = nx.from_pandas_edgelist(L, 'source', 'target')
# all shortest path
print([p for p in nx.all_shortest_paths(G, source=drug_list[0], target=drug_list[1])])
print('--')

