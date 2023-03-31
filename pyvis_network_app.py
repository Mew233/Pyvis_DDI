import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import networkx as nx
from pyvis.network import Network
import numpy as np

# Read dataset (CSV)
# df_interact = pd.read_csv('data/processed_drug_interactions.csv')
single_drug_adr = pd.read_csv('data/Single_drug_adr.csv')
ddi = pd.read_csv('data/drug_drug.csv')

# Set header title
st.title('Network Graph Visualization of Drug-Drug Interactions')

# Define list of selection options and sort alphabetically
drug_list = list(np.unique(single_drug_adr['source']))
drug_list.sort()

# Implement multiselect dropdown menu for option selection (returns a list)
selected_drugs = st.multiselect('Select drug(s) to visualize', drug_list)

# Set info message on initial site load
if len(selected_drugs) == 0:
    st.text('Choose at least 1 drug to start')

# Create network graph when user selects >= 1 item
else:

    df_select = single_drug_adr.loc[single_drug_adr['source'].isin(selected_drugs)]
    df_select = df_select.reset_index(drop=True)

    ddi_select = ddi.loc[ddi['DRUG_1_CONCEPT_NAME'].isin(selected_drugs) | \
                                ddi['DRUG_2_CONCEPT_NAME'].isin(selected_drugs)]
    ddi_select = ddi_select.reset_index(drop=True)

    # Create networkx graph object from pandas dataframe
        # Initiate PyVis network object
    got_net = Network(
                       height='400px',
                       width='100%',
                       bgcolor='#222222',
                       font_color='white'
                      )

    # set the physics layout of the network
    got_net.barnes_hut()


    # create graph using pviz network 
    edge_data = zip(df_select['source'], df_select['target'], df_select['rel'])

    for src, dst, rel in edge_data:
        #add nodes and edges to the graph
        got_net.add_node(src, src, title=src, color='#6BAEA9')

        got_net.add_node(dst, dst, title=dst, color='#F08327')

        got_net.add_edge(src, dst, color='#CDCDCD')

    # ddi
    ddi_edge_data = zip(ddi_select['DRUG_1_CONCEPT_NAME'], ddi_select['DRUG_2_CONCEPT_NAME'], ddi_select['EVENT_CONCEPT_NAME'])

    for src_1, src_2, rel in ddi_edge_data:
        got_net.add_node(src_1, src_1, title=src_1, color='#6BAEA9')
        got_net.add_node(src_2, src_2, title=src_2, color='#6BAEA9')
        got_net.add_node(rel, rel, title=rel, color='#CDCDCD')

        got_net.add_edge(src_1, src_2, color='#7575B6',value = 6)
        got_net.add_edge(src_1, rel, color='#CDCDCD')
        got_net.add_edge(src_2, rel, color='#CDCDCD')
        


    # G = nx.from_pandas_edgelist(df_select, 'DRUG_CONCEPT_NAME', 'EVENT_CONCEPT_NAME', 'EVENT_TYPE')


    # Generate network with specific layout settings
    got_net.repulsion(
                        node_distance=420,
                        central_gravity=0.33,
                        spring_length=110,
                        spring_strength=0.10,
                        damping=0.95
                       )

    # Save and read graph as HTML file (on Streamlit Sharing)
    try:
        path = '/tmp'
        got_net.save_graph(f'{path}/pyvis_graph.html')
        HtmlFile = open(f'{path}/pyvis_graph.html', 'r', encoding='utf-8')

    # Save and read graph as HTML file (locally)
    except:
        path = '/html_files'
        got_net.save_graph(f'{path}/pyvis_graph.html')
        HtmlFile = open(f'{path}/pyvis_graph.html', 'r', encoding='utf-8')

    # Load HTML file in HTML component for display on Streamlit page
    components.html(HtmlFile.read(), height=435)

# Footer
st.markdown(
    """
    <br>
    <h6>Chengqi Xu <a href="https://github.com/Mew233" target="_blank">GitHub Repo</a></h6>
    <h6>@ElementoLab, Weill Cornell Medicine</h6>
    """, unsafe_allow_html=True
    )
