import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import networkx as nx
from pyvis.network import Network
import numpy as np
from utilitis import *

# Read dataset (CSV)
# df_interact = pd.read_csv('data/processed_drug_interactions.csv')
@st.experimental_memo
def load_data():
    single_drug_adr = pd.read_csv('data/Single_drug_adr.csv',index_col=0)
    ddi = pd.read_csv('data/drug_drug.csv',index_col=0)
    dpi = pd.read_csv('data/drug_gene_all.csv')
    expert = pd.read_csv('data/expert_list.csv',index_col=0)
    ppi = pd.read_csv('data/gene_gene.csv')

    ncbi2name = pd.read_csv("data/ncbi2name.txt",sep='\t')
    expert = pd.merge(expert, ncbi2name, left_on=['NCBI_ID'], right_on=['NCBI Gene ID(supplied by NCBI)'])
    expert2 = expert[['NCBI_ID','Drug IDs','Drug IDs','Drug IDs','Approved symbol']]
    expert2.columns = dpi.columns
    dpi_L = pd.concat([dpi,expert2])
    dpi_L = dpi_L.drop_duplicates(subset=['NCBI_ID','DrugBank ID_split'])
    dpi_L = dpi_L[['drug_node_name','Gene Name']]

    return single_drug_adr, ddi, dpi_L, ppi

single_drug_adr, ddi, dpi, ppi = load_data()


# Set header title
st.title('Network Graph Visualization of Drug-Drug Interactions')

# Define list of selection options and sort alphabetically
drug_list = list(np.unique(single_drug_adr['source'])) + list(set(dpi['drug_node_name'])-\
                                                                      set(single_drug_adr['source']).intersection(set(dpi['drug_node_name'])))
# drug_list.sort()

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

    dpi_select = dpi.loc[dpi['drug_node_name'].isin(selected_drugs)]
    dpi_select = dpi_select.reset_index(drop=True)

    selected_genes = list(dpi_select['Gene Name'])
    ppi_select = ppi.loc[ppi['source'].isin(selected_genes)| \
                                ppi['target'].isin(selected_genes)]
    ppi_select = ppi_select.reset_index(drop=True)

    # Create networkx graph object from pandas dataframe
        # Initiate PyVis network object
    got_net = Network(
                       height='750px',
                       width='100%',
                       bgcolor='#222222',
                       font_color='white'
                      )

    # set the physics layout of the network
    got_net.barnes_hut()
    
    #有些drug没有side effect的信息
    if len(df_select)>0 and len(ddi_select)>0:
        add_drug_side(df_select,got_net)
        add_drug_drugddi_select(ddi_select,got_net)

    #dpi
    dpi_edge_data = zip(dpi_select['drug_node_name'], dpi_select['Gene Name'])

    for src, dst in dpi_edge_data: #gene_name, 
        #add nodes and edges to the graph
        #drug
        got_net.add_node(src, title=src, color='#6BAEA9', shape='triangle',labelHighlightBold=True)
        #protein
        got_net.add_node(dst, title=dst, color='#A0AA9B', shape='star')

        got_net.add_edge(src, dst, color='#CDCDCD')


    # Generate network with specific layout settings
    got_net.repulsion(
                        node_distance=420,
                        central_gravity=0.33,
                        spring_length=110,
                        spring_strength=0.10,
                        damping=0.95
                       )
    got_net.show_buttons(filter_=['physics'])


    
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

    # display all shortest path
    # Create networkx graph object from pandas dataframe
    genre = st.radio("Meta path exploration (please select 2 drugs)",('Genes','Side effects',))
    genre_hop = st.radio("PPI: ",('shortest path', '2-hop'),horizontal=True)

    agree = st.checkbox('Show PPI network',value=True)
    

    if len(selected_drugs) >= 2:
        if genre == 'Side effects':
            # st.write('You selected comedy.')
            G = nx.from_pandas_edgelist(df_select, 'source', 'target', 'rel')
            paths = nx.all_shortest_paths(G, source=selected_drugs[0], target=selected_drugs[1])
        elif genre == 'Genes':
            # G = nx.from_pandas_edgelist(dpi_select, 'drug_node_name', 'Gene Name')
            dpi_select.columns = ['source','target']
            L = pd.concat([dpi_select,ppi_select])
            G = nx.from_pandas_edgelist(L, 'source', 'target')
            path_1 = nx.all_shortest_paths(G, source=selected_drugs[0], target=selected_drugs[1])
            path_2 = nx.all_simple_paths(G, source=selected_drugs[0], target=selected_drugs[1],cutoff=3)
            if genre_hop == 'shortest path':
                paths = path_1
            elif genre_hop == '2-hop':
                paths = path_2

        ps = [p for p in paths]
        ps = sorted(ps, key=len, reverse=True)

        st.experimental_set_query_params(my_saved_result=ps)  # Save value
        # Retrieve app state
        app_state = st.experimental_get_query_params()  

        try:
            saved_result = app_state["my_saved_result"]
            # print([p for p in saved_result])
            # print(len(saved_result))
            st.write(saved_result)
        except:
            raise ValueError('There is no 2-hop genes between selected drugs.')

        if agree:
            create_ggi(dpi_select, G,selected_drugs)


# Footer
st.markdown(
    """
    <br>
    <h6>Chengqi Xu <a href="https://github.com/Mew233" target="_blank">GitHub Repo</a></h6>
    <h6>@ElementoLab, Weill Cornell Medicine</h6>
    """, unsafe_allow_html=True
    )
