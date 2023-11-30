# -*- coding: utf-8 -*-
"""
Module for graph building and management

Created on Thu Jun  9 10:50:12 2022

@author: Carmen Reep
"""
# All files before concatenation should have the same format:
# Edges format:
# 'subject_id',
# 'property_id',
# 'object_id',
# 'reference_uri',
# 'reference_supporting_text',
# 'reference_date',
# 'property_label',
# 'property_description',
# 'property_uri'

# Nodes format:
# 'id',
# 'semantic_groups',
# 'uri',
# 'preflabel',
# 'synonyms',
# 'description'

import pandas as pd
import os
import datetime
import requests


# VARIABLES
today = datetime.date.today()

# BUILD GRAPH

def build_edges(dgidb,monarch_disease,monarch_symptom,drugsim,input_from_file=False):
    """
    This function builds the edges graph. The user can choose to input individual networks from file or \
    from the workflow.
    :param monarch_disease: monarch disease graph edges object list
    :param monarch_symptom: monarch symptom graph edges object list
    :param dgidb: dgidb graph edges object list
    :param drugsim: drug-drug similarity object list
    :param input_from_file: False (default value) or True
    :return: edges dataframe
    """

    print('\nThe function "build_edges()" is running...')

    ## Edges

    # load networks
    print('\nPreparing networks...')
    if input_from_file:
        if isinstance(monarch_disease,str) and isinstance(monarch_symptom, str) and isinstance(dgidb,str) and isinstance(drugsim,str):
            monarch_dis_df = pd.read_csv(monarch_disease)
            monarch_symp_df = pd.read_csv(monarch_symptom)
            dgidb_df = pd.read_csv(dgidb)
            drugsim_df = pd.read_csv(drugsim)
        else:
            print("Please, if you are providing the input from file then introduce the file path to the CSV file \
            , e.g. curation=str(/home/../file_name.csv). Otherwise, provide the objects and set the 'input_from_file' \
            argument to 'False'. Thanks!")
            raise
    else:
        monarch_dis_df = pd.DataFrame(monarch_disease)
        monarch_symp_df = pd.DataFrame(monarch_symptom)
        dgidb_df = pd.DataFrame(dgidb)
        drugsim_df = pd.DataFrame(drugsim)
        
    print('Monarch:')
    print(monarch_dis_df.shape)
    print(monarch_dis_df.columns)
    print(monarch_symp_df.shape)
    print(monarch_symp_df.columns)
    
    print('DGIdb:')
    print(dgidb_df.shape)
    print(dgidb_df.columns)
    
    print('Drug similarity edges')
    print(drugsim_df.shape)
    print(drugsim_df.columns)


    # concat 1) dgidb 2) monarch 3) drugsim
    print('\nConcatenating into a graph...')
    statements = pd.concat([monarch_dis_df, monarch_symp_df, dgidb_df, drugsim_df], ignore_index=True, join="inner")
    print(statements.shape)

    # drop row duplicates
    print('\nDrop duplicated rows...')
    statements.drop_duplicates(keep='first', inplace=True)
    print(statements.shape)

    # add property_uri for those without but with a curie property_id annotated
    curie_dct = {
        'ro': 'http://purl.obolibrary.org/obo/',
        'bfo': 'http://purl.obolibrary.org/obo/',
        'geno': 'http://purl.obolibrary.org/obo/',
        'dc': 'http://purl.org/dc/elements/1.1/',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
        'skos': 'http://www.w3.org/2004/02/skos/core#',
        'pato': 'http://purl.obolibrary.org/obo/',
        'sio': 'http://semanticscience.org/resource/',
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',
        'encode': 'https://www.encodeproject.org/search/?searchTerm='
    }
    
    for i, row in statements.iterrows():
        if ':' in str(row['property_uri']):
            property_uri = row['property_uri']
        elif ':' in str(row['property_id']) and str(row['property_id']).split(':')[0].lower() == 'skos':
            property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].split(':')[1]
        elif ':' in str(row['property_id']):
            try:
                property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].replace(':',
                                                                                                                '_')
            except KeyError:
                property_uri = None
                print('There is a reference curie with and unrecognized namespace:', row['property_id'])
        else:
            property_uri = None
        statements.at[i, 'property_uri'] = property_uri

    # save graph
    print('\nSaving final graph...')
    path = os.getcwd() + "/graph"
    statements = statements[['subject_id', 'property_id', 'object_id', 'reference_uri',
                              'reference_supporting_text', 'reference_date', 'property_label',
                              'property_description', 'property_uri']]
    print(statements.shape)
    print(statements.columns)
    statements.fillna('NA').to_csv('{}/graph_edges_v{}.csv'.format(path, today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(statements.shape))
    print('* These are the edges attributes: {}'.format(statements.columns))
    print('* This is the first record:\n{}'.format(statements.head(1)))
    print('\nThe knowledge graph edges are built and saved at:'
          ' {}/graph_edges_v{}.csv\n'.format(path, today))
    print('\nFinished build_edges().\n')

    return statements


def build_nodes(statements,dgidb,monarch_disease,monarch_symptom,input_from_file=False):
    """
    This function builds the nodes graph. The user can choose to input individual networks from file or \
    from the workflow.
    :param statements: graph edges dataframe
    :param monarch_disease: monarch disease graph nodes object list
    :param monarch_symptom: monarch symptom graph nodes object list
    :param dgidb: dgidb graph nodes object list
    :param input_from_file: False (default value) or True
    :return: nodes dataframe
    """

    print('\nThe function "build_nodes()" is running...')
    # load networks
    print('\nPreparing networks...')
    if input_from_file:
        if isinstance(monarch_disease,str) and isinstance(monarch_symptom, str) and isinstance(dgidb,str):
            monarch_dis_df = pd.read_csv(monarch_disease)
            monarch_symp_df = pd.read_csv(monarch_symptom)
            dgidb_df = pd.read_csv(dgidb)

        else:
            print("Please, if you are providing the input from file then introduce the file path to the CSV file \
            , e.g. curation=str(/home/../file_name.csv). Otherwise, provide the objects and set the 'input_from_file' \
            argument to 'False'. Thanks!")
            raise
    else:
        monarch_dis_df = pd.DataFrame(monarch_disease)
        monarch_symp_df = pd.DataFrame(monarch_symptom)
        dgidb_df = pd.DataFrame(dgidb)

        
    print('Monarch:')
    print(monarch_dis_df.shape)
    print(monarch_dis_df.columns)
    print(monarch_symp_df.shape)
    print(monarch_symp_df.columns)
    
    print('DGIdb:')
    print(dgidb_df.shape)
    print(dgidb_df.columns)

    ## Annotating nodes in the graph
    print('\nAnnotating nodes in the graph...')
    # extracting nodes in the graph
    st_nodes_l = pd.concat([statements.subject_id, statements.object_id], ignore_index=True)
    st_nodes_l.drop_duplicates(inplace=True)
    st_nodes_df = pd.DataFrame({'id': st_nodes_l})
    print('graph from e', st_nodes_df.shape)

    # annotating nodes
    monarch_dis_nodes = pd.merge(monarch_dis_df, st_nodes_df, how='inner', on='id')
    monarch_symp_nodes = pd.merge(monarch_symp_df, st_nodes_df, how='inner', on='id')    
    dgidb_nodes = pd.merge(dgidb_df, st_nodes_df, how='inner', on='id')

    print('annotation check')
    print('monarch dis', monarch_dis_nodes.shape)
    print('monarch symp', monarch_symp_nodes.shape)
    print('dgidb', dgidb_nodes.shape)

    # concat all, (importantly, concatenate first curated concepts with extended definitions)
    print('\nConcatenating all nodes...')
    nodes = pd.concat([monarch_dis_nodes,monarch_symp_nodes,dgidb_nodes], ignore_index=True,
                      join="inner")
    print('graph ann', nodes.shape)
    diff = set(st_nodes_df.id) - set(nodes.id)
    print('diff', diff)

    # drop duplicated rows
    print('\nDrop duplicated rows...')
    nodes['synonyms'] = nodes.synonyms.apply(lambda x: str('|'.join(x)) if isinstance(x, list) else x)
    nodes.drop_duplicates(keep='first', inplace=True)
    print(nodes.shape)

    # drop duplicated nodes (keep first row (the curated), remove others (monarch))
    print('\nDrop duplicated nodes...')
    nodes.drop_duplicates(subset=['id'], keep='first', inplace=True)
    print(nodes.shape)

    # check
    if len(set(st_nodes_df.id)) != len(set(nodes.id)):
        print(
            '\nThere is a problem in the annotation of nodes.\nThe number of annotated nodes '
            'is different than the number of nodes in the graph.')
        print('Monarch disease nodes not in the graph: {}'.format(set(monarch_dis_df.id) - set(monarch_dis_nodes.id)))
        print('Monarch symptom nodes not in the graph: {}'.format(set(monarch_symp_df.id) - set(monarch_symp_nodes.id)))
        print('DGIdb nodes not in the graph: {}'.format(set(dgidb_df.id) - set(dgidb_nodes.id)))

    else:
        print('\nAll graph nodes are annotated.')

    # save graph nodes
    print('\nSaving final graph...')
    path = os.getcwd() + "/graph"
    nodes = nodes[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    nodes['synonyms'] = nodes.synonyms.apply(lambda x: str('|'.join(x)) if isinstance(x, list) else x)
    print(nodes.shape)
    print(nodes.columns)
    nodes.fillna('NA').to_csv('{}/graph_nodes_v{}.csv'.format(path, today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(nodes.shape))
    print('* These are the edges attributes: {}'.format(nodes.columns))
    print('* This is the first record:\n{}'.format(nodes.head(1)))
    print('\nThe knowledge graph nodes are built and saved at: '
          '{}/graph_nodes_v{}.csv\n'.format(path, today))
    print('\nFinished build_nodes().\n')

    return nodes

def run_combine_graphs(date):
    """
    This function runs the combine_graphs script and saves nodes and edges files.
    :param date: the date of creation of the disease graph
    :return: nodes and edges files in /graph folder
    """
    # path to write data
    path = os.getcwd() + "/graph"
    if not os.path.isdir(path): os.makedirs(path)

    drugsim_file = './similaritygraph/drugdrugsim_v{}.csv'.format(today)
    dgidb_edges_file = './DGIdb/DGIdb_edges_v{}.csv'.format(today)
    dgidb_nodes_file = './DGIdb/DGIdb_nodes_v{}.csv'.format(today)
    monarch_edges_disease_file = './monarch/monarch_edges_disease_v{}.csv'.format(date)
    monarch_nodes_disease_file = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    monarch_edges_symptom_file = './monarch/monarch_edges_symptom_v{}.csv'.format(today)
    monarch_nodes_symptom_file = './monarch/monarch_nodes_symptom_v{}.csv'.format(today)
    
    # build network
    statements = build_edges(dgidb=dgidb_edges_file,monarch_disease=monarch_edges_disease_file,monarch_symptom=monarch_edges_symptom_file,drugsim=drugsim_file, input_from_file=True)
    build_nodes(statements, dgidb=dgidb_nodes_file,monarch_disease=monarch_nodes_disease_file,monarch_symptom=monarch_nodes_symptom_file,input_from_file=True)
  

if __name__ == '__main__':
    run_combine_graphs()

      