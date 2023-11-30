
# -*- coding: utf-8 -*-
"""
Module for DGIdb network preparation and management

Created on Wed Apr  6 10:49:08 2022

@author: Carmen Reep
"""

import requests
import datetime
import pandas as pd
from biothings_client import get_client
import ast
import os

# timestamp
today = datetime.date.today()


def get_genes(nodes_df):
    """
    This functions finds all genes in Monarch. It returns a dict of id, name pairs of the gene IDs  
    :param monarch_nodes_csv: csv file of monarch nodes
    :return: dict of gene ids/name
    """
    print('get_genes() is running')
    
    # keep rows with 'GENES' as semantic label
    df_genes = nodes_df.loc[nodes_df['semantic_groups'] == 'gene']    
    #get names and ids of genes
    gene_name_list = df_genes['preflabel'].to_list()
    #dict of id and names
    gene_id_list = df_genes['id'].to_list()
    id_name_dict = dict(zip(gene_id_list, gene_name_list))
    
    return id_name_dict


def normalize_genes_to_graph(id_name_dict):
    """
    This function normalizes gene ID scheme. It performs Gene ID conversion \
    from curated to graph scheme. It replaces HGNC ID by human entrez, MGI ID by mouse entrez, WormBase ID \
    by worm entrez.
    :param id_name_dict: dict of gene ids/name to convert
    :return: concept dictionary (where keys are old ontology IDs and values are entrez IDs)
    """
    
    print('\n Mapping Monarch gene ontologies to entrez genes ...')
    # final dict
    concept_dct = dict()

    # input
    zfin = list()
    ensembl = list()
    hgnc = list()
    mgi = list()
    wormbase = list()
    flybase = list()
    for idx, name in id_name_dict.items():
        if ':' in idx:
            if 'flybase' in idx.split(':')[0].lower():
                flybase.append(idx.split(':')[1])
            elif 'wormbase' in idx.split(':')[0].lower():
                wormbase.append(idx.split(':')[1])
            elif 'mgi' in idx.split(':')[0].lower():
                mgi.append(idx.split(':')[1])
            elif 'hgnc' in idx.split(':')[0].lower():
                hgnc.append(idx.split(':')[1])
            elif 'ensembl' in idx.split(':')[0].lower():
                ensembl.append(idx.split(':')[1])
            elif 'zfin' in idx.split(':')[0].lower():
                zfin.append(idx.split(':')[1])
    zfin = list(set(zfin))
    ensembl = list(set(ensembl))
    hgnc = list(set(hgnc))
    mgi = list(set(mgi))
    mgi2=['MGI:'+s for s in mgi] # for api 'MGI:' is needed infront of ID
    wormbase = list(set(wormbase))
    flybase = list(set(flybase))
    # api call
    mg = get_client('gene')

    try:    
        df_mgi = mg.querymany(qterms=mgi2, scopes=['mgi'], fields='entrezgene', size=1, as_dataframe=True) #117 outof130 notfound
    except:
        df_mgi = pd.DataFrame([])
    try:
        df_flybase = mg.querymany(qterms=flybase, scopes=['FLYBASE'], fields='entrezgene', size=1, as_dataframe=True) # 2 outof3 not found
    except:
        df_flybase = pd.DataFrame([])
    try:  
        df_hgnc = mg.querymany(qterms=hgnc, scopes=['HGNC'], fields='entrezgene', size=1, as_dataframe=True)
    except:
        df_hgnc = pd.DataFrame([])
    try:  
        df_wormbase = mg.querymany(qterms=wormbase, scopes=['WormBase'], fields='entrezgene', size=1, as_dataframe=True)
    except:
        df_wormbase = pd.DataFrame([])
    try:  
        df_zfin = mg.querymany(qterms=zfin, scopes=['ZFIN'], fields='entrezgene', size=1, as_dataframe=True) #notfound
    except:
        df_zfin = pd.DataFrame([])
    try:  
        df_ensembl = mg.querymany(qterms=ensembl, scopes=['ensembl.gene'], fields='entrezgene', size=1, as_dataframe=True)
    except:
        df_ensembl = pd.DataFrame([])
    
    #build dictionaries
    ids_mgi = df_mgi.reset_index().copy()
    if 'entrezgene' in ids_mgi.columns:
        mgi2entrez = dict(zip(ids_mgi['query'], ids_mgi['entrezgene']))
        concept_dct.update(mgi2entrez)
    else:
        print('no MGI IDs can be mapped to entrez IDs')

    ids_flybase = df_flybase.reset_index().copy()
    if 'entrezgene' in ids_flybase.columns:
        flybase2entrez = dict(zip(ids_flybase['query'], ids_flybase['entrezgene']))
        for old_key in list(flybase2entrez):
            flybase2entrez['FlyBase:'+old_key] = flybase2entrez.pop(old_key) # make sure name is infront (e.g. FlyBase:)
        concept_dct.update(flybase2entrez)
    else:
        print('no flybase IDs can be mapped to entrez IDs')

    ids_hgnc = df_hgnc.reset_index().copy()
    if 'entrezgene' in ids_hgnc.columns:
        hgnc2entrez = dict(zip(ids_hgnc['query'], ids_hgnc['entrezgene']))
        for old_key in list(hgnc2entrez):
            hgnc2entrez['HGNC:'+old_key] = hgnc2entrez.pop(old_key) # make sure name is infront
        concept_dct.update(hgnc2entrez)
    else:
        print('no HGNC IDs can be mapped to entrez IDs')

    ids_wormbase = df_wormbase.reset_index().copy()
    if 'entrezgene' in ids_wormbase.columns:
        wormbase2entrez = dict(zip(ids_wormbase['query'], ids_wormbase['entrezgene']))
        for old_key in list(wormbase2entrez):
            wormbase2entrez['WormBase:'+old_key] = wormbase2entrez.pop(old_key) # make sure name is infront
        concept_dct.update(wormbase2entrez)
    else:
        print('no WormBase IDs can be mapped to entrez IDs')

    ids_zfin = df_zfin.reset_index().copy()
    if 'entrezgene' in ids_zfin.columns:
        zfin2entrez = dict(zip(ids_zfin['query'], ids_zfin['entrezgene']))
        for old_key in list(zfin2entrez):
            zfin2entrez['ZFIN:'+old_key] = zfin2entrez.pop(old_key) # make sure name is infront
        concept_dct.update(zfin2entrez)
    else:
        print('no ZFIN IDs can be mapped to entrez IDs')

    ids_ensembl = df_ensembl.reset_index().copy()
    if 'entrezgene' in ids_ensembl.columns:
        ensembl2entrez = dict(zip(ids_ensembl['query'], ids_ensembl['entrezgene']))
        for old_key in list(ensembl2entrez):
            ensembl2entrez['ENSEMBL:'+old_key] = ensembl2entrez.pop(old_key) # make sure name is infront
        concept_dct.update(ensembl2entrez)
    else:
        print('no ensembl IDs can be mapped to entrez IDs')

    # remove genes for which no entrez id was found (value is nan)
    concept_dct = {k: concept_dct[k] for k in concept_dct if not pd.isna(concept_dct[k])}
    return concept_dct


def hit_dgidb_api(node = 'genes=HTT', rows = 2000):
    """
    This function performs api calls to DGIdb to retrieve out and in edges from a query node.
    It retrieves entities plus associations via the DGIdb API service.
    The interactions endpoint can be used to return interactions for a given set of gene or drug names/symbols

    It returns edges.

    :param node: A comma delimited list of gene or drug names/symbols (string). Default: 'genes=HTT', for multiple genes, seperate with comma
    :param rows: the maximum number of results to return (integer). Default: 2000.
    :return: an api response objects.
    """
    
    print('\nThe function "hit_dgidb_api()" is running...This may take a while')
    # api address
    apilink = 'https://dgidb.org/api/v2/interactions.json?'
    
    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': rows}
    # edges: 
    edges = requests.get('{}{}'.format(apilink,node),params=parameters)
    
    print('\nFinished "hit_dgidb_api()".\n')

    return edges

def get_edges_objects(edges):
    """
    This function prepares the api object responses from DGIdb.
    It returns four lists, one for subjects, relations, objects, and references.
    Subjects, relations and objects are lists of dictionaries, where each dictionary is a node.
    References list lists strings, where each string is a chain of references for each edge.

    :param edges: DGIdb API response object
    :return: subjects, relations, objects and references lists (in this order)
    """

    # variables
    drug_dict = dict() #key=searchterm, value=list of drugs for the search term

    # compose list of dictionaries
    for associations in edges.json()['matchedTerms']: #associations

        searchterm = associations['searchTerm']
        entrez_id = str(associations['entrezId'])
        if not searchterm ==entrez_id:
            print('not the same:', searchterm, entrez_id)
            continue

        for interaction in associations['interactions']:
            if searchterm not in drug_dict:
                drug_dict[searchterm] = [interaction]
            else:
                drug_dict[searchterm].append(interaction)

    for associations in edges.json()['ambiguousTerms']: #associations

        searchterm = associations['searchTerm']
        entrez_id = str(associations['entrezId'])
        if not searchterm ==entrez_id:
            print('not the same:', searchterm, entrez_id)
            continue

        for interaction in associations['interactions']:
            if searchterm not in drug_dict:
                drug_dict[searchterm] = [interaction]
            else:
                drug_dict[searchterm].append(interaction)

    return drug_dict


def create_csv(id_name_dict, drug_dict, nodes_df, filename='DGIdb_network'):
    """
    This function saves the DGIdb network as a csv file
    :param id_name_dict: dict of gene ids/name to convert
    :param drug_dict: dict of searchterm/drugs 
    :param nodes_df: the dataframe of all monarch nodes 
    :param filename: the name of the csv file to save
    :return: location of the saved csv file
    """   
    
    df = pd.DataFrame(columns = ['target_id', 'target_name', 'target_iri', 'target_category', 'interaction_id', 'interaction_types', 'interaction_iri', 'drug_id', 
                                    'drug_name', 'drug_iri', 'drug_category', 'pm_ids', 'sources', 'score'] )
    
    ## create dict to turn activationtype into ols ontology
    interaction_dict = {}
    
    #inhibits
    interaction_dict['RO:0002408'] = ['antagonist', 'antibody', 'antisense oligonucleotide', 'blocker', 'cleavage',
                                      'inhibitor', 'inhibitory allosteric modulator', 'inverse agonist', 
                                      'negative modulator', 'partial antagonist', 'suppressor'] 
    #activates
    interaction_dict['RO:0002406'] = ['activator', 'agonist', 'chaperone', 'cofactor', 'inducer', 'partial agonist',
                                      'positive modulator', 'stimulator', 'vaccine'] 
    #regulates
    interaction_dict['RO:0011002'] = ['NA', 'None', 'n/a', 'other/unknown', 'adduct', 'allosteric modulator',
                                      'binder', 'ligand', 'modulator', 'multitarget', 'potentiator',
                                      'product of', 'substrate']      
    # turn into value:key
    interaction_dict_conv = {}
    for k,vs in interaction_dict.items():
        for v in vs:
            interaction_dict_conv[v] = k
            
    for searchterm, drug_l in drug_dict.items():
        interactionId_col = []
        interactionTypes_col = []
        interactionIri_col = []
        drugId_col = []
        drugName_col = []
        drugIri_col = []
        drugCategory_col = []
        pmids_col = []
        sources_col = []
        score_col = []
        
        for i in range(len(drug_l)):
            interactionId_col.append(drug_l[i]['interactionTypes']) # first type, later change to id
            interactionTypes_col.append(drug_l[i]['interactionTypes'])
            interactionIri_col.append(drug_l[i]['interactionTypes']) # first type, later change to iri
            drugId_col.append(drug_l[i]['drugConceptId'])
            drugName_col.append(drug_l[i]['drugName'])
            drugIri_col.append(drug_l[i]['drugConceptId']) # first id, later change to iri
            drugCategory_col.append('drug')
            pmids_col.append(drug_l[i]['pmids'])
            sources_col.append(drug_l[i]['sources'])
            score_col.append(drug_l[i]['score'])
            
        
        df_one_gene = pd.DataFrame(list(zip(interactionId_col, interactionTypes_col, interactionIri_col, drugId_col, 
                                   drugName_col, drugIri_col, drugCategory_col, pmids_col, sources_col, score_col)), 
                          columns =['interaction_id', 'interaction_types', 'interaction_iri', 'drug_id', 
                                    'drug_name', 'drug_iri', 'drug_category', 'pm_ids', 'sources', 'score'])
        
        #turn interaction id into ontology
        df_one_gene['interaction_id'] = df_one_gene['interaction_id'].apply(lambda x: x[0] if len(x)>0 else 'NA') # keep first 
        df_one_gene=df_one_gene.replace({"interaction_id": interaction_dict_conv})
        
        #create interaction iri column
        df_one_gene['interaction_iri'] = df_one_gene['interaction_id'].values #copy all ids from the just created id column
        df_one_gene['interaction_iri'] = ['http://purl.obolibrary.org/obo/' + idx.replace(':', '_') for idx in df_one_gene['interaction_iri']]

        #create drug iri column 
        iricol=[]
        for idx in df_one_gene['drug_iri']:
            iri = 'https://identifiers.org/'+idx
            iricol.append(iri)
        df_one_gene['drug_iri'] = iricol
        
        target_id_col = [searchterm]*len(df_one_gene.index) # create column of search term gene id, using the id_name_dict of genes
        target_name_col = [id_name_dict[searchterm]]*len(df_one_gene.index)
        target_iri_col = [nodes_df.loc[nodes_df['id'] == searchterm, 'uri'].iloc[0]]*len(df_one_gene.index) # create column of search term iri using nodes_df
        target_category_col = ['gene']*len(df_one_gene.index)
        df_one_gene.insert(loc=0, column='target_category', value=target_category_col)
        df_one_gene.insert(loc=0, column='target_iri', value=target_iri_col)       
        df_one_gene.insert(loc=0, column='target_name', value=target_name_col)
        df_one_gene.insert(loc=0, column='target_id', value=target_id_col)
        
        df_one_gene.interaction_types = df_one_gene.interaction_types.apply(lambda y: 'NA' if len(y)==0 else y) #replace empty lists in the column interaction_types with 'None'
        
        df = pd.concat([df, df_one_gene])

    # print output file
    path = os.getcwd() + '/DGIdb'
    if not os.path.isdir(path): os.makedirs(path)
    pd.DataFrame(df).fillna('NA').to_csv('{}/{}_v{}.csv'.format(path, filename, today), index=False)

    return print("\nSaving DGIdb network at: '{}/{}_v{}.csv'...\n".format(path, filename, today))


def read_connections(filename):
    """
    This function reads DBIdb CSV file.
    :param filename: the csv file string
    :return: dataframe
    """

    path = os.getcwd() + '/DGIdb'
    csv_path = path + '/' + filename

    network_df = pd.read_csv('{}'.format(csv_path))
    print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    print('* These are the attributes: {}'.format(network_df.columns))
    print('* This is the first record:\n{}'.format(network_df.head(1)))

    return network_df

def build_edges(edges_df):
    """
    This function builds the edges network with the graph schema.
    :param edges_df: network dataframe 
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_edges()" is running...')
    
    edges_df['sources'] = edges_df['sources'].apply(lambda x: ast.literal_eval(str(x))) #change string of list into lst
    ## variables
    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['target_id'] = tuple[0]
            record['target_name'] = tuple[1]
            record['target_iri'] = tuple[2]    
            record['target_category'] = tuple[3]  
            record['interaction_id'] = tuple[4]
            record['interaction_types'] = tuple[5]
            record['interaction_iri'] = tuple[6]
            record['drug_id'] = tuple[7]
            record['drug_name'] = tuple[8]
            record['drug_iri'] = tuple[9]
            record['drug_category'] = tuple[10]
            record['pm_ids'] = tuple[11]
            record['sources'] = tuple[12]
            record['score'] = tuple[13]
            connections_l.append(record)
        edges_df = pd.DataFrame(connections_l)

    # generate static variable: uriPrefixes_dct (url references)
    uriPrefixes_dct = {
        'pmid': 'https://pubmed.ncbi.nlm.nih.gov/' 
    }

    # provenance variables
    ref_text = 'This edge comes from the DGIdb 4.0.'  
    ref_date = today

    ## build graph schema network edges data structure and save edges file
    edges_l = list()
    for edge in edges_df.itertuples():
        
        pmid_l = edge.pm_ids
        pmid_l = pmid_l.strip('][').split(', ') #convert str representation of list to list
        pmid_l = [n.strip() for n in pmid_l]
        # create multi-term pubmed url
        if len(pmid_l[0]):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
        else:
            ref_uri = 'NA'

        # prepare edge attributes: sub_id, obj_id, rel_id, rel_label, rel_def, rel_iri
        sub_id = 'NA' if edge.target_id is None or str(edge.target_id) == 'nan' else edge.target_id
        rel_id = 'NA' if edge.interaction_id is None or str(edge.interaction_id) == 'nan' else edge.interaction_id
        obj_id = 'NA' if edge.drug_id is None or str(edge.drug_id) == 'nan' else edge.drug_id
        rel_label = 'NA' if edge.interaction_types is None or str(edge.interaction_types) == 'nan' else "|".join(eval(edge.interaction_types)) #turn list into string 
        rel_iri = 'NA' if edge.interaction_iri is None or str(edge.interaction_iri) == 'nan' else edge.interaction_iri
        rel_def = 'interaction score = '+str(edge.score) #property description: add interaction scores here
            
        # build the data structure = list of edges as list of dict, where a dict is an edge
        edge = dict()
        edge['subject_id'] = obj_id # drug
        edge['object_id'] = sub_id # gene
        edge['property_id'] = rel_id
        edge['property_label'] = rel_label
        edge['property_description'] = rel_def
        edge['property_uri'] =  rel_iri
        edge['reference_uri'] = ref_uri
        edge['reference_supporting_text'] = ref_text
        edge['reference_date'] = ref_date
        edges_l.append(edge)

    # save edges file
    df = pd.DataFrame(edges_l)
    print('df',df.shape)
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', \
             'property_label', 'property_description', 'property_uri']]
    path = os.getcwd() + '/DGIdb'
    df.fillna('NA').to_csv('{}/DGIdb_edges_v{}.csv'.format(path,today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe Monarch network edges are built and saved at: {}/DGIdb_edges_v{}.csv\n'.format(path,today))
    print('\nFinished build_edges().\n')

    return df

def build_nodes(edges_df):
    """
    This function builds the nodes network with the graph schema.
    :param edges_df: network dataframe 
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_nodes()" is running...')
    
    edges_df['sources'] = edges_df['sources'].apply(lambda x: ast.literal_eval(str(x))) #change string of list into lst
    
    # build semantic groups dictionary
    # collide concepts in a concept dict
    concept_dct = dict()

    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['target_id'] = tuple[0]
            record['target_name'] = tuple[1]
            record['target_iri'] = tuple[2]    
            record['target_category'] = tuple[3]  
            record['interaction_id'] = tuple[4]
            record['interaction_types'] = tuple[5]
            record['interaction_iri'] = tuple[6]
            record['drug_id'] = tuple[7]
            record['drug_name'] = tuple[8]
            record['drug_iri'] = tuple[9]
            record['drug_category'] = tuple[10]
            record['pm_ids'] = tuple[11]
            record['sources'] = tuple[12]
            record['score'] = tuple[13]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)

    for edge in edges_df.itertuples():
        sid = edge.target_id 
        oid = edge.drug_id 
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    print('Number of concepts: {}'.format(len(concept_dct.keys())))

    # build concept attributes dict: id integration of sub and obj IDs in a common data structure
    concept_dct = dict()

    for edge in edges_df.itertuples():
        # id: integration of sub and obj IDs in a unique data structure
        sid = edge.target_id 
        slab = edge.target_name 
        sgroup = edge.target_category
        siri = edge.target_iri
        oid = edge.drug_id 
        olab = edge.drug_name
        ogroup = edge.drug_category
        oiri = edge.drug_iri

        concept_dct[sid] = {'preflabel': slab,
                            'semantic_groups': sgroup,
                            'uri': siri,
                            'name': 'NA',
                            'synonyms': 'NA', 
                            'description': 'NA'}
        concept_dct[oid] = {'preflabel': olab,
                            'semantic_groups': ogroup,
                            'uri': oiri,
                            'name': 'NA',
                            'synonyms': 'NA', 
                            'description': 'NA'}

    # prepare data structure = [ {} ... {} ], where every {} = concept = row
    nodes_l = list()
    for concept in concept_dct:
        # define nodes (rows) for the data structure
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = concept_dct[concept]['semantic_groups']
        node['uri'] = concept_dct[concept]['uri']
        node['preflabel'] = concept_dct[concept]['preflabel']
        node['name'] = concept_dct[concept]['preflabel'] 
        node['synonyms'] = 'NA'
        node['description'] = 'NA' 
        nodes_l.append(node)

    # save nodes file
    df = pd.DataFrame(nodes_l)
    df = df[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    path = os.getcwd() + '/DGIdb'
    df.fillna('NA').to_csv('{}/DGIdb_nodes_v{}.csv'.format(path,today), index=False)

    # print info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe Monarch network nodes are built and saved at: {}/DGIdb_nodes_v{}.csv\n'.format(path,today))
    print('\nFinished build_nodes().\n')

    return df

def run_dgidb(date):
    """
    This function runs the whole DGIdb script and saves nodes and edges files.
    :param date: the date of creation of the disease graph
    :return: nodes and edges files in /DGIdb folder
    """
    
    path = os.getcwd() + '/DGIdb'
    if not os.path.isdir(path): os.makedirs(path)
    
    ## build DGIdb network
    monarch_nodes_dis_file = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    monarch_nodes_symp_file = './monarch/monarch_nodes_symptom_v{}.csv'.format(today)
    # open csv files
    nodes_dis_df = pd.read_csv(monarch_nodes_dis_file)
    nodes_symp_df = pd.read_csv(monarch_nodes_symp_file)
    # concatenate dfs
    combined_monarch_nodes_df = pd.concat([nodes_dis_df, nodes_symp_df])
    
    ## find seedlist (all genes)
    id_name_dict = get_genes(combined_monarch_nodes_df)

    ## map HGNC and other ontologies to ENTREZ such that it can be used as input for DGIdb
    id_to_entrez_dct = normalize_genes_to_graph(id_name_dict)
    entrez_to_id_dict = {y: x for x, y in id_to_entrez_dct.items()} # turn into entrez to id dict
    seed_id_lst = ','.join(list(id_to_entrez_dct.values()))
    
    ## get drugs that interact with these genes
    edges = hit_dgidb_api(node = 'genes='+seed_id_lst, rows = 2000)
    gene_drug_dict = get_edges_objects(edges)
    ## turn entrez id keys into original ids
    for old_key in list(gene_drug_dict):
        gene_drug_dict[entrez_to_id_dict[old_key]] = gene_drug_dict.pop(old_key)

    ## create network
    create_csv(id_name_dict, gene_drug_dict, nodes_df=combined_monarch_nodes_df, filename='DGIdb_network')   
    
    ## read network
    DGIdb_connections=read_connections('DGIdb_network_v{}.csv'.format(today))

    ## build nodes and edges files
    build_edges(DGIdb_connections)
    build_nodes(DGIdb_connections)


if __name__ == '__main__':
    run_dgidb()



