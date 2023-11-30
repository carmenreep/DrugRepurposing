# -*- coding: utf-8 -*-
"""
Module for Monarch network preparation and management

Created on Sat May 7 11:51:22 2022

@author: Carmen Reep
"""

import requests
import sys,os
import datetime
import pandas as pd
from biothings_client import get_client
from tqdm import tqdm


# timestamp
today = datetime.date.today()
curr_year = int(str(today)[:4])


def read_connections(filename):
    """
    This function reads monarch_connections CSV file.
    :param filename: complete path to the monarch connections csv file string
    :return: monarch edges dataframe
    """

    # monarch network
    path = os.getcwd() + '/monarch'
    csv_path = path + '/' + filename
    network_df = pd.read_csv('{}'.format(csv_path))
    print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    print('* These are the attributes: {}'.format(network_df.columns))
    print('* This is the first record:\n{}'.format(network_df.head(1)))

    return network_df


# retrieve subnetwork from Monarch knowledge graph

def hit_monarch_api(node = 'HGNC:4851', rows = 2000):
    """
    This function performs api calls to Monarch to retrieve out and in edges from a query node.
    It retrieves entities plus associations via the BioLink API service.
    It hits two endpoints:
        * association/from - for out edges
        * association/to - for in edges
    It returns out and in edges.

    :param node: node id to query (string). Default: 'HGNC:4851' (=HTT gene).
    :param rows: the maximum number of results to return (integer). Default: 2000.
    :return: two api response objects: 'out' and 'in' response objects, in this order.
    """
    
    # api address
    biolink = 'https://api.monarchinitiative.org/api/association'
    
    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': rows}
    # out edges: from/
    r_out = requests.get('{}/from/{}'.format(biolink,node),params=parameters)

    # in edges: to/
    r_in = requests.get('{}/to/{}'.format(biolink,node),params=parameters)

    return r_out, r_in

def get_disease_name_id(disease_input_ID = 'OMIM:143100'):
    """
    This function finds the disease name and disease ID of the input OMIM ID.
    It returns two strings, name and ID respectively.

    :param disease_ID: The input ID of the disease
    :return: two strings, name and ID respectively
    """
    # api address
    biolink = 'https://api.monarchinitiative.org/api/association'
    
    # parameters
    parameters = {'fl_excludes_evidence': False, 'rows': 1}
    # out edges: from/
    r_out_disease = requests.get('{}/from/{}'.format(biolink,disease_input_ID),params=parameters)

    for association in r_out_disease.json()['associations']:
        disease_name = association['subject']['label']
        disease_id = association['subject']['id']
    
        
    return disease_name, disease_id


def get_edges_objects(r_out, r_in):
    """
    This function prepares the api object responses from Monarch.
    It returns four lists, one for subjects, relations, objects, and references.
    Subjects, relations and objects are lists of dictionaries, where each dictionary is a node.
    References list lists strings, where each string is a chain of references for each edge.

    :param r_out: BioLink API 'out' response object
    :param r_in: BioLink API 'in' response object
    :return: subjects, relations, objects and references lists (in this order)
    """

    # variables
    sub_l = list()
    rel_l = list()
    obj_l = list()
    ref_l = list()

    # compose list of dictionaries
    for associations in [r_out.json()['associations'], r_in.json()['associations']]:
        for association in associations:
            pub_l = list()
            sub_l.append(association['subject'])
            rel_l.append(association['relation'])
            obj_l.append(association['object'])
            # add references to each association as a list of strings
            if association['publications']:
                for publication in association['publications']:
                    pub_l.append(publication['id'])
            else:
                pub_l.append('NA')
            ref_l.append('|'.join(pub_l))

    return sub_l, rel_l, obj_l, ref_l


def get_edges(sub_l, rel_l, obj_l, ref_l, attribute='id'):
    """
    This function builds edges using a user-specified attribute for each node.
    It returns a set of edges, where edges are tuples.

    :param sub_l: subjects (objects) list from the get_edges_objects() function
    :param rel_l: relations (objects) list from the get_edges_objects() function
    :param obj_l: objects (objects) list from the get_edges_objects() function
    :param ref_l: references (strings) list from the get_edges_objects() function
    :param attribute: object attribute, default 'id'
    :return: edges (as tuples) set
    """
    edges = set()
    # compose tuple
    for i in range(len(sub_l)):
        sub = sub_l[i][attribute]
        rel = rel_l[i][attribute]
        obj = obj_l[i][attribute]
        ref = ref_l[i]
        edges.add((sub, rel, obj, ref))

    return edges


def keep_edges(keep, new):
    """
    This function adds edges from a new set to a keep set.
    :param keep: edges set
    :param new: edges set
    :return: updated edges set
    """

    for edge in new:
        keep.add(edge)

    return keep


def keep_nodes(keep, edges, seed):
    """
    This function collects nodes from the edges that are not in the nodes query list to keep nodes set.
    It filters out: PMID nodes, and nodes related by provenance:
        * None
        * 'dc:source'
        * 'IAO:0000136' or 'is about'
        * 'IAO:0000142' or 'mentions'
    i.e., not biologically related
    It returns a set of nodes.

    :param keep: nodes set
    :param edges: edges set
    :param seed: query nodes list
    :return: updated nodes set
    """

    for (sub, rel, obj, ref) in edges:
        if 'PMID' in sub or 'PMID' in obj:
            continue
        if rel == None:
            rel = 'None'
        if 'dc:source' in rel:
            continue
        if 'IAO:0000136' in rel:  # is about
            continue
        if 'IAO:0000142' in rel:  # mentions
            continue
        if sub not in seed:
            keep.add(sub)
        if obj not in seed:
            keep.add(obj)

    return keep


def get_neighbours(seed):
    """
    This function gets the first layer of neighbours and relations.
    :param seed: query nodes list
    :return: nodes set, edges set (in this order)
    """

    keepNodes = set()
    keepEdges = set()
    seedNodes = set(seed)
    for node in tqdm(seedNodes):
        try:
            r_out, r_in = hit_monarch_api(node)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            keepEdges = keep_edges(keepEdges, edges)
            keepNodes = keep_nodes(keepNodes, edges, seedNodes)

        except (ValueError, KeyError):
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            print(node)

    return keepNodes, keepEdges


def filter_edges(nodes, edges):
    """
    This function filters down edges with both nodes in a nodes list.

    :param nodes: nodes list
    :param edges: edges set
    :return: filtered edges set
    """
    nodes = set(nodes)
    keep = set()
    for (start, pred, stop, ref) in edges:
        if {start, stop} <= nodes: # is {..} a subset of {..}
            keep.add((start, pred, stop, ref))

    return keep


def add_attributes(sub_l, rel_l, obj_l, edges):
    """
    This function adds 'label', 'iri', 'category' attribute to each entity in the edge.
    :param sub_l: subjects (object) list
    :param rel_l: relations (object) list
    :param obj_l: objects (object) list
    :param edges: edges set
    :return: metaedges set
    """

    metaedges = set()
    for (sub_id, rel_id, obj_id, refs) in edges:
        for i in range(len(sub_l)):
            if sub_l[i]['id'] == sub_id and rel_l[i]['id'] == rel_id and obj_l[i]['id'] == obj_id:
                metaedges.add((sub_l[i]['id'],
                               sub_l[i]['label'],
                               sub_l[i]['iri'],
                               sub_l[i]['category'][0], # add [0] because category is in a list
                               rel_l[i]['id'],
                               rel_l[i]['label'],
                               rel_l[i]['iri'],
                               obj_l[i]['id'],
                               obj_l[i]['label'],
                               obj_l[i]['iri'],
                               obj_l[i]['category'][0],
                               refs)
                              )
                break
    return metaedges


def keep_node_type(edges, seed, nodeType='ortho'):
    """
    This function keeps specific node types for objects in edges.

    :param edges: edges set
    :param seed: the query nodes list
    :param nodeType: Introduce node type to keep (string): 'ortho' for orthologs or 'pheno' \
    for phenotypes/diseases, default is 'ortho'
    :return: nodes set
    """

    propertyList = ['RO:HOM0000017', 'RO:HOM0000020']
    if nodeType == 'pheno':
        propertyList = ['RO:0002200', 'RO:0002607', 'RO:0002326', 'GENO:0000840']

    keep = set()
    for (sub, rel, obj, ref) in edges:
        if rel == None:
            continue
        if rel in propertyList:
            if sub not in seed:
                keep.add(sub)
            if obj not in seed:
                keep.add(obj)

    return keep


def get_connections(nodes):
    """
    This function returns associations retrieved from Monarch among a list of query nodes.
    :param nodes: the query nodes list
    :return: edges set
    """""

    keep = set()
    for node in tqdm(nodes):
        try:
            r_out, r_in = hit_monarch_api(node, 1000)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            filteredEdges = filter_edges(nodes, edges)
            metaFilteredEdges = add_attributes(sub_l, rel_l, obj_l, filteredEdges)
            keep = keep_edges(keep, metaFilteredEdges)

        except (ValueError, KeyError):
            pass
        except:
            print('error: {}'.format(sys.exc_info()[0]))
            print(node)

    return keep


# NETWORK MANAGEMENT FUNCTIONS

def get_neighbours_list(seed_list):
    """
    This function returns the first explicit layer of neighbours from a list of query nodes.
    :param seed_list: biomedical entities list, where each entity is the identifier string like 'HGNC:4851'
    :return: neighbours list
    """

    # print executing function
    print('\nThe function "get_neighbours_list()" is running. Its runtime may take some minutes. '
          'If you interrupt the process, you will lose all the nodes retrieved '
          'and you should start over the execution of this function.')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)
    print('\nFinished get_neighbours_list().\n')

    return list(neighbours)


def get_orthopheno_list(seed_list):
    """
    This function returns orthologs-phenotypes nodes in ortho-pheno relationships for a list of query genes.
    :param seed_list: gene list, where each gene is the identifier string like 'HGNC:4851'
    :return: orthopheno list
    """

    # print executing function
    print('\nThe function "get_orthopheno_list()" is running. Its runtime may take some hours. '
          'If you interrupt the process, you will lose all the nodes retrieved '
          'and you should start over the execution of this function.')

    # get first layer of neighbour nodes
    neighbours, relations = get_neighbours(seed_list)

    # keep orthologs in the first layer
    orthologs = keep_node_type(relations, seed_list)

    # get second layer from orthologs
    neighbours, relations = get_neighbours(orthologs)

    # keep phenotypes in the second layer
    phenotypes = keep_node_type(relations, orthologs, 'pheno')

    nodes = set()
    nodes.update(orthologs, phenotypes)
    print('\nFinished get_orthopheno_list().\n')

    return list(nodes)


def extract_edges(gene_list):
    """
    This function returns the Monarch network from a list of query nodes. It retrieves connectivity from Monarch, i.e. \
    edges from Monarch between query nodes.
    :param gene_list: gene list
    :return: edges (as tuples) set
    """

    # print executing function
    print('\nThe function "extract_edges()" is running. Its runtime may take some hours. '
          'If you interrupt the process, you will lose all the edges retrieved '
          'and you should start over the execution of this function.')

    # set network nodes: gene list provided by the user
    nodes = set(gene_list)

    # get connections
    network = get_connections(nodes)
    print('\nFinished extract_edges(). To save the retrieved Monarch edges use the function "print_network()".\n')

    return network


def _print_network2(network, filename):
    """This function saves the Monarch expanded network into a CSV file. this function save connections file format into
    get-monarch-connections/"""

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
        f.write(
            'subject_id,subject_label,subject_uri,subject_category,relation_id,relation_label,relation_uri,object_id,object_label,object_uri,object_category,reference_id_list\n')
        for edge in network:
            edge = ['None' if t is None else '"{}"'.format(str(t)) for t in edge]
            f.write('{}\n'.format(','.join(edge)))

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))


def print_network(network, filename):
    """
    This function saves the Monarch network into a CSV file. It only prints connections file format only.
    :param network: monarch edges dataframe
    :param filename: file name without extension string, e.g. 'monarch_connections'
    :return: None object
    """

    # transform set of tuples to list of dictionaries
    edges = list()
    for tuple in network:
        row = dict()
        row['subject_id'] = tuple[0]
        row['subject_label'] = tuple[1]
        row['subject_uri'] = tuple[2]
        row['subject_category'] = tuple[3]
        row['relation_id'] = tuple[4]
        row['relation_label'] = tuple[5]
        row['relation_uri'] = tuple[6]
        row['object_id'] = tuple[7]
        row['object_label'] = tuple[8]
        row['object_uri'] = tuple[9]
        row['object_category'] = tuple[10]
        row['reference_id_list'] = tuple[11]
        edges.append(row)
    
    df=pd.DataFrame(edges).fillna('None')
    
    # change every 'model' node that does not have 'Coriell' as prefix to genotype
    df.loc[(df['subject_category'] == 'model') &  (df['subject_id'].str.contains('zfin|mgi|flybase|wormbase', case=False)), 'subject_category'] = 'genotype'
    df.loc[(df['object_category'] == 'model') &  (df['object_id'].str.contains('zfin|mgi|flybase|wormbase', case=False)), 'object_category'] = 'genotype'

    # remove every row that has Coriell in either subject or object id
    df.drop(df[df['object_id'].str.contains('coriell', case=False)].index, inplace=True)
    df.drop(df[df['subject_id'].str.contains('coriell', case=False)].index, inplace=True)

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    df.to_csv('{}/{}_v{}.csv'.format(path, filename, today), index=False)

    return print("\nSaving Monarch edges at: '{}/{}_v{}.csv'...\n".format(path, filename, today))


def print_nodes(nodes, filename):
    """
    This function saves Monarch nodes into a CSV file.
    :param nodes: nodes list
    :param filename: file name without path and extension, e.g. 'monarch_nodes'
    :return: None object
    """

    # print output file
    path = os.getcwd() + '/monarch'
    if not os.path.isdir(path): os.makedirs(path)
    with open('{}/{}_v{}.csv'.format(path, filename, today), 'w') as f:
        f.write('{}\n'.format('\n'.join(list(nodes))))

    return print("\nFile '{}/{}_v{}.csv' saved.".format(path, filename, today))


def expand_edges(seed_list):
    """
    This function returns the Monarch network expanded with the first layer of neighbors from a list of query nodes.
    This function builds monarch 1shell network.
    This function receives a list of nodes and returns a network or edges from Monarch.

    :param seed_list: the query nodes list
    :return: edges set
    """

    # get 1shell list of nodes or neighbors
    neighbours, relations = get_neighbours(seed_list)

    # network nodes:  seed + 1shell
    nodes = set(seed_list).union(neighbours)

    # get connections for network nodes
    network = get_connections(nodes)

    return network


def orthopheno_expand_edges(seed_list):
    """
    This function returns the Monarch network expanded with the orthologs and the ortholog associated phenotypes from
     a list of query nodes.
    This function builds monarch 1shell-animal network.
    This function receives a list of nodes and returns a network or edges from Monarch.

    :param seed_list: the query nodes list
    :return: edges set
    """

    # get ortholog-phenotypes list
    orthophenoList = get_orthopheno_list(seed_list)

    # network nodes:  seed + neighbors + orthologs + phenotypes
    nodes = set(seed_list).union(set(orthophenoList))

    # get connections for network nodes
    network = get_connections(nodes)

    return network


# BUILD NETWORK

def build_edges(edges_df, edges_fname):
    """
    This function builds the edges network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_edges()" is running...')
    ## variables
    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            record['subject_uri'] = tuple[2]
            record['subject_category'] = tuple[3]
            record['relation_id'] = tuple[4]
            record['relation_label'] = tuple[5]
            record['relation_uri'] = tuple[6]
            record['object_id'] = tuple[7]
            record['object_label'] = tuple[8]
            record['object_uri'] = tuple[9]
            record['object_category'] = tuple[10]
            record['reference_id_list'] = tuple[11]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)


    # generate static variable: uriPrefixes_dct (url references)
    uriPrefixes_dct = {
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',  
        'react': 'http://reactome.org/content/detail/',  
        'zfin': 'http://zfin.org/',
        'go_ref': 'http://purl.obolibrary.org/obo/go/references/',  
        'mgi': 'http://www.informatics.jax.org/accession/MGI:',  
        'flybase': 'http://flybase.org/reports/',
        'wormbase': 'http://www.wormbase.org/resources/paper/',
        'hpo': 'http://compbio.charite.de/hpoweb/showterm?id=HP:',
        'isbn-10': 'ISBN-10:',
        'isbn-13': 'ISBN-13:',
        'mondo': 'http://purl.obolibrary.org/obo/MONDO_',  
        'rgd': 'https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=', 
        'omim': 'http://omim.org/entry/',  
        'sgd_ref': 'https://db.yeastgenome.org/reference/',  
        'genereviews': 'https://www.ncbi.nlm.nih.gov/books/',  
        'omia': 'http://omia.angis.org.au/',  
        'hgnc': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:', 
    }

    # generate static variable: dbPrefixes_dct (source/database references)
    dbPrefixes_dct = {
        'na': 'NA',
        'nan': 'NA',
        'mgi': 'http://www.informatics.jax.org/',
        'fb': 'http://flybase.org/',
        'rgd': 'http://rgd.mcw.edu/',
        'zfin': 'http://zfin.org/',
        'sgd': 'https://www.yeastgenome.org/',
        'hgnc': 'https://www.genenames.org/',
        'xenbase': 'http://www.xenbase.org/'
    }

    # provenance variables
    ref_text = 'This edge comes from the Monarch Initiative {}.'.format(curr_year) 
    ref_date = today

    ## build graph schema network edges data structure and save edges file
    # prepare dataframe = [ {} ... {} ], where every row = {} = concept
    edges_l = list()
    for edge in edges_df.itertuples():
        # edge or row is a tuple (named and ordered attributes)
        # edge.reference_id_list >> can be 1) np.nan (float type) or 2) str without "|" 3) str with "|"
        ref_s = str(edge.reference_id_list)

        # prepare reference_uri_list attribute
        ref_uri_l = list()
        # expand to uri or NA
        pmid_l = list()
        # reference_id list iteration
        for ref in ref_s.strip().split('|'):
            # NA or database
            if ':' not in ref:
                try:
                    ref_uri = dbPrefixes_dct[ref.lower()]
                except KeyError:
                    print("Warning:")
                    print('Detected a new reference database in Monarch not yet implemented in this module. '
                          'The new source should be added to the dictionary of databases.'
                          'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                    print("In the build_edges() method, update 'dbPrefixes_dct' dictionary with '{}'".format(ref))
                    print('The edge that includes this new reference database is {}'.format(edge))
                    print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                          'until the dictionary is updated.')
                ref_uri_l.append(ref_uri)
            # publication uri: pubmed_id or url
            else:
                pref, uriId = ref.split(':')
                # separate pmid from non pmid and detect:
                # pubmed_id
                if ref.startswith('PMID'):
                    pmid_l.append(uriId)
                # url
                elif ref.lower().startswith('http'):
                    ref_uri_l.append(ref)
                else:
                    try:
                        ref_uri = uriPrefixes_dct[pref.lower()] + uriId
                    except KeyError:
                        print("Warning:")
                        print('Detected a new reference source in Monarch not yet implemented in this module. '
                              'The new source should be added to the dictionary of sources.'
                              'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                        print("In the build_edges() method, update 'uriPrefixes_dct' dictionary with '{}'".format(pref))
                        print('The edge that includes this new reference source is {}'.format(edge))
                        print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                              'until the dictionary is updated.')
                    ref_uri_l.append(ref_uri)
        # create multi-term pubmed url
        if len(pmid_l):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
            ref_uri_l.append(ref_uri)
        ref_uri_list = '|'.join(ref_uri_l)

        # prepare edge attributes: sub_id, obj_id, rel_id, rel_label, rel_def, rel_uri
        sub_id = 'NA' if edge.subject_id is None or str(edge.subject_id) == 'nan' or str(edge.subject_id) == 'None' else edge.subject_id
        rel_id = 'NA' if edge.relation_id is None or str(edge.relation_id) == 'nan'  or str(edge.relation_id) == 'None' else edge.relation_id
        obj_id = 'NA' if edge.object_id is None or str(edge.object_id) == 'nan'  or str(edge.object_id) == 'None' else edge.object_id
        rel_label = 'NA' if edge.relation_label is None or str(edge.relation_label) == 'nan'  or str(edge.relation_label) == 'None' else edge.relation_label
        rel_uri = 'NA' if edge.relation_uri is None or str(edge.relation_uri) == 'nan'  or str(edge.relation_uri) == 'None' else edge.relation_uri
        rel_def = 'NA'
        
        # make sure there are no 'NA' relation ids
        if rel_id == 'NA':
            if edge.object_category == 'genotype': # if object is genotype, relation is 'has genotype'
                rel_id = 'GENO:0000222'
                rel_label = 'has_genotype'
                rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000222'
            else:
                rel_id = 'RO:0002610' # else (gene, model objects), relation is 'correlated with'
                rel_label = 'correlated with'
                rel_uri = 'http://purl.obolibrary.org/obo/RO_0002610'
        
        # change every 'contributes to' to 'contributes to condition'
        if rel_id == 'RO:0002326':
            rel_id = 'RO:0003304'
            rel_label = 'contributes to condition'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003304'
        
        # change every 'is causal germline mutation in' and 'is causal germline mutation partially giving rise to' and 'pathogenic for condition' to 'causes condition'
        if rel_id == 'RO:0004013' or rel_id == 'RO:0004016' or rel_id == 'GENO:0000840':
            rel_id = 'RO:0003303'
            rel_label = 'causes condition'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003303'

        # change every 'in orthology relationship with' to 'in 1 to 1 orthology relationship with'
        if rel_id == 'RO:HOM0000017':
            rel_id = 'RO:HOM0000020'
            rel_label = 'in 1 to 1 orthology relationship with'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_HOM0000020'  
        
        # if 'genotype' -->'has phenotype' --> 'disease', change 'has phenotype' to 'has role in modelling'
        if edge.subject_category == 'genotype' and edge.object_category == 'disease' and rel_id == 'RO:0002200':
            rel_id = 'RO:0003301'
            rel_label = 'has role in modeling'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003301'
        
        # change every 'is reference allele of' to 'is allele of'
        if rel_id == 'GENO:0000610':
            rel_id = 'GENO:0000408'
            rel_label = 'is_allele_of'
            rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000408'
        
        # build the data structure = list of edges as list of dict, where a dict is an edge
        edge = dict()
        edge['subject_id'] = sub_id
        edge['object_id'] = obj_id
        edge['property_id'] = rel_id
        edge['property_label'] = rel_label
        edge['property_description'] = rel_def
        edge['property_uri'] = rel_uri
        edge['reference_uri'] = ref_uri_list
        edge['reference_supporting_text'] = ref_text
        edge['reference_date'] = ref_date
        edges_l.append(edge)

    # save edges file
    df = pd.DataFrame(edges_l)
    print('df',df.shape)
    path = os.getcwd() + '/monarch'
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', \
             'property_label', 'property_description', 'property_uri']]
    df.fillna('NA').to_csv('{}/{}_v{}.csv'.format(path,edges_fname,today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe Monarch network edges are built and saved at: {}/monarch_edges_v{}.csv\n'.format(path,today))
    print('\nFinished build_edges().\n')

    return df


def build_nodes(edges_df, nodes_fname):
    """
    This function builds the nodes network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_nodes()" is running...')
    # build semantic groups dictionary
    # collide concepts in a concept dict
    concept_dct = dict()

    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            record['subject_uri'] = tuple[2]
            record['subject_category'] = tuple[3]
            record['relation_id'] = tuple[4]
            record['relation_label'] = tuple[5]
            record['relation_uri'] = tuple[6]
            record['object_id'] = tuple[7]
            record['object_label'] = tuple[8]
            record['object_uri'] = tuple[9]
            record['object_category'] = tuple[10]
            record['reference_id_list'] = tuple[11]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)

    for edge in edges_df.itertuples():
        sid = edge.subject_id 
        oid = edge.object_id 
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    print('Number of concepts: {}'.format(len(concept_dct.keys())))


    # build concept attributes dict: id integration of sub and obj IDs in a common data structure
    concept_dct = dict()

    # read edges from variable
    for edge in edges_df.itertuples():
        #fields = list(edge_tuple)
        # id: integration of sub and obj IDs in a unique data structure
        sid = edge.subject_id 
        slab = edge.subject_label 
        sgroup = edge.subject_category
        suri = edge.subject_uri
        oid = edge.object_id 
        olab = edge.object_label 
        ogroup = edge.object_category
        ouri = edge.object_uri
        
        # if 'genotype' -->'is allele of' --> 'gene', change 'genotype' to 'variant'
        if sgroup == 'genotype' and ogroup == 'gene' and edge.relation_id == 'GENO:0000408':
            sgroup = 'variant'
        
        # build the concept data structure
        concept_dct[sid] = {'preflabel': slab,
                            #'name': 'NA',
                            'semantic_groups': sgroup,
                            'uri': suri,
                            'synonyms': 'NA', 
                            'description': 'NA'}
        concept_dct[oid] = {'preflabel': olab,
                            #'name': 'NA',
                            'semantic_groups': ogroup,
                            'uri': ouri,
                            'synonyms': 'NA', 
                            'description': 'NA'}

    # biothings: annotate name,synonyms,description to genes
    print('\nAdding BioThings annotation: name, synonyms, description...')
    # input: (preflabel) symbol,alias
    symbols = list()
    for concept in concept_dct:
        if isinstance(concept_dct[concept]['semantic_groups'], list):
            for label in concept_dct[concept]['semantic_groups']:
                if 'gene' in label:
                    symbols.append(concept_dct[concept]['preflabel'])
        else:
            if 'gene' in concept_dct[concept]['semantic_groups']:
                symbols.append(concept_dct[concept]['preflabel'])
    
    print('symbols:', len(symbols))

    # query biothings
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='name,alias,summary', size=1, as_dataframe=True)

    # dictionary: {symbol:name}
    ids = (df.reset_index().rename(columns={'query': 'symbol'}))
    ids['synonyms'] = ids.alias.apply(lambda x: x if str(x) != 'nan' else 'NA')
    ids['description'] = ids.summary.apply(lambda x: x if str(x) != 'nan' else 'NA')
    monarch_s2n = dict(zip(ids.symbol, ids.name))
    monarch_s2s = dict(zip(ids.symbol, ids.synonyms))
    monarch_s2d = dict(zip(ids.symbol, ids.description))

    # prepare data structure = [ {} ... {} ], where every {} = concept = row
    nodes_l = list()
    for concept in concept_dct:
        # define nodes (rows) for the data structure
        preflabel = concept_dct[concept]['preflabel']
        concept_dct[concept]['synonyms'] = monarch_s2s[preflabel] if preflabel in monarch_s2s.keys() else 'NA'
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = concept_dct[concept]['semantic_groups']
        node['uri'] = concept_dct[concept]['uri']
        node['preflabel'] = preflabel
        node['name'] = monarch_s2n[preflabel] if preflabel in monarch_s2n.keys() else preflabel
        node['synonyms'] = '|'.join(list(concept_dct[concept]['synonyms'])) if isinstance(
            concept_dct[concept]['synonyms'], list) else concept_dct[concept]['synonyms']
        node['description'] = monarch_s2d[preflabel] if preflabel in monarch_s2d.keys() else 'NA'
        nodes_l.append(node)

    # save nodes file
    df = pd.DataFrame(nodes_l)
    df = df[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    path = os.getcwd() + '/monarch'
    df.fillna('NA').to_csv('{}/{}_v{}.csv'.format(path,nodes_fname, today), index=False)

    # print info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe Monarch network nodes are built and saved at: {}/monarch_nodes_v{}.csv\n'.format(path,today))
    print('\nFinished build_nodes().\n')

    return df

def get_symptoms_disease(disease_id, monarch_edges_csv, monarch_nodes_csv):
    """
    This function finds all symptoms of a disease.
    :param disease_id: id of the disease of which symptoms and drugs should be found (e.g. 'MONDO:0007739' for HD)  
    :param monarch_edges_csv: name of the edges csv file from Monarch
    :param monarch_nodes_csv: name of the nodes csv file from Monarch
    :return: a list of symptom names
    """
    # open csv file
    fname_edges = monarch_edges_csv
    fname_nodes = monarch_nodes_csv
    edges_df = pd.read_csv(fname_edges)
    nodes_df = pd.read_csv(fname_nodes)

    # find all nodes that are one step away from disease and have as relation 'has phenotype' (='RO:0002200')
    df_symptoms = edges_df.loc[(edges_df['subject_id'] == disease_id) & (edges_df['property_id'] == 'RO:0002200')]
    
    symptoms_id_lst = df_symptoms['object_id'].tolist()
    #get names of symptoms
    symptoms_name_lst = nodes_df.loc[nodes_df['id'].isin(symptoms_id_lst), 'preflabel'].to_list()
    return symptoms_name_lst

def symptom_list_specified_folder(input_folder = 'Huntington disease (2022-06-23)'):
    """
    This function finds the symptom list of the disease specified by user 
    clicking a folder and changes the cwd to the folder specified.
    :param input_folder: The input folder of the disease
    :return: symptom name list
    """
    date_str = input_folder[-11:-1] # get the date from the disease name
    date = datetime.datetime.strptime(date_str, '%Y-%m-%d').date()
    #change directory to the folder
    os.chdir(os.getcwd() +'/drugapp/data/' + input_folder)

    # find the disease ID that was previously written in a txt file
    with open(os.getcwd()+"/disease_id.txt", 'r') as file:
        file_text = file.read()
        input_number, disease_id, disease_name = file_text.split(';')

    monarch_fname_edges = './monarch/monarch_edges_disease_v{}.csv'.format(date)
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    
    # get list of symptoms for disease to ask user
    symptoms_name_lst=get_symptoms_disease(disease_id, monarch_fname_edges, monarch_fname_nodes)
    return symptoms_name_lst, date, input_number

def symptom_list_today():
    """
    This function finds the symptom list of the disease specified by user on 
    that same day.
    :return: symptom name list, date of creation of files required, disease name
    """
    date=today
    
    # find the disease ID that was previously written in a txt file
    with open(os.getcwd()+"/disease_id.txt", 'r') as file:
        file_text = file.read()
        input_number, disease_id, disease_name = file_text.split(';')

    disease_name_date = disease_name+' ({})'.format(date)
    
    monarch_fname_edges = './monarch/monarch_edges_disease_v{}.csv'.format(date)
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    
    # get list of symptoms for disease to ask user
    symptoms_name_lst=get_symptoms_disease(disease_id, monarch_fname_edges, monarch_fname_nodes)
    
    return symptoms_name_lst, date, disease_name_date

def run_monarch(input_number = '143100'):
    """
    This function runs the whole Monarch script and saves nodes and edges files.

    :param input_number: The input phenotype MIM number of the disease
    :return: nodes and edges files in /monarch folder
    """
    
    # turn input number into input ID
    input_id = 'OMIM:'+input_number
        
    # get the disease name (str)
    disease_name, disease_id = get_disease_name_id(disease_input_ID = input_id)
    
    # make folder for HD and change cwd to this folder
    new_path = os.getcwd() + '/drugapp/data/'+disease_name+' ({})'.format(today)
    os.makedirs(new_path)
    os.chdir(new_path)

    # save disease_ID as text file
    text_file = open(os.getcwd()+"/disease_id.txt", "w")
    text_file.write(input_number + ';' + disease_id + ';'+ disease_name)
    text_file.close()
    
    #build monarch network
    seedList = [input_id]

    neighbourList = get_neighbours_list(seedList)
    orthophenoList = get_orthopheno_list(seedList) 
    geneList = sum([seedList,neighbourList,orthophenoList], [])
    network = extract_edges(geneList) 

    # save network
    print_network(network, 'monarch_orthopeno_network_disease')
    
    file = 'monarch_orthopeno_network_disease_v{}.csv'.format(today)
    monarch_connections = read_connections(file)
    build_edges(monarch_connections, edges_fname = 'monarch_edges_disease')
    build_nodes(monarch_connections, nodes_fname = 'monarch_nodes_disease')

def run_monarch_symptom(input_symptom, date):
    """
    This function runs the whole Monarch script using the disease phenotype MIM number
    and symptom ID as seeds and saves nodes and edges files.

    :param input_symptom: The input symptom
    :return: nodes and edges files in /monarch folder
    """

    # turn input symptom name into ID
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    nodes_df = pd.read_csv(monarch_fname_nodes)
    symptom_id = nodes_df.loc[nodes_df['preflabel'] == input_symptom, 'id'].iloc[0]
    
    #build monarch network
    seedList = [symptom_id]

    neighbourList = get_neighbours_list(seedList)
    orthophenoList = get_orthopheno_list(seedList) 
    geneList = sum([seedList,neighbourList,orthophenoList], [])
    network = extract_edges(geneList) 

    # save network
    print_network(network, 'monarch_orthopeno_network_symptom')
    
    file = 'monarch_orthopeno_network_symptom_v{}.csv'.format(today)
    monarch_connections = read_connections(file)
    build_edges(monarch_connections, edges_fname='monarch_edges_symptom')
    build_nodes(monarch_connections, nodes_fname = 'monarch_nodes_symptom')    



if __name__ == '__main__':
    
    input_number = '143100' # input of user
    run_monarch(input_number)



