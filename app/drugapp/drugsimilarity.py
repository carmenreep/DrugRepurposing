# -*- coding: utf-8 -*-
"""
Module to create the drug-drug similarity graph.

Created on Thu Jun 16 08:21:14 2022

@author: Carmen Reep
"""

from biothings_client import get_client
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger, DataStructs
import datetime
import logging
import numpy as np
import os
import pandas as pd
import requests

today = datetime.date.today()

# Mute RDKit logger
RDLogger.logger().setLevel(RDLogger.CRITICAL)
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def drug_smiles_conversion(dgidb_nodes_file_str):
    """
    This function performs drug ID conversion from chembl and wikidata ID to smiles chemical structure notation.
    :param dgidb_nodes_file_str: string of file name containing all drugs to converse
    :return: concept dictionary (where keys are smiles structures and values are old ontology IDs)
    """

    print('\n Mapping DGIdb drug ontologies to smiles structures ...')
    
    # open csv file
    dgidb_nodes_file_str = './DGIdb/'+ dgidb_nodes_file_str
    dgidb_nodes_df = pd.read_csv(dgidb_nodes_file_str)
    
    drug_lst = dgidb_nodes_df.loc[dgidb_nodes_df['semantic_groups'] == 'drug', 'id']
    
    # final dict
    concept_dct = dict()

    # input
    chembl = list()
    wikidata = list()
    for idx in drug_lst:
        if ':' in idx:
            if 'chembl' in idx.split(':')[0].lower():
                chembl.append(idx.split(':')[1])
            elif 'wikidata' in idx.split(':')[0].lower():
                wikidata.append(idx.split(':')[1])

    chembl = list(set(chembl))
    wikidata = list(set(wikidata))

    # api call
    mc = get_client('chem')

    df_chembl = mc.querymany(qterms=chembl, scopes=['chembl.molecule_chembl_id'], size=1, fields=['chembl.smiles'],as_dataframe=True) #9 out of1586 not found
    df_wikidata = mc.querymany(qterms=wikidata, fields='chembl.smiles', size=1, as_dataframe=True) #8 out of 8 not found

    #build dictionaries
    ids_chembl = df_chembl.reset_index().copy()
    if 'chembl.smiles' in ids_chembl.columns:
        # turn ids into uris
        ids_chembl['query'] = 'chembl:' + ids_chembl['query']
        smiles2chembl = dict(zip(ids_chembl['chembl.smiles'], ids_chembl['query'])) 
        concept_dct.update(smiles2chembl)
    else:
        print('no chembl IDs can be mapped to smiles')

    ids_wikidata = df_wikidata.reset_index().copy()
    if 'chembl.smiles' in ids_wikidata.columns:
        # turn ids into uris
        ids_chembl['query'] = 'wikidata:' + ids_chembl['query']
        smiles2wikidata = dict(zip(ids_wikidata['chembl.smiles'],ids_wikidata['query']))
        concept_dct.update(smiles2wikidata)
    else:
        print('no wikidata IDs can be mapped to smiles')

    # remove drugs for which no smiles was found (value is nan)
    concept_dct = {k: concept_dct[k] for k in concept_dct if not pd.isna(concept_dct[k])}
    return concept_dct

def get_mols(smiles_list):
    """
    This function converts a list of smiles into a list of RDKit molecules 
    :param smiles_list: the list of smiles strings
    :return: a list of RDKit molecules
    """   
    for i in smiles_list:
        try:
            mol = Chem.MolFromSmiles(i) 
            if mol is not None:
                yield mol
                
        except Exception as e:
            logger.warning(e)

def get_fingerprints(mols, radius=2, length=4096):
    """
    This function converts molecules to ECFP bitvectors.
    :param mols: RDKit molecules
    :param radius: ECFP fingerprint radius
    :param length: number of bits
    :return: a list of fingerprints
    """
    return [AllChem.GetMorganFingerprintAsBitVect(m, radius, length) for m in mols]


def calculate_internal_pairwise_similarities(smiles_list):
    """
    This function computes the pairwise similarities of the provided list of smiles against itself.
    It calculates the molecular similarity using Tanimoto coefficient.
    (the higher the coefficent, the more similar the two chemicals are based on molecular structure)
    :param smiles_list: the list of smiles
    :return: symmetric matrix of pairwise similarities. Diagonal is set to zero. 
    """
    if len(smiles_list) > 10000:
        logger.warning(f'Calculating internal similarity on large set of '
                        f'SMILES strings ({len(smiles_list)})')

    mols = get_mols(smiles_list)
    fps = get_fingerprints(mols)
    nfps = len(fps)

    similarities = np.zeros((nfps, nfps))

    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        similarities[i, :i] = sims
        similarities[:i, i] = sims

    return similarities 

def create_drug_sim_graph(smiles_conc_dict, K=10):
    '''This function creates the KNN similarity drug graph using Tanimo similarity on smiles structures 
    and saves it as a csv file.
    :param smiles_conc_dict: dictionary that converts ensembl IDs to smiles strings
    :param K: the number of nearest neighbours to select
    :return: csv file with columns subj(=drug1), obj (=drug2), weight (=similarity score)
    '''
    # get list of all smiles strings in the dictionary
    smiles_lst = list(smiles_conc_dict.keys())
    # remove nan values
    smiles_lst = [x for x in smiles_lst if pd.isnull(x) == False] 
    # convert list of smiles to ids using the smiles_concept_dict
    id_lst = list(smiles_conc_dict.get(item,item)  for item in smiles_lst) 
    
    # calculate the similarity between all molecules in this list 
    similarity_matr = calculate_internal_pairwise_similarities(smiles_lst)
    similarity_matr = similarity_matr * -1 # to make sure most similar is smallest
    
    sims_np_sort = np.argpartition(similarity_matr, K, axis=1) #indices of smallest similarities come first
    
    drug_similarity_graph = pd.DataFrame(columns=['subject_id','property_id', 'object_id','reference_uri',
                                                  'reference_supporting_text', 'reference_date', 'property_label',
                                                  'property_description', 'property_uri'])
    
    for i in range(similarity_matr.shape[0]):
        for j in sims_np_sort[i, :K]:
            subj = id_lst[i]
            obj = id_lst[j]
            weight = similarity_matr[i,j] *-1
        
            prop = 'CHEMINF:000481'
            prop_label = 'similar to'
            prop_uri = 'http://semanticscience.org/resource/CHEMINF_000481' # similar to in chemical inf ontology
            ref_text = 'This edge comes from calculating the Tanimo similarity score on the smiles structures'
            if weight !=0:
                new_row = {'subject_id':subj, 'property_id': prop, 'object_id':obj, 'reference_uri': 'NA',
                           'reference_supporting_text': ref_text, 'reference_date':today, 
                           'property_label': prop_label, 'property_description': 'NA', 'property_uri': prop_uri}
                drug_similarity_graph = drug_similarity_graph.append(new_row, ignore_index=True)
    
    # save output file
    path = os.getcwd() + '/similaritygraph'
    if not os.path.isdir(path): os.makedirs(path)
    drug_similarity_graph.to_csv('{}/{}_v{}.csv'.format(path, 'drugdrugsim', today), index=False)
    
    return drug_similarity_graph

def run_drugsimilarity():
    """
    This function runs the drugsimilarity script and saves nodes and edges files.

    :return: nodes and edges files in /similaritygraph folder
    """

    # create drug similarity graph
    dgidb_nodes_file_str = 'DGIdb_nodes_v{}.csv'.format(today) # drug file
    smiles_concept_dict = drug_smiles_conversion(dgidb_nodes_file_str)
    # get and save the drug similarity graph
    create_drug_sim_graph(smiles_concept_dict, K=10) 


if __name__ == '__main__':
    run_drugsimilarity()

 