from requests_html import HTMLSession
from lxml.html import fromstring
import requests, re, os, argparse, time, json
import urllib.parse
import urllib.request
import pandas as pd
import numpy as np

def get_cids(smiles):
    session = HTMLSession()
    new = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' + smiles + '/cids/XML'
    response = session.get(new)
    cid = response.html.find('CID')
    try:
        return cid[0].text
    except IndexError:
        return 'no entry found'

def get_assay_ids(cid):
    session = HTMLSession()
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/' + cid + '/aids/XML?aids_type=active'
    response = session.get(url)
    assay_ids = response.html.find('AID')
    if len(assay_ids) > 0:
        aids = []
        for aid in assay_ids:
            aids.append(aid.text)
        df = pd.DataFrame({'aid' : aids})
        return df
    else:
        df = pd.DataFrame()
        return df

def get_gene_name(aid):
    session = HTMLSession()
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/' + aid + '/targets/ProteinGI,ProteinName,GeneID,GeneSymbol/XML'
    response = session.get(url)
    try:
        gene_name = response.html.find('GeneSymbol')[0].text
        return gene_name
    except IndexError:
        return 'no entry found'

def get_uniprot_entry(gene_name):
    session = HTMLSession()
    uniprot = 'https://www.uniprot.org/uniprot/?query='
    organism = '&fil=organism%3A"Homo+sapiens+%28Human%29+%5B9606%5D"'
    sort = '&sort=score'
    cols = '&columns=id,entry name,genes,protein names&format=tab'
    response = requests.get(uniprot+gene_name+organism+sort+cols).text
    rows = re.split('\n',response)
    try:
        row = re.split('\t',rows[1])
        uniprot_name = row[1]
        uniprot_title = row[3]
        return uniprot_name, uniprot_title
    except IndexError:
        return 'no entry found'

def get_reactome_data(uniprot_id):
    uniprot_url = 'https://www.uniprot.org/uploadlists/'
    reactome_url = 'https://reactome.org/ContentService/data/discover/'
    params = {
    'from': 'ACC+ID',
    'to': 'REACTOME_ID',
    'format': 'tab',
    'query': uniprot_id
    }
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(uniprot_url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read().decode('utf-8')
    rows = re.split('\n',response)[1:-1]
    items = []
    for row in rows:
        item = re.split('\t',row)[1]
        items.append(item)
    reactome_names = []
    for reactom_id in items:
        req = urllib.request.Request(reactome_url+reactom_id)
        with urllib.request.urlopen(req) as f:
            resp = f.read().decode('utf-8')
            res = json.loads(resp)
        reactome_name = res['name']
        reactome_names.append(reactome_name)
    return reactome_names

if __name__ == '__main__':

    start = time.time()
    parser = argparse.ArgumentParser(description='generate bioactivity graph')
    parser.add_argument('-in',
        '--input',
        type=str,
        metavar='',
        required=True,
        help='csv-table in the format "name ; PubChem compound ID" of n compounds')
    parser.add_argument('-out',
        '--output',
        type=str,
        metavar='',
        required=True,
        help='sif-file input for Cytoscape')
    args = parser.parse_args()
    with open (args.input, 'r') as fin:
        data = pd.read_csv(fin, sep=';', names=['name','smiles'])
    rowcount = data['name'].count()
    print('     Found {} molecules in "{}"\n'.format(rowcount, args.input))
    df = pd.DataFrame(columns=['compound_name','aid','gene_name'])
    for index, row in data.iterrows():
        compound = row['name']
        smiles = row['smiles']
        cid = get_cids(smiles)
        if cid != 'no entry found':
            bioactivity_data = get_assay_ids(cid)
            if not bioactivity_data.empty:
                bioactivity_data['gene_name'] = bioactivity_data['aid'].apply(get_gene_name)
                bioactivity_data = bioactivity_data[bioactivity_data['gene_name'].notnull()]
                bioactivity_data['temp'] = bioactivity_data['gene_name'].apply(get_uniprot_entry)
                bioactivity_data = bioactivity_data.assign(**pd.DataFrame(bioactivity_data['temp'].values.tolist()).add_prefix('val_'))
                bioactivity_data = bioactivity_data[bioactivity_data['gene_name'] != 'no entry found']
                bioactivity_data['compound_name'] = compound
                bioactivity_data = bioactivity_data[['compound_name','aid','gene_name','val_0','val_1']]
                df = pd.concat([df,bioactivity_data],sort=True)
        print('     finished retireving literature data from {} ({} of {})'.format(compound,index+1,rowcount))
    df = df[df['val_0'].notnull()]
    df = df.groupby(['compound_name','val_0','val_1'])[['gene_name']].count()
    df['type'] = 'literature'
    df = df.reset_index(level=[0,1,2])
    df.columns = ['compound_name','uniprot_entry','protein_name','edge_weight','type']
    df = df[df['edge_weight'] != 'no entry found']
    df['reactome_pathways'] = df['uniprot_entry'].apply(get_reactome_data)
    df = df[['compound_name','type','uniprot_entry','edge_weight','protein_name','reactome_pathways']]
    df.to_csv(args.output,sep=';')
    end = time.time()
    print('')
    print('     find your results in "{}"'.format(args.output))
    print('     execution time was {} minutes'.format(round((end-start)/60, 2)))
