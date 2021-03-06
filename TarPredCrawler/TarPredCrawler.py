### IMPORT LIBRARIES ###
# make sure libraries are installed on your PC
# install libraries via 'pip install xxx'

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import UnexpectedAlertPresentException
import argparse
import requests
import numpy as np
import pandas as pd
from lxml.html import fromstring
from lxml import etree
from sklearn import preprocessing

### DEFINE FUNCTIONs ###

def SwissCrawler (smiles, CpdName):
    SwissUrl = 'http://www.swisstargetprediction.ch/index.php' # URL of prediction page
    platform = 'SwissTargetPrediction'
    driver.get(SwissUrl) # browsing to prediction page
    SearchField = driver.find_element_by_name('smiles')  # locate smiles input fied
    SearchField.send_keys(smiles) # send smiles string to smiles input field
    SearchField.submit() # submit filled form
    try:
        WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.XPATH, '//*[@id="resultTable"]/tbody'))) # wait for the result page to be loaded
        CurrUrl = driver.current_url # get url of result page
        df = pd.read_html(CurrUrl) # move result table into pandas dataframe
        df = df[0] # eliminating all but the result table
        cols = [col for col in df.columns if col in [df.columns[2],df.columns[5]]] # keeping only 2 columns
        df = df[cols]
        dfKeep = df[df.columns[1]] > 0.09999 # remove all below probability of 0.12
        df = df[dfKeep]
        df.insert(0,'compound',CpdName) # insert compound name in new column
        df.insert(1,'platform',platform) # insert platform name in new column
        df = df.rename(columns={"Uniprot ID": "uniprotID", "Probability*": "prob"}) # rename headers
        newCol = []
        def get_uniprot_name(entry):
            resp = requests.get('https://www.uniprot.org/uniprot/' + str(entry) + '.xml')
            try:
                html = fromstring (resp.content)
                Entr = html.xpath('//entry/name/text()')[0]
                return Entr
            except etree.ParserError:
                Entr = 'no_entry_found_in_uniprot'
                return Entr
        for entry in df.uniprotID.values:
            if entry.count(' ') == 0:
                Entr = get_uniprot_name(entry)
                newCol.append(Entr)
                #print(entry)
            else:
                new_lst = entry.split(' ')
                new_Entr = []
                for i in new_lst:
                    Entr = get_uniprot_name(i)
                    new_Entr.append(Entr)
                new = '|'.join(new_Entr)
                newCol.append(new)
                #print(new)
        #print(len(newCol))
        #print(df.shape)
        df['UniProt_name'] = newCol
        df = df.drop('uniprotID', axis=1) # UniProt entry names assigned and entry numbers dropped
        return df
    except TimeoutException:
        CurrUrl = driver.current_url
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        temp = [CpdName, platform, driver.current_url]
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'result page reached timeout',
            'prob':CurrUrl,},
            ignore_index=True)
        return df
    except UnexpectedAlertPresentException:
        alert = driver.switch_to.alert
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'error message',
            'prob':alert.text},
            ignore_index=True)
        alert.accept()
        return df

def SEACrawler (smiles, CpdName):
    SEAUrl = 'http://sea.bkslab.org/' # URL of prediction page
    platform = 'SEA'
    driver.get(SEAUrl) # browsing to prediction page
    SearchField = driver.find_element_by_name('query_custom_targets_paste')  # locate smiles input fied
    SearchField.send_keys(smiles) # send smiles string to smiles input field
    SearchField.submit() # submit filled form
    try:
        WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.XPATH, '/html/body/div/div/table/tbody'))) # wait for the result page to be loaded
        CurrUrl = driver.current_url # get url of result page
        df = pd.read_html(CurrUrl) # move result table into pandas dataframe
        df = df[0] # eliminating all but the result table
        cols = [col for col in df.columns if col in ['Target Key','P-Value']] # keeping only 2 columns
        df = df[cols] # only keep 2 columns
        df = df.drop([0],axis='rows') # first row always contains missing values
        df.insert(0,'compound',CpdName) # insert compound name in new column
        df.insert(1,'platform',platform) # insert platform name in new column
        df = df.rename(columns={"Target Key": "targetkey", "P-Value": "prob"}) # rename headers
        TargetKeys = df.targetkey.values
        genes = []
        for key in TargetKeys:
            gene = key[:-2]
            genes.append(gene)
        df = df.drop('targetkey', axis=1)
        df['UniProt_name'] = genes #string manipulation to yield uniprot names
        return df
    except TimeoutException:
        CurrUrl = driver.current_url
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        temp = [CpdName, platform, driver.current_url]
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'result page reached timeout',
            'prob':CurrUrl,},
            ignore_index=True)
        return df
    except UnexpectedAlertPresentException:
        alert = driver.switch_to.alert
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'error message',
            'prob':alert.text},
            ignore_index=True)
        alert.accept()
        return df
    #finally:
    #    return df

def SuperPredCrawler (smiles, CpdName):
    platform = 'SuperPred'
    SPUrl = 'http://prediction.charite.de/index.php?site=chemdoodle_search_target'
    driver.get(SPUrl)
    SearchField = driver.find_element_by_xpath('//*[@id="contentwrapper"]/div/div[1]/div[2]/form/table/tbody/tr[2]/td[2]/input[1]')
    SearchField.send_keys(smiles)
    smiles_button = driver.find_elements_by_xpath('//*[@id="contentwrapper"]/div/div[1]/div[2]/form/table/tbody/tr[2]/td[2]/input[2]')[0]
    smiles_button.click()
    submit_button = driver.find_element_by_xpath('//*[@id="contentwrapper"]/div/div[4]/table/tbody/tr[2]/td/input[1]')
    submit_button.click()
    try:
        WebDriverWait(driver, 60).until(EC.element_to_be_clickable((By.XPATH, '//*[@id="hits"]/center/div/form/input[@type="submit"]')))
        redirect = driver.find_element_by_xpath('//*[@id="hits"]/center/div/form/input[@type="submit"]')
        redirect.click()
        targets = driver.find_elements_by_xpath('//*[@id="contentwrapper"]/div[3]/a/div[1]/div/table/tbody/tr/td[2]/select//option')
        preds = []
        for tars in targets:
            targetsExtr = tars.get_attribute('value')
            preds.append(targetsExtr)
        df = pd.DataFrame(preds, columns=['temp'])
        if df['temp'].count() > 0:
            df = df['temp'].str.split(',', n=1, expand=True)
            df.columns = ['UniProt_name','prob']
            df.insert(0,'compound',CpdName) # inroduce two columns with compound name and name of the tool
            df.insert(1,'platform',platform)
            df = df[['compound','platform','prob','UniProt_name']]
        else:
            df.insert(0,'compound',CpdName) # inroduce two columns with compound name and name of the tool
            df.insert(1,'platform',platform)
            df.insert(2,'prob','error message')
            df.insert(3,'UniProt_name','no predictions found')
        return df
    except TimeoutException:
        CurrUrl = driver.current_url
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        temp = [CpdName, platform, driver.current_url]
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'result page reached timeout',
            'prob':CurrUrl},
            ignore_index=True)
        return df
    except UnexpectedAlertPresentException:
        alert = driver.switch_to.alert
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'error message',
            'prob':alert.text},
            ignore_index=True)
        alert.accept()
        return df
    #finally:
    #    return df

def EndocrineDisruptomeCrawler (smiles, CpdName):
    platform = 'Endocrine Disruptome'
    EDurl = 'http://endocrinedisruptome.ki.si/prediction.html'
    driver.get(EDurl)
    SearchField = driver.find_element_by_xpath('//*[@id="id_smiles"]')
    SearchField.send_keys(smiles)
    SubmitButton = driver.find_elements_by_xpath('//*[@id="selected-button"]')[0]
    SubmitButton.click()
    try:
        WebDriverWait(driver, 5).until(EC.presence_of_element_located((By.XPATH, '//*[@id="top"]/div/div[5]/div[1]'))) # wait for the result page to be loaded
        CurrUrl = driver.current_url # get url of result page
        r = requests.get(CurrUrl) # retrieve result page
        html = fromstring(r.content) # convert to html to make it searchable with XPath
        extract = html.xpath('//*[@id="top"]/div/div[5]/div[1]//div[@class="row show-grid"]/div[@style="border-style:solid;border-width:1px;border-color:white; background:#2ecc71"]/strong/a/text()')
        # this XPath yield all "green" fields
        df = pd.DataFrame(extract, columns=['temp']) # move extracted data into dataframe
        df = df['temp'].str.split(':', n = 1, expand = True) # split column in 'target' and 'DockingScore'
        df.columns = ['target','DockingScore'] # rename headers
        df.target.str.decode("utf-8") # decode bytes objects in strings using utf-8 encoding
        df.target = df.target.str.strip() # stripping whitespaces from strings (not sure if they appear all over the html table, but for AR they do...)
        df.DockingScore = df.DockingScore.str.strip()
        df.insert(0,'compound',CpdName) # inroduce two columns with compound name and name of the tool
        df.insert(1,'platform',platform)
        df.target = df.target.replace({
        'AR' : 'ANDR_HUMAN',
        'AR an.' : 'ANDR_HUMAN',
        'ER α' : 'ESR1_HUMAN',
        'ER α an.' : 'ESR1_HUMAN',
        'ER β' : 'ESR2_HUMAN',
        'ER β an.' : 'ESR2_HUMAN',
        'GR' : 'GCR_HUMAN',
        'GR an.' : 'GCR_HUMAN',
        'LXR α' : 'NR1H3_HUMAN',
        'LXR β' : 'NR1H2_HUMAN',
        'PPAR α' : 'PPARA_HUMAN',
        'PPAR β' : 'PPARD_HUMAN',
        'PPAR γ' : 'PPARG_HUMAN',
        'RXR α' : 'RXRA_HUMAN',
        'TR α' : 'THA_HUMAN',
        'TR β' : 'THB_HUMAN'
        }) # changes target names to UniProt nomenclature
        df = df.rename(columns={'target': 'UniProt_name', 'DockingScore':'prob'}) # rename headers
    except TimeoutException:
        CurrUrl = driver.current_url
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        temp = [CpdName, platform, driver.current_url]
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'result page reached timeout',
            'prob':CurrUrl,},
            ignore_index=True)
    except UnexpectedAlertPresentException:
        alert = driver.switch_to.alert
        df = pd.DataFrame(columns=['compound','platform','UniProt_name','prob'])
        df = df.append({'compound':CpdName, # creates row with name; platform; url of result table
            'platform':platform,
            'UniProt_name':'error message',
            'prob':alert.text},
            ignore_index=True)
        alert.accept()
        return df
    #finally:
    #    return df

### INITIALIZE SCRIPT ###

if __name__ == '__main__':
    print('\n\n')
    print('TarPredCrawler initialized...')
    print('\n\n')

    ### PROCESS INPUT ###

    parser = argparse.ArgumentParser(description='Crawl throug 4 Target Prediction Servers')
    parser.add_argument('-in',
        '--input',
        type=str,
        metavar='',
        required=True,
        help='csv-table in the format "name ; smiles-code" of n compounds')
    parser.add_argument('-out',
        '--output',
        type=str,
        metavar='',
        required=True,
        help='csv-table populated with processed results')
    args = parser.parse_args()

    ### READ IN INPUT ###

    with open (args.input, 'r') as fin:
        data = pd.read_csv(fin, sep=';', names=['name','smiles'])
    rowcount = data['name'].count()
    print('     Found {} molecules in "{}"\n'.format(rowcount, args.input))
    print('     Start screening:')

    ### START CRAWLING ###

    options = webdriver.ChromeOptions()
    options.add_experimental_option('excludeSwitches', ['enable-logging'])
    driver = webdriver.Chrome('C:\\Users\\c7401370\\Documents\\GIT\\chromedriver.exe', options=options)

    cols = ['compound','platform','prob','UniProt_name']

    results = pd.DataFrame(columns=cols)
    results.to_csv(args.output,sep=';') # write empty dataframe as csv file

    for index, row in data.iterrows():
        CpdName = row['name']
        smiles = row['smiles']
        SwissResult = SwissCrawler(smiles, CpdName)
        SEAResult = SEACrawler(smiles, CpdName)
        SuperPredResult = SuperPredCrawler(smiles, CpdName)
        #EndocrineDisruptomeResult = EndocrineDisruptomeCrawler(smiles, CpdName)
        with open (args.output,'a',newline='') as f:
            SwissResult.to_csv(f,sep=';',header=False)
            SEAResult.to_csv(f,sep=';',header=False)
            SuperPredResult.to_csv(f,sep=';',header=False)
            #EndocrineDisruptomeResult.to_csv(f,sep=';',header=False)
        print('         screened {} of {} molecules ({})'.format(index+1, rowcount, CpdName))

    driver.quit()
    print('')
    print('     Finished Analysis')
    print('     Results are now available in "{}"'.format(args.output))
