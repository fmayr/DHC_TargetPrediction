### IMPORT LIBRARIES ###
# make sure libraries are installed on your PC
# install libraries via 'pip install xxx'

from selenium import webdriver
import time, requests
import pandas as pd
from lxml.html import fromstring
from selenium.webdriver.support.ui import WebDriverWait

### DEFINE VARIABLES ###

smiles = 'CN1CC[C@]23C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5'
CpdName = 'pomposide I'

driver = webdriver.Chrome('C:\\Users\\c7401370\\Desktop\\ScreeningSlaveWeb\\chromedriver.exe')

### DEFINE FUNCTIONs ###

def SwissCrawler (smiles):
    SwissUrl = 'http://www.swisstargetprediction.ch/index.php' # URL of prediction page
    platform = 'SwissTargetPrediction'
    driver.get(SwissUrl) # browsing to prediction page
    SearchField = driver.find_element_by_name('smiles')  # locate smiles input fied
    SearchField.send_keys(smiles) # send smiles string to smiles input field
    SearchField.submit() # submit filled form
    time.sleep(10)
    ResultUrl = driver.current_url
    r = requests.get(ResultUrl) # get-request on results page
    df = pd.read_html(r.text) # move result table into pandas dataframe
    df = df[0] # eliminating all but the result table
    cols = [col for col in df.columns if col in [df.columns[2],df.columns[5]]] # keeping only 2 columns
    df = df[cols]
    dfKeep = df[df.columns[1]] > 0.0999999 # remove all below probability of 0.1
    df = df[dfKeep]
    df.insert(0,'compound',CpdName) # insert compound name in new column
    df.insert(1,'platform',platform) # insert platform name in new column
    df = df.rename(columns={"Uniprot ID": "uniprotID", "Probability*": "STP_prob"}) # rename headers
    return df

def SEACrawler (smiles):
    SEAUrl = 'http://sea.bkslab.org/'
    platform = 'SEA'
    driver.get(SEAUrl) # browsing to prediction page
    SearchField = driver.find_element_by_name('query_custom_targets_paste') # locate smiles input fied
    SearchField.send_keys(smiles) # send smiles string to smiles input field
    SearchField.submit() # submit filled form
    ResultUrl = driver.current_url
    r = requests.get(ResultUrl) # get URL of the result page after being redirect by JS script
    time.sleep(40)
    df = pd.read_html(r.text) # move result table into pandas dataframe
    df = df[0] # eliminating all but the result table
    cols = [col for col in df.columns if col in ['Target Key','P-Value']] # keeping only 2 columns
    df = df[cols] # only keep 2 columns
    df = df.drop([0],axis='rows') # first row always contains missing values
    df.insert(0,'compound',CpdName) # insert compound name in new column
    df.insert(1,'platform',platform) # insert platform name in new column
    df = df.rename(columns={"Target Key": "targetkey", "P-Value": "SEA_prob"}) # rename headers
    TargetKeys = df.targetkey.values
    genes = []
    for key in TargetKeys:
        gene = key[:-2]
        genes.append(gene)
    df = df.drop('targetkey', axis=1)
    df['gene_uniprot'] = genes #string manipulation to yield uniprot names

    #driver.quit()
    return df

def SuperPredCrawler (smiles, SleepTime):
    SPUrl = 'http://prediction.charite.de/index.php?site=chemdoodle_search_target'
    SPUrlResult = 'http://prediction.charite.de/index.php?site=pred_results'
    driver.get(SPUrl)
    SearchField = driver.find_element_by_xpath('//*[@id="contentwrapper"]/div/div[1]/div[2]/form/table/tbody/tr[2]/td[2]/input[1]')
    SearchField.send_keys(smiles)
    smiles_button = driver.find_elements_by_xpath('//*[@id="contentwrapper"]/div/div[1]/div[2]/form/table/tbody/tr[2]/td[2]/input[2]')[0]
    smiles_button.click()
    time.sleep(50)
    r = requests.get(SPUrlResult)

    #submit_button = driver.find_elements_by_name('start')
    #submit_button.click()
    #time.sleep(15)
    #submit_button.click()

    cookies = {'enwiki_session': '17ab96bd8ffbe8ca58a78657a918558'}
    r = requests.post('http://wikipedia.org', cookies=cookies)

    print(r.text)

    print('yess')




    print('************************************************')

def EndocrineDisruptomeCrawler (smiles):
    platform = 'Endocrine Disruptome'
    EDurl = 'http://endocrinedisruptome.ki.si/prediction.html'
    driver.get(EDurl)
    SearchField = driver.find_element_by_xpath('//*[@id="id_smiles"]')
    SearchField.send_keys(smiles)
    SubmitButton = driver.find_elements_by_xpath('//*[@id="selected-button"]')[0]
    SubmitButton.click()
    EDresultUrl = driver.current_url
    time.sleep(300)
    print('finished waiting')
    r = requests.get(EDresultUrl) # retrieve result page
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

    df = df.rename(columns={'target': 'gene_uniprot', 'DockingScore':'ED_prob'}) # rename headers
    driver.quit()
    return df

### CALL FUNCTION ###

#SwissResult = SwissCrawler(smiles)
#print(SwissResult)
#
#SEAResult = SEACrawler(smiles)
#print(SEAResult)

#SuperPredCrawler (smiles, SleepTime)

EndocrineDisruptomeResult = EndocrineDisruptomeCrawler (smiles)
print (EndocrineDisruptomeResult)
