### IMPORT LIBRARIES ###
# make sure libraries are installed on your PC
# install libraries via 'pip install xxx'

from selenium import webdriver
import time, requests
import pandas as pd

from selenium.webdriver.support.ui import WebDriverWait

### DEFINE VARIABLES ###

smiles = 'CN1CC[C@]23C4=C5C=CC(O)=C4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5'
CpdName = 'pomposide I'

driver = webdriver.Chrome('C:\\Users\\c7401370\\Desktop\\ScreeningSlaveWeb\\chromedriver.exe')

### DEFINE FUNCTIONs ###

def SwissCrawler (smiles):
    SwissUrl = 'http://www.swisstargetprediction.ch/index.php'
    platform = 'SwissTargetPrediction'
    driver.get(SwissUrl)
    SearchField = driver.find_element_by_name('smiles')
    SearchField.send_keys(smiles)
    SearchField.submit()
    time.sleep(10)
    ResultUrl = driver.current_url
    r = requests.get(ResultUrl)
    df = pd.read_html(r.text)
    df = df[0]
    cols = [col for col in df.columns if col in [df.columns[2],df.columns[5]]]
    df = df[cols]
    dfKeep = df[df.columns[1]] > 0.0999999
    df = df[dfKeep]
    df.insert(0,'compound',CpdName)
    df.insert(1,'platform',platform)
    df = df.rename(columns={"Uniprot ID": "uniprotID", "Probability*": "probability"}) # rename headers
    return df

def SEACrawler (smiles):
    SEAUrl = 'http://sea.bkslab.org/'
    platform = 'SEA'
    driver.get(SEAUrl)
    SearchField = driver.find_element_by_name('query_custom_targets_paste')
    SearchField.send_keys(smiles)
    SearchField.submit()
    ResultUrl = driver.current_url
    time.sleep(40)
    r = requests.get(ResultUrl)
    df = pd.read_html(r.text)

    df = df[0]
    cols = [col for col in df.columns if col in ['Target Key','P-Value']]
    df = df[cols] # only keep 2 columns
    df = df.drop([0],axis='rows') # first row always contains missing values
    df.insert(0,'compound',CpdName)
    df.insert(1,'platform',platform) # adding 2 columns with name of the tool and the compound
    df = df.rename(columns={"Target Key": "targetkey", "P-Value": "probability"}) # rename headers
    TargetKeys = df.targetkey.values
    genes = []
    for key in TargetKeys:
        gene = key[:-2]
        genes.append(gene)
    df = df.drop('targetkey', axis=1)
    df['gene_uniprot'] = genes #string manipulation to yield uniprot names

    driver.quit()
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
    EDurl = 'http://endocrinedisruptome.ki.si/prediction.html'
    driver.get(EDurl)
    SearchField = driver.find_element_by_xpath('//*[@id="id_smiles"]')
    SearchField.send_keys(smiles)
    SubmitButton = driver.find_elements_by_xpath('//*[@id="selected-button"]')[0]
    SubmitButton.click()
    #time.sleep(600)
    ResultUrl = driver.current_url
    #r = requests.get(ResultUrl)
    #df = pd.read_html(r.text)
    #df = df[0]
    print(ResultUrl)
    print('yessssss')

### CALL FUNCTION ###

#SwissResult = SwissCrawler(smiles)
#print(SwissResult)
#
#SEAResult = SEACrawler(smiles)
#print(SEAResult)

#SuperPredCrawler (smiles, SleepTime)

EndocrineDisruptomeCrawler (smiles)

#print (df)
