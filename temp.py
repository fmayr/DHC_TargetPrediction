import pandas as pd
import os
from sklearn import preprocessing
from selenium import webdriver
import requests
from lxml.html import fromstring


loc = 'C:\\Users\\c7401370\\Desktop\\ScreeningSlaveWeb'
os.chdir(loc)

CpdName = 'pomposide I'
platform = 'SwissTargetPrediction'

EDurl = 'http://endocrinedisruptome.ki.si/prediction.html'
EDresultUrl = 'http://endocrinedisruptome.ki.si/docking/eciaqhpuqf/'

r = requests.get(EDresultUrl) # retrieve result page
html = fromstring(r.content) # convert to html to make it searchable with XPath
extract = html.xpath('//*[@id="top"]/div/div[5]/div[1]//div[@class="row show-grid"]/div[@style="border-style:solid;border-width:1px;border-color:white; background:#2ecc71"]/strong/a/text()') # alle Rezeptoren mit dockingscore, die grün sind
df = pd.DataFrame(extract, columns=['temp']) # move extracted data into dataframe
df = df['temp'].str.split(':', n = 1, expand = True) # split column in 'target' and 'DockingScore'
df.columns = ['target','DockingScore'] # rename headers

EDequi = pd.read_csv('EndocrineDisruptomeEquivalence.txt', sep=',', header=None)
EDequi.columns = ['target','UniProt']

result = pd.merge(df, EDequi, on='target', how='left')



#print(df)
print(result)


