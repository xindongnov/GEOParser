import os
import sys
import csv
import time
import argparse
import requests
import traceback
import pandas as pd
import urllib.request
import xml.etree.ElementTree as ET

from math import ceil
from Bio import Entrez
from Bio import Medline
from bs4 import BeautifulSoup
from fake_useragent import UserAgent


def getXML_Paper(accession):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=xml"% accession
    try:
        with urllib.request.urlopen(URL) as response:
            content = response.read()               
            docString = content.decode('utf-8')
        soup = BeautifulSoup(docString, 'xml')
        article_title = soup.find('ArticleTitle').get_text()
    except:
        print("Network fluctuations. Wait for 5s.")
        time.sleep(5)
        try:
            with urllib.request.urlopen(URL) as response:
                content = response.read()               
                docString = content.decode('utf-8')
            soup = BeautifulSoup(docString, 'xml')
            article_title = soup.find('ArticleTitle').get_text()
        except:
                print("Exception in user code: XML")
                print(f"{accession} is error")
                print('-' * 60)
                traceback.print_exc(file=sys.stdout)
                print('-' * 60)
    return article_title

def download_paper(pmc_id,save_path):
    ua = UserAgent()
    headers={"User-Agent":ua.random,'Connection':'close'}
    sess = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries = 20)
    sess.mount('https://', adapter)
    link = 'https://pmc.ncbi.nlm.nih.gov/articles/' + pmc_id
    '''using session instead to avoid for SSL EOF error'''
    try:
        response = sess.get(url=link, headers=headers, params={'param':'1'}).text
        '''avoid for some link without pdf'''
        soup = BeautifulSoup(response, 'html.parser')
        pdf_link = soup.find('meta', attrs={'name': 'citation_pdf_url'})['content']
        if pdf_link:
            try :
                pdf = sess.get(url=pdf_link, headers=headers, params={'param':'1'}, stream="TRUE")
                with open(save_path, mode='wb') as f:
                    for data in pdf.iter_content():
                        f.write(data)
            except:
                time.sleep(5)
                print("Network fluctuation, please rest for 5 seconds.")
                pdf = sess.get(url=pdf_link, headers=headers, params={'param':'1'}, stream="TRUE")
                with open(save_path, mode='wb') as f:
                    for data in pdf.iter_content():
                        f.write(data)
    except:
        time.sleep(5)
        print("Network fluctuation, please rest for 5 seconds.")
        response = sess.get(url=link, headers=headers, params={'param':'1'}).text
        '''avoid for some link without pdf'''
        soup = BeautifulSoup(response, 'html.parser')
        pdf_link = soup.find('meta', attrs={'name': 'citation_pdf_url'})['content']
        if pdf_link:
            try:
                pdf = sess.get(url=pdf_link, headers=headers, params={'param':'1'}, stream="TRUE")
                with open(save_path, mode='wb') as f:
                    for data in pdf.iter_content():
                        f.write(data)
            except:
                time.sleep(5)
                print("Network fluctuation, please rest for 5 seconds.")
                pdf = sess.get(url=pdf_link, headers=headers, params={'param':'1'}, stream="TRUE")
                with open(save_path, mode='wb') as f:
                    for data in pdf.iter_content():
                        f.write(data)
    return pdf_link

def get_species(gseid):
    ua = UserAgent()
    headers={"User-Agent":ua.random,'Connection':'close'}
    sess = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries = 20)
    sess.mount('https://', adapter)
    link = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + gseid
    '''using session instead to avoid for SSL EOF error'''
    try:
        response = sess.get(url=link, headers=headers, params={'param':'1'}).text
        '''avoid for some link without pdf'''
        soup = BeautifulSoup(response, 'html.parser')
        a_tags = soup.find_all('a')
        # 提取每个<a>标签的文本内容
        for a_tag in a_tags:
            # 检查href属性是否符合特定模式
            if '/Taxonomy/Browser/wwwtax.cgi?mode=Info' in a_tag['href']:
                organism_tag = a_tag.text.strip()
    except:
        time.sleep(5)
        print("Network fluctuation, please rest for 5 seconds.")
        response = sess.get(url=link, headers=headers, params={'param':'1'}).text
        '''avoid for some link without pdf'''
        soup = BeautifulSoup(response, 'html.parser')
        a_tags = soup.find_all('a')
        # 提取每个<a>标签的文本内容
        for a_tag in a_tags:
            # 检查href属性是否符合特定模式
            if '/Taxonomy/Browser/wwwtax.cgi?mode=Info' in a_tag['href']:
                organism_tag = a_tag.text.strip()
    return organism_tag
                


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pmid_file',default="/fse/home/wangyongyan/parser_tail/GSE_meta_output.csv" ,type=str,help='Files containing PMIDs')
    parser.add_argument('-np', '--name_of_pmid_column', type=str,default="GSE_Pubmed_ID",required=False, help='The column name of the column in which the pmid is located in the csv file')
    parser.add_argument('-n', '--name_of_GSE_column', type=str,default="GSE",required=False, help='The column name of the column in which the pmid is located in the csv file')
    parser.add_argument('-f', '--folder_path', type=str, default = 'load_paper',help='path to store paper')
    parser.add_argument('-o', '--paper_output_csv', type=str, default='paper_output.csv',help='path to some parts of meta')
    parser.add_argument('-e', '--paper_error_csv', type=str, default='paper_error.csv',help='path to error of meta')


    args = parser.parse_args()
    pmid_file = args.pmid_file
    name_of_pmid_column = args.name_of_pmid_column
    name_of_GSE_column = args.name_of_GSE_column
    folder_path = args.folder_path
    paper_output_csv = args.paper_output_csv
    table = pd.read_csv(pmid_file, sep=",")
    series_list = table[name_of_pmid_column].tolist()
    result = table[table[name_of_pmid_column].isin(series_list)][[name_of_GSE_column, name_of_pmid_column]]
    paper_error_csv = args.paper_error_csv

    with open(paper_output_csv, 'a', newline='') as csvfile:
        fieldnames = ["GSE", 'PMID', 'Paper','PMC', 'DOI','Journal', 'Sub_Journal', "Paper_load_url"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader() 
    if not os.path.exists(folder_path):
            os.makedirs(folder_path)
    for i in result.index:
        start = time.time()
        gseid = result.loc[i, name_of_GSE_column]
        pmid = result.loc[i, name_of_pmid_column]
        if isinstance(pmid, str):
            pmid = pmid.strip("[]")
            pmid = pmid.strip("'")
            if ',' in pmid:
                pmid_list = pmid.split("', '")
                pmid_list = [id.strip("'") for id in pmid_list]
            else :
                pmid_list = [pmid]
        else:
            print("{} is no PMC".format(gseid))
            continue
        ua = UserAgent()
        for pmid in pmid_list :
            try:
                handle = Entrez.efetch(db="pubmed", id=pmid)
                # 解析XML数据
                tree = ET.parse(handle)
            except:
                handle = Entrez.efetch(db="pubmed", id=pmid)
                # 解析XML数据
                tree = ET.parse(handle)
            pmc_id = None
            article_info = {} 
            article_info["ArticleTitle"] = getXML_Paper(pmid)
            for elem in tree.iter():
                if elem.tag == "ELocationID":
                    article_info["ELocationID"] = elem.text
                if elem.tag == "Title":
                    article_info["Title"] = elem.text
                if elem.tag == "ISOAbbreviation":
                    article_info["ISOAbbreviation"] = elem.text
            for article_id in tree.iter('ArticleId'):
                    if article_id.attrib.get('IdType') == 'pmc':
                        if 'PMC' in article_id.text:
                            pmc_id = article_id.text
                            break  # 找到第一个符合条件的即可
            if pmc_id is not None :
                try:
                    species = get_species(gseid)
                    save_path = folder_path + '/' + '{}.pdf'.format(pmid)
                    pdf_link = download_paper(pmc_id, save_path)
                    with open(paper_output_csv, 'a', newline='') as csvfile:
                        fieldnames = ["GSE", 'Species', 'PMID', 'Paper','PMC', 'DOI','Journal', 'Sub_Journal', "Paper_load_url"]
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writerow({'GSE': gseid, 'Species': species,'PMID': pmid, 'Paper': article_info["ArticleTitle"], 
                                        'PMC':pmc_id, 'DOI': article_info["ELocationID"], 'Journal':article_info["Title"],
                                        'Sub_Journal':article_info["ISOAbbreviation"], 'Paper_load_url':pdf_link}) 
                    print('{} is finished.'.format(gseid))
                except:
                    with open(paper_error_csv, 'a', newline='') as csvfile:
                        fieldnames = ["GSE", 'PMID']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writerow({'GSE': gseid, 'PMID': pmid})
                    print('{} is error.'.format(gseid))
            end = time.time()
            print('Total time:', end - start)

if __name__ == "__main__":
    main()          

            