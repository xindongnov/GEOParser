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



def getXML_Paper(pmid):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&retmode=xml"% pmid
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
                print(f"{pmid} is error")
                print('-' * 60)
                traceback.print_exc(file=sys.stdout)
                print('-' * 60)
    return article_title


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
        fieldnames = ["GSE", 'Species', 'PMID', 'Paper','PMC', 'DOI','Journal', 'Sub_Journal', "Paper_load_url"]
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
        # 获取PMID为pmid的论文，返回medline格式的文本
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
                    with open(paper_output_csv, 'a', newline='') as csvfile:
                        fieldnames = ["GSE", 'Species', 'PMID', 'Paper','PMC', 'DOI','Journal', 'Sub_Journal', "Paper_load_url"]
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writerow({'GSE': gseid, 'Species': species,'PMID': pmid, 'Paper': article_info["ArticleTitle"], 
                                        'PMC':pmc_id, 'DOI': article_info["ELocationID"], 'Journal':article_info["Title"],
                                        'Sub_Journal':article_info["ISOAbbreviation"], 'Paper_load_url':None}) 
                    print('{} is finished.'.format(gseid))
                except:
                    with open(paper_error_csv, 'a', newline='') as csvfile:
                        fieldnames = ["GSE", 'PMID']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writerow({'GSE': gseid, 'PMID':pmid})
                    print('{}_{} is error.'.format(gseid, pmc_id))
                ##解压
            time.sleep(3)
            end = time.time()
            print('Total time:', end - start)
    
    
if __name__ == "__main__":
    main()


## read pdf  
# search_text = "Genome-Scale Oscillations in DNA Methylation during Exit from Pluripotency"
# pdf_file = '/fse/home/wangyongyan/parser_tail/load_paper/GSE75973/PMC6066359/main.pdf'
# with open(pdf_file, 'rb') as file:
#         # 创建PDF阅读器对象
#         pdf_reader = PyPDF2.PdfReader(file)
#         # 遍历PDF的每一页
#         for page_num in range(len(pdf_reader.pages)):
#             # 获取当前页的文本内容
#             page_text = pdf_reader.pages[page_num].extract_text()        
#             # 检查搜索文本是否在当前页的文本中
#             if search_text in page_text:
#                 print(f"Found '{search_text}' in {pdf_file} on page {page_num + 1}")
#                 break
#             else:
#                 print(f"Did not find '{search_text}' in {pdf_file}")
