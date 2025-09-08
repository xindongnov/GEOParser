import random
import math
import urllib.request
import sys
import os
from bs4 import BeautifulSoup
from Bio import Medline
import pandas as pd
import requests
import re
import socket
import time
import xml.etree.ElementTree as ET
from Bio import Entrez
from math import ceil
import argparse
import csv
import tarfile

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
            for elem in tree.iter():
                if elem.tag == "ArticleTitle":
                    article_info["ArticleTitle"] = elem.text
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
            pmc_url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=%s' %(pmc_id)
            link_href = None
            # 发起GET请求获取文章内容
            try:
                response = requests.get(pmc_url)
                if response.status_code == 200:
                # 获取文章内容
                    article_content = response.content
                    root = ET.fromstring(article_content)
                    # 查找链接（href 属性）
                    for record in root.findall('.//record'):
                        link_href = record.find('link').get('href')
            except:
                time.sleep(5)
                print("Network disturbances, try again.")
                response = requests.get(pmc_url)
                if response.status_code == 200:
                # 获取文章内容
                    article_content = response.content
                    root = ET.fromstring(article_content)
                    # 查找链接（href 属性）
                    for record in root.findall('.//record'):
                        link_href = record.find('link').get('href')
            if link_href:
                if not os.path.exists(os.path.join(folder_path, gseid)):
                    os.makedirs(os.path.join(folder_path, gseid))
                file_path = os.path.join(folder_path, gseid)
                file_name = link_href.split('/')[-1]
                os.system('wget %s -O %s/%s' % (link_href, file_path, file_name))
                # os.system(f"wget {link_href} -O {file_path}")
                print("{}'s paper is download.".format(gseid))
                with open(paper_output_csv, 'a', newline='') as csvfile:
                    fieldnames = ["GSE", 'PMID', 'Paper','PMC', 'DOI','Journal', 'Sub_Journal', "Paper_load_url"]
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writerow({'GSE': gseid, 'PMID': pmid, 'Paper': article_info["ArticleTitle"], 
                                    'PMC':pmc_id, 'DOI': article_info["ELocationID"], 'Journal':article_info["Title"],
                                    'Sub_Journal':article_info["ISOAbbreviation"], 'Paper_load_url': link_href}) 
                ##解压
                file_path_2 = file_path +'/'+ file_name
                try:
                    with tarfile.open(file_path_2, 'r:gz') as tar:
                        tar.extractall(path=file_path)
                        print(f"File {file_path_2} has been extracted to {file_path}")
                    if os.path.exists(file_path_2):
                        os.remove(file_path_2)
                        print(f"File {file_path_2} has been deleted.")
                    directory_path = file_path +'/'+ pmc_id # 替换为你的目录路径
                    # 遍历目录中的所有文件
                    for filename in os.listdir(directory_path):
                        # 获取文件的完整路径
                        file_path = os.path.join(directory_path, filename)
                        # 检查是否为文件
                        if os.path.isfile(file_path):
                            # 获取文件扩展名
                            file_extension = os.path.splitext(filename)[1].lower()
                            # 如果文件扩展名不是.pdf，则删除该文件
                            if file_extension != '.pdf':
                                os.remove(file_path)
                                print(f"Deleted: {file_path}")
                            else:
                                print(f"Kept: {file_path}")
                    time.sleep(3)
                except:
                    with open(paper_error_csv, 'a', newline='') as csvfile:
                        fieldnames = ["GSE"]
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writerow({'GSE': gseid})
            else:
                print("{}'s paper is not free.".format(gseid))
                time.sleep(3)
            end = time.time()
            print('Total time:', end - start)
    
    
if __name__ == "__main__":
    main()
    
