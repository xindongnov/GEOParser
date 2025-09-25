import argparse
import parser_function as p_f
import csv
import os
import pandas as pd
import time
import re
from ruamel.yaml import YAML
from tqdm import tqdm

yaml=YAML(typ='safe')


def main():
    parser = argparse.ArgumentParser()
    # data
    parser.add_argument('-o', '--gse_output_csv', type=str, default='GSE_meta_output.csv',
                        help='path to some parts of meta')
    parser.add_argument('-g', '--gse_gsm_output_csv', type=str, default='GSE_GSM_URL_output.csv',
                        help='path to the URL of GSE and GSM ')
    parser.add_argument('-n', '--need_parser_GSE', type=str, required=False, 
                        help='a one-column file contains GSE numbers that are needed ')
    parser.add_argument('-el', '--exclude',  required=False,
                        help="a one-column file contains GSE number that has been parsed")
    parser.add_argument('-d', '--date_region', required=False,
                        help="crawl GSDs for a specific time range") ##'2023/01/02-2023/02/04'
    parser.add_argument('-p', '--geo_path', required=False,default='geo_gse',
                        help="path to xml files")
    parser.add_argument('-k', '--keywords_yaml', type=str, required=False, default='keywords.yml',
                        help='a yaml file contains key words for different data types')
    parser.add_argument('-t', '--ttype', required=True, default='sc-rna-seq',
                        help="which type of data you want to crawl, seperate by ','")
    parser.add_argument("-ef",'--error_file_csv',required=False,default='error_file_name.csv',
                        help='log for error files')
    parser.add_argument('-ol','--only_GSE',action="store_true",
                        help='just only save information we need')
    
    
    args = parser.parse_args()
    gse_output_csv = args.gse_output_csv
    gse_gsm_output_csv = args.gse_gsm_output_csv
    need_parser_GSE = args.need_parser_GSE
    exludeFile = args.exclude
    date_region = args.date_region
    geo_path = args.geo_path
    keywords_yaml = args.keywords_yaml
    ttype = args.ttype
    error_file_name_path = args.error_file_csv
    only_GSE = args.only_GSE

    ttype = ttype.split(',')

    keywords = yaml.load(open(keywords_yaml))
    # Write headers to CSV files
    with open(gse_output_csv, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        csvwriter.writerow(["GSE", "GSE_Link", "GSE_Title", "GSE_Pubmed_ID", "GSE_Last_Update_Date", "GSE_Release_Date", "GSE_URL",'Matched_key','Matched_key_word'])
    
    with open(error_file_name_path, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        csvwriter.writerow(['gse','gsm'])

    with open(gse_gsm_output_csv, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t')
        csvwriter.writerow(["GSE", "GSM", "GSM_Organism", "GSM_URL"])
    start = time.time()
    if need_parser_GSE:
        with open(need_parser_GSE, mode='r', newline='', encoding='utf-8') as file:
            df = pd.read_csv(file)
            gdsSamples = df.iloc[:,0].tolist()  ##GSE_id
        for gseid in gdsSamples:
            p_f.get_meta_url(gseid = gseid,
                             gse_output_csv=gse_output_csv,
                             gse_gsm_output_csv=gse_gsm_output_csv,
                             error_file_name_path = error_file_name_path,
                             keywords = keywords,
                             ttype = ttype,
                             path = geo_path,
                             only_GSE = only_GSE)
            print("GSE:{} is finished.".format(gseid))
            end = time.time()
            print('Total time:', end - start)
    else:
        if exludeFile:   # GSE_id
            with open(exludeFile, mode='r', newline='', encoding='utf-8') as file:
                df = pd.read_csv(file)
                exludeFile_list = df.iloc[:,0].tolist()
            gdsSamples = p_f.getGDSSamples(date_region)
            for gds in gdsSamples:
                gseid = p_f.gse_idToAcc(gds)
                if gseid not in exludeFile_list:
                    p_f.get_meta_url(gseid = gseid,
                                    gse_output_csv=gse_output_csv,
                                    gse_gsm_output_csv=gse_gsm_output_csv,
                                    error_file_name_path = error_file_name_path,
                                    keywords = keywords,
                                    ttype = ttype,
                                    path = geo_path,
                                    only_GSE = only_GSE)
                    end = time.time()
                    print('Total time:', end - start)
        else:
            gdsSamples = p_f.getGDSSamples(date_region)
            gdsSamples_len = len(gdsSamples)
            print("gdsSamples_len:{}".format(gdsSamples_len))
            for i in range(gdsSamples_len):
                with tqdm(total=gdsSamples_len) as pbar:
                    gseid = p_f.gse_idToAcc(gdsSamples[i])
                    pbar.set_description(f"Processing {gseid}")
                    # print('{} is starting'.format(gseid))
                    p_f.get_meta_url(gseid = gseid,
                                    gse_output_csv=gse_output_csv,
                                    gse_gsm_output_csv=gse_gsm_output_csv,
                                    error_file_name_path = error_file_name_path,
                                    keywords = keywords,
                                    ttype = ttype,
                                    path = geo_path,
                                    only_GSE = only_GSE)
                    pbar.update(i)
                    # end = time.time()
                    # print('Total time:', end - start)
                    
if __name__ == "__main__":
    main()
    

    
    