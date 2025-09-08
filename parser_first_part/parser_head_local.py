import os
import csv
import time
import argparse
import parser_function as p_f

def main():
    parser = argparse.ArgumentParser()
    # data
    parser.add_argument('-o', '--gse_output_csv', type=str, default='GSE_meta_output.csv',
                        help='path to some parts of meta')
    parser.add_argument('-g', '--gse_gsm_output_csv', type=str, default='GSE_GSM_URL_output.csv',
                        help='path to the URL of GSE and GSM ')
    parser.add_argument('-p', '--geo_path', required=False, default='geo_gse',
                        help="path to xml files")
    parser.add_argument('-t', '--ttype', required=False, default='',
                        help="which type of data you want to crawl")
    parser.add_argument("-ef",'--error_file_csv',required=False, default='error_file_name.csv',
                        help='log for error files')
    parser.add_argument('-ol','--only_GSE',action="store_true",
                        help='just only save information we need')

    args = parser.parse_args()
    gse_output_csv = args.gse_output_csv
    gse_gsm_output_csv = args.gse_gsm_output_csv
    geo_path = args.geo_path
    ttype = args.ttype
    error_file_name_path = args.error_file_csv
    only_GSE = args.only_GSE
    print(only_GSE)

    dirs =os.listdir(geo_path)
    for folder in dirs:
        if folder.startswith('GSE'):  ### 
            try:
                    # if folder in GSEIDs: ## if want to get the specific GSEIDs(such as ST's GSE, use this)
                folder_path = os.path.join(geo_path, folder) 
                gse_xml_path = os.path.join(folder_path, "%s.xml" % folder) 
                p_f.get_meta_url(folder,gse_output_csv,gse_gsm_output_csv,path = geo_path,error_file_name_path = error_file_name_path,ttype=ttype,only_GSE = only_GSE)
                print("GSE:{} is finished.".format(folder))
            except:
                print("GSE:{} is failed.".format(folder))
                break
        

if __name__ == "__main__":
    main()




    

    
    