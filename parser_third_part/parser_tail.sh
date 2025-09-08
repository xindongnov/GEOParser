#!/bin/bash
url_csv='/fse/home/wangyongyan/Project/parser_code/final_parser_code/parser_third_part/need_down_GSE_URL_df.csv'
folder_path='/fse/DC/human/third_unprocessed'
name_of_URL_column='GSE_URL'
name_of_GSE_column='GSE'

python parser_tail.py --url_csv $url_csv --folder_path $folder_path --name_of_URL_column $name_of_URL_column --name_of_GSE_column $name_of_GSE_column