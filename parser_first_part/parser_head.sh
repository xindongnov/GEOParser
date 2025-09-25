gse_output_csv='/home/xdong/Projects/parser/parse_result_rri_20250901/gse_output.tsv'
gse_gsm_output_csv='/home/xdong/Projects/parser/parse_result_rri_20250901/gse_gsm_output.tsv'
date_region='2005/01/01-2025/09/01'
geo_path='/home/xdong/Projects/parser/GEO'
ttype="RRI"
error_file_csv='/home/xdong/Projects/parser/parse_result_rri_20250901/error_file.tsv'
# need_parser_GSE='/fse/home/wangyongyan/Project/parser_code/GEO_XML_2024.01_2025.04.09/2024.01_2025.04.09_GSE_meta_output_v2.csv'

python parser_head.py --geo_path $geo_path --ttype $ttype --gse_output_csv $gse_output_csv \
                    --gse_gsm_output_csv $gse_gsm_output_csv  --date_region $date_region \
                    --keywords_yaml keywords.yml --ttype $ttype --only_GSE
                    # --only_GSE
                    # --need_parser_GSE $need_parser_GSE
