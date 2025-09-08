gse_output_csv='/fse/home/dongxin/parse_result_clip_20250901/gse_output.csv'
gse_gsm_output_csv='/fse/home/dongxin/parse_result_clip_20250901/gse_gsm_output.csv'
date_region='2005/01/01-2025/09/01'
geo_path='/fse/home/dongxin/parse_result_clip_20250901/error_file.csv'
ttype="clip"
error_file_csv='/fse/home/dongxin/parse_result_clip_20250901/error_file.csv'
# need_parser_GSE='/fse/home/wangyongyan/Project/parser_code/GEO_XML_2024.01_2025.04.09/2024.01_2025.04.09_GSE_meta_output_v2.csv'

python parser_head.py --geo_path $geo_path --ttype $ttype --gse_output_csv $gse_output_csv \
                    --gse_gsm_output_csv $gse_gsm_output_csv  --date_region $date_region \
                    --error_file_csv $error_file_csv --only_GSE
                    # --need_parser_GSE $need_parser_GSE
