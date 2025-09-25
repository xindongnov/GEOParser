gse_output_csv='/fse/home/dongxin/parse_result_rri_20250901/gse_output.csv'
gse_gsm_output_csv='/fse/home/dongxin/parse_result_rri_20250901/gse_gsm_output.csv'
geo_path='/fse/home/dongxin/GEO'
error_file_csv='/fse/home/dongxin/parse_result_rri_20250901/error_file.csv'
ttype='RRI'
# need_parser_GSE='/fse/home/wangyongyan/Project/parser_code/result/2024.01_2025.04.09/Stereo/2024.01_2025.04.09_GSE_meta_output_v2.csv'

python parser_head_local.py --geo_path $geo_path --gse_output_csv $gse_output_csv \
                    --gse_gsm_output_csv $gse_gsm_output_csv  --error_file_csv $error_file_csv \
                    --keywords_yaml keywords.yml --ttype $ttype --only_GSE
