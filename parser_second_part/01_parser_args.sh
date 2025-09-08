save_path='/fse/home/wangyongyan/Project/parser_code/result/2024_01_01_2025_08_17/parser'
existed_GSE='/fse/home/wangyongyan/Project/parser_code/result/2024_01_01_2025_08_17/2024.01.01_2025.08.17_GSE_meta_output.csv'
root_folder='/fse/home/wangyongyan/Project/parser_code/GEO_XML_2024.01_2025.08.17/geo_gse'

python parser_args.py --save_path $save_path --existed_GSE $existed_GSE --root_folder $root_folder 

# python parser_args.py --save_path $save_path --root_folder $root_folder 