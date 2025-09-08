#!/bin/bash

root_folder='/fse/home/wangyongyan/Project/parser_code/result/2024_01_01_2025_08_17/parser'
save_path='/fse/home/wangyongyan/Project/parser_code/result/2024_01_01_2025_08_17/rectify'
word_bank='/fse/home/wangyongyan/Project/parser_code/final_parser_code/parser_second_part/word_bank'

python parser_rectify_meta.py --root_folder $root_folder --save_path $save_path --word_bank $word_bank