import os
import re
import csv
import time
import xml.etree.ElementTree as ET
import pandas as pd
import subprocess
import argparse

def run_parser_script(script_name, root_folder, csv_file_name, no_meta_file_name, error_file_name,existed_GSE,save_path):
    try:
        if existed_GSE != '':
            command = [
                'nohup', 'python', f'parser_{script_name}.py',
                '--root_folder', root_folder,
                '--all_meta_file_path', csv_file_name,
                '--no_meta_file_path', no_meta_file_name,
                '--error_file_path', error_file_name,
                '--existed_GSE', existed_GSE,
                '--only',
                '--save_path',save_path
            ]
            command_str = ' '.join(command) + f' > {script_name}_nohup.log '
            subprocess.run(command_str, check=True,shell= True)
        else:
            command = [
                'nohup', 'python', f'parser_{script_name}.py',
                '--root_folder', root_folder,
                '--all_meta_file_path', csv_file_name,
                '--no_meta_file_path', no_meta_file_name,
                '--error_file_path', error_file_name,
                '--only',
                '--save_path',save_path
            ]
            command_str = ' '.join(command) + f' > {script_name}_nohup.log '
            subprocess.run(command_str, check=True,shell= True)
    except subprocess.CalledProcessError as e:
        print(f"Error while executing parser_{script_name}.py: {e}")

def run_parser_script_no_GSM(script_name, root_folder, csv_file_name, no_meta_file_name, error_file_name,existed_GSE,save_path):
    try:
        if existed_GSE != '':
            command = [
                'nohup', 'python', f'parser_{script_name}.py',
                '--root_folder', root_folder,
                '--all_meta_file_path', csv_file_name,
                '--no_meta_file_path', no_meta_file_name,
                '--error_file_path', error_file_name,
                '--existed_GSE', existed_GSE,
                '--save_path',save_path
            ]
            command_str = ' '.join(command) + f' > {script_name}_nohup.log '
            subprocess.run(command_str, check=True,shell= True)
        else:
            command = [
                'nohup', 'python', f'parser_{script_name}.py',
                '--root_folder', root_folder,
                '--all_meta_file_path', csv_file_name,
                '--no_meta_file_path', no_meta_file_name,
                '--error_file_path', error_file_name,
                '--save_path',save_path
            ]
            command_str = ' '.join(command) + f' > {script_name}_nohup.log '
            subprocess.run(command_str, check=True,shell= True)
    except subprocess.CalledProcessError as e:
        print(f"Error while executing parser_{script_name}.py: {e}")
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser()
    # data
    parser.add_argument('-a', '--all_meta_file_path', type=str, default='meta.csv', help='path to all_meta')
    parser.add_argument('-n', '--no_meta_file_path', type=str, default='no_meta_file.csv', help='path to no_meta_file')
    parser.add_argument('-e', '--error_file_path', type=str, default='error_file.csv', help='path to erroe_file')
    parser.add_argument('-p', '--need_parser_parameter', nargs='+', type=str, default=['all'], help="choose needed parameter to parser")
    parser.add_argument('-r', '--root_folder', type=str, help='path to dirs of GSE ')
    parser.add_argument('-d', '--existed_GSE', type=str, default='', help='path to existed GSE file ')
    parser.add_argument('-ol', '--only', action="store_true", help='just only save information we need')
    parser.add_argument('-sp', '--save_path', required=False, default='', help='the path to save result')

    args = parser.parse_args()

    csv_file_name = args.all_meta_file_path
    no_meta_file_name = args.no_meta_file_path
    error_file_name = args.error_file_path
    need_parser_parameter = args.need_parser_parameter
    root_folder = args.root_folder
    existed_GSE = args.existed_GSE
    only_need = args.only
    save_path = args.save_path
    
    os.makedirs(save_path, exist_ok=True)

    logging.debug(f"only_need: {only_need}")
    logging.debug(f"need_parser_parameter: {need_parser_parameter}")
    logging.debug(f"root_folder: {root_folder}")
    logging.debug(f"csv_file_name: {csv_file_name}")
    logging.debug(f"no_meta_file_name: {no_meta_file_name}")
    logging.debug(f"error_file_name: {error_file_name}")
    logging.debug(f"existed_GSE: {existed_GSE}")
    logging.debug(f"save_path: {save_path}")

    if only_need:
        logging.debug("only_need is True, running run_parser_script")
        if "all" in need_parser_parameter:
            scripts_to_run = ['age', 'race', 'platform', 'sex', 'tissue', 'disease', 'treatment']
            for script in scripts_to_run:
                run_parser_script(script, root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
        else:
            if "age" in need_parser_parameter:
                run_parser_script('age', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "race" in need_parser_parameter:
                run_parser_script('race', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "platform" in need_parser_parameter:
                run_parser_script('platform', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "sex" in need_parser_parameter:
                run_parser_script('sex', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "tissue" in need_parser_parameter:
                run_parser_script('tissue', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "disease" in need_parser_parameter:
                run_parser_script('disease', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "treatment" in need_parser_parameter:
                run_parser_script('treatment', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
    else:
        logging.debug("only_need is False, running run_parser_script_no_GSM")
        if "all" in need_parser_parameter:
            scripts_to_run = ['age', 'race', 'platform', 'sex', 'tissue', 'disease', 'treatment']
            for script in scripts_to_run:
                run_parser_script_no_GSM(script, root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
        else:
            if "age" in need_parser_parameter:
                run_parser_script_no_GSM('age', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "race" in need_parser_parameter:
                run_parser_script_no_GSM('race', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "platform" in need_parser_parameter:
                run_parser_script_no_GSM('platform', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "sex" in need_parser_parameter:
                run_parser_script_no_GSM('sex', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "tissue" in need_parser_parameter:
                run_parser_script_no_GSM('tissue', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "disease" in need_parser_parameter:
                run_parser_script_no_GSM('disease', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)
            if "treatment" in need_parser_parameter:
                run_parser_script_no_GSM('treatment', root_folder, csv_file_name, no_meta_file_name, error_file_name, existed_GSE, save_path)

if __name__ == "__main__":
    main()