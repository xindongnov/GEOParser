import os
import re
import csv
import time
import xml.etree.ElementTree as ET
import pandas as pd
import argparse
from itertools import islice
def read_xmlfile(xml_path):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        return root
    except ET.ParseError as e:
        print(f"[WARN] 跳过损坏XML {xml_path}: {e}")
        return None
# 提取所有不同的元素标签
#unique_tags = set()
#def extract_unique_tags(element):
#    unique_tags.add(element.tag)
#    for child in element:
#        extract_unique_tags(child)
#extract_unique_tags(xml_root)

##GSE
GSE_Element_names = ["Summary","Overall-Design",
                     "Title","Supplementary-Data",
                     "Last-Update-Date","Submission-Date",
                     "Release-Date","Pubmed-ID"]

def get_element_data(xml_root, element_name, attribute_name):
    elements = xml_root.findall(f".//{{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}}{element_name}")
    element_dict = {}
    element_list = []
    for element in elements:
        attribute_value = element.text.strip() if element.text else None
        if attribute_value:
            element_list.append(attribute_value)
    if element_list:
        element_dict[attribute_name] = element_list
    return element_dict

def get_GSE_dict(xml_root):
    GSE_dict = {}
    GSE_dict.update(get_element_data(xml_root, 'Summary', 'GSE_Summary'))
    GSE_dict.update(get_element_data(xml_root, 'Overall-Design', 'GSE_Overall_Design'))
    GSE_dict.update(get_element_data(xml_root, 'Title', 'GSE_Title'))
    GSE_dict.update(get_element_data(xml_root, 'Supplementary-Data', 'GSE_URL'))
    GSE_dict.update(get_element_data(xml_root, 'Last-Update-Date', 'GSE_Last_Update_Date'))
    GSE_dict.update(get_element_data(xml_root, 'Submission-Date', 'GSE_Submission_Date'))
    GSE_dict.update(get_element_data(xml_root, 'Release-Date', 'GSE_Release_Date'))
    GSE_dict.update(get_element_data(xml_root, 'Pubmed-ID', 'GSE_Pubmed_ID'))
    return GSE_dict


##GSM
GSM_Element_names = ["Title","Accession","Characteristics",
                     "Organism","Last-Update-Date","Source",
                     "Release-Date","Submission-Date",
                     "Supplementary-Data","Treatment-Protocol",'Extract-Protocol']

def get_GSM_dict(xml_root):
    GSM_dict = {}
    GSM_dict.update(get_element_data(xml_root, 'Title', 'GSM_Title'))
    GSM_dict.update(get_element_data(xml_root, 'Accession', 'GSM_Accession'))
    GSM_dict.update(get_element_data(xml_root, 'Characteristics', 'GSM_Characteristics'))
    GSM_dict.update(get_element_data(xml_root, 'Organism', 'GSM_Organism'))
    GSM_dict.update(get_element_data(xml_root, 'Supplementary-Data', 'GSM_GSM_URL'))
    GSM_dict.update(get_element_data(xml_root, 'Last-Update-Date', 'GSM_Last_Update_Date'))
    GSM_dict.update(get_element_data(xml_root, 'Submission-Date', 'GSM_Submission_Date'))
    GSM_dict.update(get_element_data(xml_root, 'Source', 'GSM_Source'))
    GSM_dict.update(get_element_data(xml_root, 'Release-Date', 'GSM_Release_Date'))
    GSM_dict.update(get_element_data(xml_root, 'Treatment-Protocol', 'GSM_Treatment_Protocol'))
    GSM_dict.update(get_element_data(xml_root, 'Extract-Protocol', 'GSM_Extract-Protocol'))
    return GSM_dict

##characteristics
def get_characteristics_dict(xml_root):
    characteristics = xml_root.findall(".//{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics")
    characteristics_dict = {}
# 遍历 characteristics 元素列表
    for characteristic in characteristics:
    # 获取 Characteristic 的属性名和属性值
        attribute_name = characteristic.attrib.get('tag',"Unkown")
        attribute_value = characteristic.text.strip()
    # 将属性名和属性值添加到字典中
        if attribute_name and attribute_value is not None:
            characteristics_dict[attribute_name.lower()] = attribute_value
        else:
            characteristics_dict[attribute_name.lower()] = None
    return characteristics_dict

def find_characteristics_key(dict,key_list):
    characteristics_dict = dict
    keys_to_check = key_list
    characteristics_need={}
    for key in keys_to_check:
        for characteristics in characteristics_dict.keys():
            if key in characteristics:
                characteristics_need[characteristics] = characteristics_dict[characteristics]
    if characteristics_need is not None:
        return characteristics_need
    else:
        return None
    

##set
def generate_ordered_ngrams(text, n):
    # 将文本按空格分割成单词列表
    words = text.lower().split()
    # 使用itertools.islice和zip生成n-gram
    n_grams = [' '.join(w) for w in zip(*[islice(words, i, None) for i in range(n)])]
    return n_grams

def make_key_value(species_bank):
    species_dict ={}
    #使用zip函数将两列组合成一个元组列表
    #然后使用dict函数将元组列表转换为字典
    for i in range(species_bank.shape[0]):
        second_column = species_bank.iloc[i, 1]
        third_column = species_bank.iloc[i, 2]  # 第三列
        species_dict[second_column] = third_column
    return species_dict

##set
def find_key_words(dict,word_bank_name_list,word_bank_list):
    match_key_word = []
    match_key_class = []
    matched_key_word = []
    word_bank_name_list = word_bank_name_list
    word_bank_list = word_bank_list
    i = 0
    for word_bank in word_bank_list:
        word_bank_name = word_bank_name_list[i]
        longest_word = len(max(word_bank, key=len).split())
        word_bank_set = set(word.lower() for word in set(word_bank))
        for key in dict.keys():
            text_all_set = []
            text = dict[key]
            if key in ['GSM_Treatment_Protocol','GSM_Extract-Protocol']:
                for j in range(longest_word+1):
                    text_set_j = set(generate_ordered_ngrams(text[0],j))
                    text_all_set.append(text_set_j)
                for text_set_j in text_all_set:
                    intersection_set = word_bank_set & text_set_j
                    for words in intersection_set:
                        match_key_word.append(words)
                        match_key_class.append(word_bank_name)
                        matched_key_word.append(key)
            else:
                text_list = [word.lower() for phrase in text for word in phrase.split()]
                for key_text in text_list:
                    if key_text in word_bank_set:
                        match_key_word.append(key_text)
                        match_key_class.append(word_bank_name)
                        matched_key_word.append(key)  
        i = i+1  
    return matched_key_word,match_key_class,match_key_word

def find_key_by_value(dict_of_elements, search_value):
    for key, value in dict_of_elements.items():
        if value == search_value:
            return key
    return None  # 如果没有找到匹配的值
def find_ko_gene_symbol(dict,word_bank,match_key_meta,species):
    match_key_word = []
    word_bank_key_set = set(word_bank.keys())
    text = dict[match_key_meta]
    text_list =[word for phrase in text for word in phrase.split()]
    text_list_cleaned = [re.sub(r'[^a-zA-Z0-9\s]', '', word) for word in text_list]
    text_list_cleaned = [word for word in text_list_cleaned if word]
    for words in text_list:
        if words in word_bank_key_set:
            match_key_word.append(words)
    if len(match_key_word)==0 and species == "['Mus musculus']" :
        for words in text_list_cleaned:
            pattern = rf'\b{words}\b'
            match_list = [value for value in word_bank.values() if re.search(pattern, value, re.IGNORECASE)]
            for match_value in match_list:
                if len(match_value) == len(words):
                    key_words = find_key_by_value(word_bank, match_value)
                    match_key_word.append(key_words)
    return match_key_word
            
            
def find_all_meta(folder_path,file,key_list,word_bank_name_list,word_bank_list):
    if file.startswith('GSM') and file.endswith('.xml'):
        gsm = file.replace('.xml', '')
        xml_file_path = os.path.join(folder_path, file)
        gsm_xml_root = read_xmlfile(xml_file_path)
        GSM_meta_dict=get_GSM_dict(gsm_xml_root)
        characteristics_dict = get_characteristics_dict(gsm_xml_root)
        char_key = ','.join(list(find_characteristics_key(characteristics_dict,key_list).keys()))
        result_dict = find_characteristics_key(characteristics_dict,key_list).values()
        char = ','.join(list(str(value) for value in result_dict  if value is not None))
                #GSE_GSM_dict ={}
                #GSE_GSM_dict.update(GSE_meta_dict)
                #GSE_GSM_dict.update(GSM_meta_dict)
        match_key,key_bank,key_word=find_key_words(GSM_meta_dict,word_bank_name_list,word_bank_list)
        return gsm,char_key,char,match_key,key_bank,key_word,GSM_meta_dict
    



def main():
    parser = argparse.ArgumentParser()
    # data
    parser.add_argument('-a', '--all_meta_file_path', type=str, default='meta.csv',help='path to all_meta')
    parser.add_argument('-n', '--no_meta_file_path', type=str, default='no_meta_file.csv',help='path to no_meta_file')
    parser.add_argument('-e', '--error_file_path', type=str, default='error_file.csv',help='path to erroe_file')
    parser.add_argument('-r', '--root_folder', type=str,default='/fse/DC/sc_geo_gse/',help='path to dirs of GSE ')
    parser.add_argument('-d', '--existed_GSE', type=str,default = '',help='path to existed GSE file ')
    parser.add_argument('-ol','--only',action="store_true",help='just only save information we need')
    parser.add_argument('-sp','--save_path',required=False,default='',help='the path to save result')

    args = parser.parse_args()
    
    csv_file_name = 'treatment_' + args.all_meta_file_path
    no_meta_file_name = 'treatment_' +  args.no_meta_file_path
    error_file_name  = 'treatment_' + args.error_file_path
    root_folder = args.root_folder
    existed_GSE = args.existed_GSE
    only_need = args.only
    save_path = args.save_path
##char_key_list
    Treatment_key_list = ['treat']
    ##treatment
    with open('./word_bank/drug.txt','r',encoding='ISO-8859-1')as f:
        content = f.read()
    drug_bank=content.strip().split('\n')
    drug_bank=[drug.lower()for drug in drug_bank]
    drug_bank.remove('')#drug匹配不区分大小写
    df = pd.read_excel("./word_bank/cytokine.xlsx", engine='openpyxl') 
    cyto_bank = df.iloc[:, 0].tolist()
    with open('./word_bank/gene.txt', 'r',encoding='ISO-8859-1') as f:
        content = f.read()
    gene_bank = content.strip().split('\n')
    with open('./word_bank/ctrl.txt','r',encoding='ISO-8859-1')as f:
        content = f.read()
    ctrl_bank=content.strip().split('\n')
    human_gene_bank=pd.read_csv('./word_bank/human_data_protein_coding.tsv',sep='\t')
    #human_gene_bank = human_gene_bank.iloc[:, 1].tolist()
    mouse_gene_bank=pd.read_csv('./word_bank/mouse_data_protein_coding.tsv',sep='\t')
    #mouse_gene_bank = mouse_gene_bank.iloc[:, 1].tolist()
    # species_bank = pd.read_csv('/fse/home/wangyongyan/new_meta_parser_7.0/age/age_meta.csv',sep=',')
    # species_dict = make_key_value(species_bank)
    species_dict = []
    human_gene_dict = make_key_value(human_gene_bank)
    keys = list(human_gene_dict.keys())
    human_gene_dict.pop(keys[-1])
    mouse_gene_dict = make_key_value(mouse_gene_bank)
    existed_GSE_choose = False
    if existed_GSE != '':
        df = pd.read_csv(existed_GSE)
        gse_list = df.iloc[:, 0].tolist()
        existed_GSE_choose = True
    if save_path == '':
        if not os.path.exists('treatment'):
            os.makedirs("treatment")
            csv_file_name_path = os.path.join("treatment", csv_file_name)
            no_meta_file_name_path = os.path.join("treatment", no_meta_file_name)
            error_file_name_path = os.path.join("treatment", error_file_name)
    else:
        if not os.path.exists(os.path.join(save_path,"treatment")):
            os.makedirs(os.path.join(save_path,"treatment"))
        csv_file_name_path = os.path.join(save_path,"treatment", csv_file_name)
        no_meta_file_name_path =os.path.join(save_path,"treatment",  no_meta_file_name)
        error_file_name_path = os.path.join(save_path,"treatment",  error_file_name)
    if not os.path.isfile(csv_file_name_path):
        with open(csv_file_name_path, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                      "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                      'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                                    ] 
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader() 
    dirs =os.listdir(root_folder)
    for folder in dirs:
        start = time.time()
        counts = 1
        if folder.startswith('GSE'):  ### 
                    # if folder in GSEIDs: ## if want to get the specific GSEIDs(such as ST's GSE, use this)
            folder_path = os.path.join(root_folder, folder) 
            gse = folder
            files = os.listdir(folder_path)
            if existed_GSE_choose :
                if  gse in gse_list :
                    GSE_file = [file for file in files if file.startswith('GSE')]
                    xml_file_path = os.path.join(folder_path, GSE_file[0])
                    gse_xml_root= read_xmlfile(xml_file_path)
                    if gse_xml_root is None:
                        continue
                    GSE_meta_dict=get_GSE_dict(gse_xml_root)
                    GSE_drug_match_key,GSE_drug_key_bank,GSE_drug_key_word=find_key_words(GSE_meta_dict,["drug_bank","ctrl_bank"],[drug_bank,ctrl_bank])
                    GSE_gene_match_key,GSE_gene_key_bank,GSE_gene_key_word=find_key_words(GSE_meta_dict,["cyto_bank","gene_bank"],[cyto_bank,gene_bank])
                    with open(csv_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                    fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                    "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                    'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                                ]                            
                                    writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                    writer.writerow({'gse': gse, 'gsm': gse, 
                                                                                    'char_treatment_key':'','char_treatment':'',
                                                                                    'drug_ctrl_treatment_match_key':GSE_drug_match_key,'drug_ctrl_treatment_key_bank':GSE_drug_key_bank,'drug_ctrl_treatment_key_word':GSE_drug_key_word,
                                                                                    'gene_cyto_treatment_match_key':GSE_gene_match_key,'gene_cyto_treatment_key_bank':GSE_gene_key_bank,'gene_cyto_treatment_key_word':GSE_gene_key_word,
                                                                                    'gene_symbol':''
                                                                                    })
                    print("{} is starting ".format(gse))
                    if not only_need and len(files) > 1:
                        for file in files:
                            try:
                                gene_symbol = []
                                if file.startswith('GSM') and file.endswith('.xml'):
                                    gsm,char_treatment_key,char_treatment,drug_ctrl_treatment_match_key,drug_ctrl_treatment_key_bank,drug_ctrl_treatment_key_word,GSM_meta_dict = find_all_meta(folder_path,file,Treatment_key_list,["drug_bank","ctrl_bank"],[drug_bank,ctrl_bank])
                                    _,_,_,gene_cyto_treatment_match_key,gene_cyto_treatment_key_bank,gene_cyto_treatment_key_word,_ = find_all_meta(folder_path,file,Treatment_key_list,["cyto_bank","gene_bank"],[cyto_bank,gene_bank])
                                    for words_gene in [gene.lower() for gene in gene_bank]:
                                        if words_gene in gene_cyto_treatment_key_word:
                                            print(words_gene)
                                            position = gene_cyto_treatment_key_word.index(words_gene)
                                            match_key_meta = gene_cyto_treatment_match_key[position]
                                            species = species_dict[gsm]
                                            if species == "['Mus musculus']":
                                                gene_symbol.extend(find_ko_gene_symbol(GSM_meta_dict,mouse_gene_dict,match_key_meta,species))
                                            else:
                                                gene_symbol.extend(find_ko_gene_symbol(GSM_meta_dict,human_gene_dict,match_key_meta,species))
                                    #if char_treatment_key != '' or char_treatment != '' or drug_ctrl_treatment_match_key != [] or drug_ctrl_treatment_key_bank != [] or drug_ctrl_treatment_key_word != [] or gene_cyto_treatment_match_key != [] or gene_cyto_treatment_key_bank != [] or gene_cyto_treatment_key_word != []:
                                    with open(csv_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                            fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                        "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                        'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                                        ]                            
                                            writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                            writer.writerow({'gse': gse, 'gsm': gsm, 
                                                                                            'char_treatment_key':char_treatment_key,'char_treatment':char_treatment,
                                                                                            'drug_ctrl_treatment_match_key':drug_ctrl_treatment_match_key,'drug_ctrl_treatment_key_bank':drug_ctrl_treatment_key_bank,'drug_ctrl_treatment_key_word':drug_ctrl_treatment_key_word,
                                                                                            'gene_cyto_treatment_match_key':gene_cyto_treatment_match_key,'gene_cyto_treatment_key_bank':gene_cyto_treatment_key_bank,'gene_cyto_treatment_key_word':gene_cyto_treatment_key_word,
                                                                                            'gene_symbol':gene_symbol})
                            except Exception as e:
                                print(f"Error parsing {gse}: {e}")
                                with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                                fieldnames = ['gse','gsm']
                                                writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                                writer.writerow({'gse': gse,'gsm':gsm})
                                with open(csv_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                                fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                                "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                                'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                                            ]                            
                                                writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                                writer.writerow({'gse': gse, 'gsm': gsm, 
                                                                                                'char_treatment_key':'','char_treatment':'',
                                                                                                'drug_ctrl_treatment_match_key':'','drug_ctrl_treatment_key_bank':'','drug_ctrl_treatment_key_word':'',
                                                                                                'gene_cyto_treatment_match_key':'','gene_cyto_treatment_key_bank':'','gene_cyto_treatment_key_word':'',
                                                                                                'gene_symbol':''})
            else:
                GSE_file = [file for file in files if file.startswith('GSE')]
                xml_file_path = os.path.join(folder_path, GSE_file[0])
                gse_xml_root= read_xmlfile(xml_file_path)
                if gse_xml_root is None:
                    continue
                GSE_meta_dict=get_GSE_dict(gse_xml_root)
                print("{} is starting ".format(gse))
                GSE_drug_match_key,GSE_drug_key_bank,GSE_drug_key_word=find_key_words(GSE_meta_dict,["drug_bank","ctrl_bank"],[drug_bank,ctrl_bank])
                GSE_gene_match_key,GSE_gene_key_bank,GSE_gene_key_word=find_key_words(GSE_meta_dict,["cyto_bank","gene_bank"],[cyto_bank,gene_bank])
                with open(csv_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                            ]                            
                                writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                writer.writerow({'gse': gse, 'gsm': gse, 
                                                                                'char_treatment_key':'','char_treatment':'',
                                                                                'drug_ctrl_treatment_match_key':GSE_drug_match_key,'drug_ctrl_treatment_key_bank':GSE_drug_key_bank,'drug_ctrl_treatment_key_word':GSE_drug_key_word,
                                                                                'gene_cyto_treatment_match_key':GSE_gene_match_key,'gene_cyto_treatment_key_bank':GSE_gene_key_bank,'gene_cyto_treatment_key_word':GSE_gene_key_word,
                                                                                'gene_symbol':''})
                if not only_need and len(files) > 1:
                    for file in files:
                        try:
                            gene_symbol = []
                            if file.startswith('GSM') and file.endswith('.xml'):
                                print(file)
                                gsm,char_treatment_key,char_treatment,drug_ctrl_treatment_match_key,drug_ctrl_treatment_key_bank,drug_ctrl_treatment_key_word,_ = find_all_meta(folder_path,file,Treatment_key_list,["drug_bank","ctrl_bank"],[drug_bank,ctrl_bank])
                                _,_,_,gene_cyto_treatment_match_key,gene_cyto_treatment_key_bank,gene_cyto_treatment_key_word,GSM_meta_dict = find_all_meta(folder_path,file,Treatment_key_list,["cyto_bank","gene_bank"],[cyto_bank,gene_bank])
                                for words_gene in [gene.lower() for gene in gene_bank]:
                                    if words_gene in gene_cyto_treatment_key_word:
                                        position = gene_cyto_treatment_key_word.index(words_gene)
                                        match_key_meta = gene_cyto_treatment_match_key[position]
                                        species = species_dict[gsm]
                                        if species == "['Mus musculus']":
                                            gene_symbol.extend(find_ko_gene_symbol(GSM_meta_dict,mouse_gene_dict,match_key_meta,species))
                                        else:
                                            gene_symbol.extend(find_ko_gene_symbol(GSM_meta_dict,human_gene_dict,match_key_meta,species))
                                with open(csv_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                        fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                        "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                        'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                                    ]                            
                                        writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                        writer.writerow({'gse': gse, 'gsm': gsm, 
                                                                                        'char_treatment_key':char_treatment_key,'char_treatment':char_treatment,
                                                                                        'drug_ctrl_treatment_match_key':drug_ctrl_treatment_match_key,'drug_ctrl_treatment_key_bank':drug_ctrl_treatment_key_bank,'drug_ctrl_treatment_key_word':drug_ctrl_treatment_key_word,
                                                                                        'gene_cyto_treatment_match_key':gene_cyto_treatment_match_key,'gene_cyto_treatment_key_bank':gene_cyto_treatment_key_bank,'gene_cyto_treatment_key_word':gene_cyto_treatment_key_word,
                                                                                        'gene_symbol':gene_symbol})
                        except Exception as e:
                            print(f"Error parsing {gsm}: {e}")
                            with open(error_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                            fieldnames = ['gse','gsm']
                                            writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                            writer.writerow({'gse': gse,'gsm':gsm})
                            with open(csv_file_name_path, 'a', newline='', encoding='utf-8') as csvfile:
                                            fieldnames = ['gse', 'gsm', 'char_treatment_key','char_treatment',
                                                            "drug_ctrl_treatment_match_key",'drug_ctrl_treatment_key_bank','drug_ctrl_treatment_key_word',
                                                            'gene_cyto_treatment_match_key','gene_cyto_treatment_key_bank','gene_cyto_treatment_key_word','gene_symbol'
                                        ]                            
                                            writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
                                            writer.writerow({'gse': gse, 'gsm': gsm, 
                                                                                            'char_treatment_key':'','char_treatment':'',
                                                                                            'drug_ctrl_treatment_match_key':'','drug_ctrl_treatment_key_bank':'','drug_ctrl_treatment_key_word':'',
                                                                                            'gene_cyto_treatment_match_key':'','gene_cyto_treatment_key_bank':'','gene_cyto_treatment_key_word':'',
                                                                                            'gene_symbol':''})
            print("{} is finishing ".format(gse))
        end = time.time()
        print('Total time:', end - start)
        

    
if __name__ == "__main__":
    main()