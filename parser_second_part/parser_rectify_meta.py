import os
import argparse
import pandas as pd

from rectify_funtion import *



def main():
    parser = argparse.ArgumentParser()
    # data
    parser.add_argument('-r', '--root_folder', type=str,help='path to dirs of meta imformation')
    parser.add_argument('-sp','--save_path',required=False,default='',help='the path to save result')
    parser.add_argument('-w', '--word_bank', type=str,help='path to dirs of word_bank')
    
    args = parser.parse_args()
    
    root_folder = args.root_folder
    save_path = args.save_path
    word_bank_path = args.word_bank
    store_path = os.path.dirname(save_path)
    os.makedirs(save_path, exist_ok=True)
    print(f'{save_path} is created...')
    
    ## path
    disease_meta_path = os.path.join(root_folder,'disease/disease_meta.csv')
    age_meta_path = os.path.join(root_folder,'age/age_meta.csv')
    
    sex_meta_path = os.path.join(root_folder,'sex/sex_meta.csv')
    race_meta_path = os.path.join(root_folder,'race/race_meta.csv')
    platform_meta_path = os.path.join(root_folder,'platform/platform_meta.csv')
    tissue_meta_path = os.path.join(root_folder,'tissue/tissue_meta.csv')
    treatment_meta_path = os.path.join(root_folder,'treatment/treatment_meta.csv')
    
    ## save path
    ## age
    os.makedirs(os.path.join(save_path,'age'), exist_ok=True)
    age_gse_gsm_out_csv = os.path.join(save_path,'age/age_gse_gsm_out.csv')
    age_gse_out_csv = os.path.join(save_path,'age/age_gse_out.csv')
    age_mouse_right_csv = os.path.join(save_path,'age/mouse_accuracy.csv')
    age_human_right_csv = os.path.join(save_path, 'age/human_accuracy.csv')
    
    print('Human_age_correct is starting...')
    human_age_correct(
    age_meta_path=age_meta_path,
    output_path=os.path.join(save_path,'age'))
    print('Human_age_correct is finished...')
    
    print('Mouse_age_correct is starting...')
    mouse_age_correct(
    age_meta_path=age_meta_path,
    output_path=os.path.join(save_path,'age'))
    print('Mouse_age_correct is finished...')
    
    
    print('rectify_age_meta is starting...')
    age_gse_gsm_df, age_gse_df = rectify_age_meta(
    mouse_right_csv=age_mouse_right_csv,
    human_right_csv=age_human_right_csv,
    gse_gsm_out_csv=age_gse_gsm_out_csv,
    gse_out_csv=age_gse_out_csv
)
    print('rectify_age_meta is finshed...')
    ## sex
    os.makedirs(os.path.join(save_path,'sex'), exist_ok=True)
    corr_name = os.path.join(save_path,'sex/sex_corresponding_df.csv')
    sex_gse_gsm_out_csv = os.path.join(save_path,'sex/sex_gse_gsm_out.csv')
    sex_gse_out_csv = os.path.join(save_path,'sex/sex_gse_out.csv')

    print('rectify_sex_meta is starting...')
    sex_gse_gsm_df, sex_gse_df = rectify_sex_meta(
    in_csv=sex_meta_path,
    corr_name=corr_name,
    sex_gse_gsm_out_csv=sex_gse_gsm_out_csv,
    sex_gse_out_csv = sex_gse_out_csv
    )
    print('rectify_sex_meta is finshed...')
    
    ## platform
    os.makedirs(os.path.join(save_path,'paltform'), exist_ok=True)
    platform_level_path = os.path.join(word_bank_path,'platform_mapping_list.tsv')
    gse_gsm_platform_output_csv =os.path.join(save_path,'paltform/platform_gse_gsm_out.csv')
    gse_platform_output_csv = os.path.join(save_path,'paltform/platform_gse_out.csv')
    
    print('rectify_platform_meta is starting...')
    platform_gse_gsm_df, platform_gse_df = rectify_platform_meta(gse_gsm_platform_path = platform_meta_path,
                            platform_level_path = platform_level_path,
                            gse_gsm_platform_output_csv = gse_gsm_platform_output_csv,
                            gse_platform_output_csv =   gse_platform_output_csv)
    print('rectify_platform_meta is finshed...')
    

 

    
    ## treatment
    os.makedirs(os.path.join(save_path,'treatment'), exist_ok=True)
    gse_gsm_treatment_output_csv = os.path.join(save_path,'treatment/gse_gsm_treatment_out.csv')
    gse_treatment_output_csv = os.path.join(save_path,'treatment/gse_treatment_out.csv')
    
    drug_list_path = os.path.join(word_bank_path,'drug.txt')
    gene_key_path = os.path.join(word_bank_path,'gene.txt')
    control_path = os.path.join(word_bank_path,'ctrl.txt')
    cyto_key_path = os.path.join(word_bank_path,'cytokine.xlsx')

    
    print('rectify_treatment_meta is starting...')
    treat_gse_gsm_df, treat_gse_df = rectify_treatment_meta(gse_gsm_treat_path= treatment_meta_path,
                            drug_list_path = drug_list_path,
                            gene_key_path = gene_key_path,
                            control_path =  control_path,
                            cyto_key_path = cyto_key_path,
                            gse_gsm_treatment_output_csv = gse_gsm_treatment_output_csv,
                            gse_treatment_output_csv =  gse_treatment_output_csv
                            )
        
    print('rectify_treatment_meta is finshed...')
    
    
    ## disease 
    os.makedirs(os.path.join(save_path,'disease'), exist_ok=True
                )
    disease_type_path = os.path.join(word_bank_path,'DiseaseType.csv')
    disease_name_type_path = os.path.join(word_bank_path,'DiseaseType_dictionary_1-multi.tsv')
    disease_gse_gsm_out_csv = os.path.join(save_path,'disease/disease_gse_gsm_out.csv')
    disease_gse_out_csv = os.path.join(save_path,'disease/disease_gse_out.csv')

    
    print('rectify_disease_meta is starting...')
    disease_gse_gsm_df, disease_gse_df = rectify_disease_meta(gse_gsm_disease_path = disease_meta_path,
                            disease_type_path = disease_type_path,
                            disease_name_type_path = disease_name_type_path,
                            gse_gsm_disease_output_csv =   disease_gse_gsm_out_csv,
                            gse_disease_output_csv = disease_gse_out_csv)
    print('rectify_disease_meta is finshed...')
    
    ## tissue
    os.makedirs(os.path.join(save_path,'tissue'), exist_ok=True)
    tissue_level_path = os.path.join(word_bank_path,'Simplify_new.txt')
    cancer_cellline_path = os.path.join(word_bank_path,'cancer_cellline.txt')
    tissue_gse_gsm_out_csv=os.path.join(save_path,'tissue/gse_gsm_tissue_meta_rectify_test.csv')
    tissue_gse_out_csv= os.path.join(save_path,'tissue/gse_tissue_meta_rectify_test.csv')

    print('rectify_tissue_meta is starting...')
    tissue_gse_gsm_df, tissue_gse_df = rectify_tissue_meta(gse_gsm_tissue_path = tissue_meta_path,
                            tissue_level_path = tissue_level_path,
                            cancer_cellline_path = cancer_cellline_path,
                            gse_gsm_tissue_output_csv = tissue_gse_gsm_out_csv,
                            gse_tissue_output_csv = tissue_gse_out_csv
                            )

    print('rectify_tissue_meta is finshed...')
  
    
        
    ## race 
    print('rectify_race_meta is starting...')
    os.makedirs(os.path.join(save_path,'race'), exist_ok=True)
    race_gse_gsm_df = pd.read_csv(race_meta_path, sep=',', index_col=False)
    race_gse_df = race_gse_gsm_df[race_gse_gsm_df['gsm'].str.contains('GSE')]
    race_gse_gsm_df.to_csv(os.path.join(save_path,'race/race_gse_gsm_df.csv'),index=False)
    race_gse_df.to_csv(os.path.join(save_path,'race/race_gse_df.csv'),index=False)
    print('rectify_race_meta is finshed...')
    
    ##
    sex_cols = ['gse', 'gsm', 'sex_final']
    platform_cols = ["gse","gsm","platform_final"]
    treat_cols = ['gse', 'gsm', 'treat_type', 'treat_final'] 
    disease_cols = ['gse', 'gsm', "Disease_Global_type","Disease_Anatomical_type"]
    tissue_cols = ['gse', 'gsm', 'tissue_final', 'tissue_type']
    race_cols = ['gse', 'gsm', 'char_race']
    
    ## merge
    gse_gsm_list = [sex_gse_gsm_df, platform_gse_gsm_df, treat_gse_gsm_df, disease_gse_gsm_df, tissue_gse_gsm_df,race_gse_df]
    gse_list = [ sex_gse_df, platform_gse_df, treat_gse_df, disease_gse_df, tissue_gse_df,race_gse_df]
    need_cols = [sex_cols, platform_cols, treat_cols, disease_cols, tissue_cols, race_cols]
    gse_gsm_merged_df = age_gse_gsm_df.copy()
    gse_merged_df = age_gse_df.copy()
    
    print('Merge gse_gsm_list is starting...')
    for i in range(len(gse_gsm_list)):
        gse_gsm_merged_df = pd.merge(gse_gsm_merged_df, gse_gsm_list[i][need_cols[i]], on=['gse', 'gsm'], how='left')
        
    print('Merge gse_list is starting...')
    for i in range(len(gse_list)):
        gse_merged_df = pd.merge(gse_merged_df, gse_list[i][need_cols[i]], on=['gse', 'gsm'], how='left')
    
    gse_gsm_merged_df.to_csv(os.path.join(store_path,'gse_gsm_meta_rectify.csv'), index=False)
    gse_merged_df.to_csv(os.path.join(store_path,'gse_meta_rectify.csv'), index=False)
    
    return gse_gsm_merged_df, gse_merged_df

if __name__ == '__main__':
    main()
    
