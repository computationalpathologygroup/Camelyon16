# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 14:13:11 2017

@author: Babak Ehteshami Bejnordi
"""

import csv
import numpy as np
from datetime import datetime
from sklearn.metrics import roc_auc_score
import random
import math
import scipy.stats as st 
    
def readCSVtoDictionary(path):
    reader = csv.reader(open(path))
    
    GT_info = {}
    for row in reader:
        key = row[0]
        if key in GT_info:
            # implement your duplicate row handling here
            pass
        GT_info[key] = row[1:]  
    return GT_info

def get_size_type_lesions(gt_info):
    lesion_only_type = {}
    lesion_only_size = {}
    for case in gt_info:
        if gt_info[case][0].lower() == 'tumor' and 'Test_049' not in case:
            if gt_info[case][1].lower() == 'idc':
                lesion_only_type[case] = 'idc'
            else: 
                lesion_only_type[case] = 'non-idc'
            if gt_info[case][2].lower() == 'macro':
                lesion_only_size[case] = 'macro'
            else: 
                lesion_only_size[case] = 'micro'
    return lesion_only_type, lesion_only_size

def np_random_shuffle(a):
    random.seed(datetime.now())
    keys = a.keys()
    random.shuffle(keys)
    b = dict(zip(keys, a.values()))
    return b

def compute_auc_diffs(lesion_only_type, lesion_only_size, gt_info, pathol_probs_dic, do_permute = 1):
    pathol_probs_micro = []
    pathol_probs_macro = []
    pathol_probs_idc = []
    pathol_probs_nonidc = []
    gt_labels_micro = []
    gt_labels_macro = []
    gt_labels_idc = []
    gt_labels_nonidc = []
    if do_permute:
        lesion_only_type_perm = np_random_shuffle(lesion_only_type)
        lesion_only_size_perm = np_random_shuffle(lesion_only_size)
    else:
        lesion_only_type_perm = lesion_only_type
        lesion_only_size_perm = lesion_only_size

    for case in gt_info:
        if 'Test_049' not in case:
            if gt_info[case][0].lower() == 'normal':
                pathol_probs_micro.append(float(pathol_probs_dic[case][0]))
                pathol_probs_macro.append(float(pathol_probs_dic[case][0]))
                pathol_probs_idc.append(float(pathol_probs_dic[case][0]))
                pathol_probs_nonidc.append(float(pathol_probs_dic[case][0]))               
                gt_labels_micro.append(0)
                gt_labels_macro.append(0)
                gt_labels_idc.append(0)
                gt_labels_nonidc.append(0)
            else:
                if lesion_only_size_perm[case] == 'macro':                   
                    pathol_probs_macro.append(float(pathol_probs_dic[case][0]))
                    gt_labels_macro.append(1)
                else:
                    pathol_probs_micro.append(float(pathol_probs_dic[case][0]))
                    gt_labels_micro.append(1)
                if lesion_only_type_perm[case] == 'idc':                   
                    pathol_probs_idc.append(float(pathol_probs_dic[case][0]))
                    gt_labels_idc.append(1)
                else:
                    pathol_probs_nonidc.append(float(pathol_probs_dic[case][0]))
                    gt_labels_nonidc.append(1)
    
    auc_mic = roc_auc_score(gt_labels_micro, pathol_probs_micro)
    auc_mac = roc_auc_score(gt_labels_macro, pathol_probs_macro)
    auc_idc = roc_auc_score(gt_labels_idc, pathol_probs_idc)
    auc_nonidc = roc_auc_score(gt_labels_nonidc, pathol_probs_nonidc)
    return auc_mic, auc_mac, auc_idc, auc_nonidc

def compute_pValue(diff_lesion_size, diff_lesion_type,
                   all_mic_mac_aucs, all_idc_nonidc_aucs, num_permutations):
    p_value_s = np.sum(diff_lesion_size < map(abs, all_mic_mac_aucs)) / float(num_permutations)
    p_value_t = np.sum(diff_lesion_type < map(abs, all_idc_nonidc_aucs)) / float(num_permutations)
    return p_value_s, p_value_t    

if __name__ == "__main__":
    roc_pathol = r'...path_to_scores_of_pathologists.csv'
    pathol_probs_dic = readCSVtoDictionary(roc_pathol)
        
    gt_path = r'...path_to_reference_standard_data.csv'
    gt_info = readCSVtoDictionary(gt_path)

    permuted_labels = np.random.randint(2, size=len(gt_info))
    
    lesion_only_type, lesion_only_size = get_size_type_lesions(gt_info)

    num_permutations = 1000
    all_mic_mac_aucs = []
    all_idc_nonidc_aucs = []
    for i in range(num_permutations):
        auc_mic, auc_mac, auc_idc, auc_nonidc = compute_auc_diffs(lesion_only_type, \
            lesion_only_size, gt_info, pathol_probs_dic, do_permute = 1)
        all_mic_mac_aucs.append(auc_mac - auc_mic)
        all_idc_nonidc_aucs.append(auc_idc - auc_nonidc)                    
    
    # diffs of the auc for the distributions
    diff_mean_size = np.mean(all_mic_mac_aucs)
    diff_std_size = np.std(all_mic_mac_aucs)
    diff_mean_type = np.mean(all_idc_nonidc_aucs)
    diff_std_type = np.std(all_idc_nonidc_aucs)
    
    # compute true AUCs diffs  
    auc_mic_real, auc_mac_real, auc_idc_real, auc_nonidc_real = compute_auc_diffs(lesion_only_type, \
        lesion_only_size, gt_info, pathol_probs_dic, do_permute = 0)
    diff_lesion_size = auc_mac_real - auc_mic_real
    diff_lesion_type = auc_idc_real - auc_nonidc_real    
    
    pval_s, pval_t = compute_pValue(diff_lesion_size, diff_lesion_type,
                                    all_mic_mac_aucs, all_idc_nonidc_aucs,
                                    num_permutations)
                                    
    # bonferonni correctced *12
    print "the p-value: Size difference: ", pval_s*12
    print "the p-value: Type difference: ", pval_t*12