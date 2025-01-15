"""
Description:
    This file contains functions for loading the wild type and mutant background embryo data,
    including gap gene expression datasets and pair-rule gene expression datasets.

Key Functions:
    - load_data: Load raw data from a specified path.
"""

from data._preprocessing import *


def get_wt_pr_data():
    pair_rule_wt_path = MUTANT_PAIR_RULE_DIR_PATH + WT_PAIR_RULE_DATA + '.mat'
    wt_pair_rule_data = get_normalized_decode_data([pair_rule_wt_path],
                                                   WT_PAIR_RULE_NON_GENE_LIST,
                                                   ALL_PAIR_RULE_GENES,
                                                   PAIR_RULE_LOWER_TIME,
                                                   PAIR_RULE_UPPER_TIME)[PAIR_RULE_GENES]
    return wt_pair_rule_data


def load_wt_gap_test_data():
    wt_gap_gene_data = get_normalized_decode_data(WT_EMBRYO_FILE_NAMES, NON_GENES_LIST,
                                                  GAP_GENES, LOW_TIME, UPPER_TIME)
    #38,48
    return wt_gap_gene_data

def load_all_wt_droso_train_data(save=False, save_dir='', save_name=''): #combine_all_wt_38_48_without_test
    all_wt_datas = []
    for i, filename in enumerate(WT_EMBRYO_FILE_NAMES[:-1]):
        mutant_type_stained_with_wt = filename.split('dorsal_')[1].split('.mat')[0]
        if mutant_type_stained_with_wt == 'wt_bcd_nos_tsl':  # all morphogen signals knocked out, acts as a control and not as wt data
            continue
        data = load_data(filename, NON_GENES_LIST, GAP_GENES)
        wt_data = data[data['genotype'] == 1]
        if mutant_type_stained_with_wt=='wt': #leave out 38 wt embryos 40-44 as test set to compare to paper
            wt_data = wt_data[(wt_data['age']>UPPER_TIME) | (wt_data['age']<LOW_TIME)]
        all_wt_datas.append(wt_data)
    combined_df = pd.concat(all_wt_datas, ignore_index=True)
    normalized_combined_df = normalize_all_training_data(combined_df[GAP_GENES])
    if save:
        normalized_combined_df.to_csv(os.path.join(save_dir, save_name))
    return normalized_combined_df[GAP_GENES], combined_df[NON_GENES_LIST]


def get_wt_pr_data():
    pair_rule_wt_path = MUTANT_PAIR_RULE_DIR_PATH + WT_PAIR_RULE_DATA + '.mat'
    wt_pair_rule_data = get_normalized_decode_data([pair_rule_wt_path],
                                                   WT_PAIR_RULE_NON_GENE_LIST,
                                                   ALL_PAIR_RULE_GENES,
                                                   PAIR_RULE_LOWER_TIME,
                                                   PAIR_RULE_UPPER_TIME)[PAIR_RULE_GENES]
    return wt_pair_rule_data

def get_mean_wt_pr_per_position(edge_trim=None):
    wt_pr_data = get_wt_pr_data()
    full_wt_pr_data = reformat_exp_data_to_arr(wt_pr_data, edge_trim)
    mean_exp_wt_pr_per_position = np.mean(full_wt_pr_data, axis=0)
    mean_exp_wt_per_positions_normed = ((mean_exp_wt_pr_per_position - np.min(mean_exp_wt_pr_per_position, axis=0))
                                        / (np.max(mean_exp_wt_pr_per_position, axis=0) - np.min(
                mean_exp_wt_pr_per_position, axis=0)))
    return mean_exp_wt_per_positions_normed


def reformat_exp_data_to_arr(data, edge_trim=None):
    data_arr = np.array(data.values.tolist()).transpose(0, 2, 1)
    if edge_trim is not None:
        data_arr = data_arr[:, edge_trim:-edge_trim, :]
    return data_arr

##MUTANT data
def get_mutant_pr_data(mutant_type:str):
    mutant_path = MUTANT_PAIR_RULE_DIR_PATH + mutant_type + '.mat'
    mutant_pr_data = get_normalized_decode_data([mutant_path],
                                                PAIR_RULE_GENES_NON_GENE_LIST,
                                                ALL_PAIR_RULE_GENES,
                                                PAIR_RULE_LOWER_TIME,
                                                PAIR_RULE_UPPER_TIME)[PAIR_RULE_GENES]
    return mutant_pr_data



def get_mutant_gap_data(mutant_type:str):
    mutant_path = MUTANT_GAP_GENE_DIR_PATH + mutant_type + '.mat'
    mutant_gap_data = get_normalized_decode_data_mutant([mutant_path],NON_GENES_LIST,GAP_GENES,
                                             38,
                                             49)
    return mutant_gap_data

def get_normalized_decode_data_mutant(file_names, non_genes_list, genes, low_time, upper_time):
    data_decode = load_data(file_names[0], non_genes_list, genes)
    mutant_data = time_filter_data_mutant(data_decode, genes, low_time,
        upper_time)

    simulateneous_wt = filter_simultaneous_wt(data_decode, genes, low_time, upper_time)
    min_mean_exp, max_mean_exp = min_and_max_mean_gene_expression(simulateneous_wt,
                                                                  genes)
    data_decode = normalize_gene_exp(mutant_data, min_mean_exp, max_mean_exp, genes)
    return data_decode

def time_filter_data_mutant(data,genes, lower_thresh,upper_thresh):
    data = data.loc[
        (data['age'] >= lower_thresh) & (data['age'] <= upper_thresh)]
    data = data.loc[data['genotype'] == 2]
    return data[genes]

def filter_simultaneous_wt(data,genes, lower_thresh,upper_thresh):
    data = data.loc[
        (data['age'] >= lower_thresh) & (data['age'] <= upper_thresh)]
    data = data.loc[data['genotype'] == 1]
    return data[genes]


def normalize_all_training_data(training_data):
    min_mean_exp, max_mean_exp = min_and_max_mean_gene_expression(training_data, GAP_GENES)
    data_decode = normalize_gene_exp(training_data, np.array(min_mean_exp), np.array(max_mean_exp), GAP_GENES)
    return data_decode
