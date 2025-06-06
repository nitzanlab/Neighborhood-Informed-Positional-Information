from data.WildTypeDrosophila import *
from test_results_analysis.TestResults import *
from src._utils import *

def calculate_wt_droso_decoding_maps_one_gene_group(encode_genes):
    """
    This function calculates the wild type Drosophila embryos decoding map
    when decoding given one group of gap genes, a subset of the four gap genes.
    """
    wt_droso = WildTypeDrosophilaData(training=True)
    test_data = load_wt_gap_test_data()

    #train and test cell-independent
    decoding_sc = wt_droso.train_and_test_sc(test_data, encode_genes)
    wt_droso_results_sc = TestResults(decoding_sc, 'sc_wt', wt_droso.means_sc, wt_droso.std_sc, encode_genes,
                                      edge_trim=EDGE_TRIM)
    wt_droso_results_sc.save(DROSO_RES_DIR)

    #train and test neighborhood-informed
    decoding_wn = wt_droso.train_and_test_wn(test_data, encode_genes)
    wt_droso_results_wn = TestResults(decoding_wn, 'wn_wt', wt_droso.means_wn, wt_droso.covs_wn, encode_genes,
                                      edge_trim=EDGE_TRIM)
    wt_droso_results_wn.save(DROSO_RES_DIR)

def calculate_all_gene_subset_decoding_maps_WT():
    """
    This function calculates the decoding maps for all of the wild type embryos given
    each subset of the four gap genes
    :return:
    """
    gap_genes_subsets = get_all_gap_gene_subsets_list()
    wt_droso = WildTypeDrosophilaData(training=True, edge_trim=EDGE_TRIM)
    test_data = load_wt_gap_test_data()
    for gap_genes_subset in gap_genes_subsets:
        print(gap_genes_subset)
        decoding_sc = wt_droso.train_and_test_sc(test_data, gap_genes_subset)
        wt_droso_results_sc = TestResults(decoding_sc, 'sc_wt', wt_droso.means_sc, wt_droso.std_sc, gap_genes_subset,
                                          edge_trim=EDGE_TRIM)
        wt_droso_results_sc.save(DROSO_RES_DIR)

        decoding_wn = wt_droso.train_and_test_wn(test_data, gap_genes_subset)
        wt_droso_results_wn = TestResults(decoding_wn, 'wn_wt', wt_droso.means_wn, wt_droso.covs_wn, gap_genes_subset,
                                          edge_trim=EDGE_TRIM)
        wt_droso_results_wn.save(DROSO_RES_DIR)

def calculate_all_mutant_decoding_maps():
    """
    This function calculates the decoding map for all of the mutant embryos of one
    mutant background
    """
    wt_droso = WildTypeDrosophilaData(training=True, edge_trim=EDGE_TRIM)
    wt_droso.train_sc(GAP_GENES)
    wt_droso.train_wn(GAP_GENES)

    for mutant_type in MUTANT_TYPES:
        gap_data = get_mutant_gap_data(mutant_type)
        decoding_wn = wt_droso.test_wn(gap_data, GAP_GENES)
        decoding_sc = wt_droso.test_sc(gap_data, GAP_GENES)
        results_wn = TestResults(decoding_wn, f'{mutant_type}_wn', wt_droso.means_wn, wt_droso.covs_wn, GAP_GENES,
                                 edge_trim=EDGE_TRIM)
        results_wn.save(DROSO_RES_DIR)
        results_sc = TestResults(decoding_sc, f'{mutant_type}_sc', wt_droso.means_wn, wt_droso.covs_wn, GAP_GENES,
                                 edge_trim=EDGE_TRIM)
        results_sc.save(DROSO_RES_DIR)