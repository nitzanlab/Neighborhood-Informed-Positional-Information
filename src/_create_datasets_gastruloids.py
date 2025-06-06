from data.Gastruloids import *
from test_results_analysis.TestResults import *
from src._utils import *

def calculate_gastru_decoding_maps_one_gene_group(encode_genes, data_path:str):
    """
    """
    gastru = Gastruloids(data_path=data_path,training=True)
    gastu_test_data = load_gastruloid_data(data_path)

    #train and test cell-independent
    decoding_sc = gastru.train_and_test_sc(gastu_test_data, encode_genes)
    gasrtu_results_sc = TestResults(decoding_sc, 'gastru', gastru.means_sc, gastru.std_sc, encode_genes,
                                      edge_trim=EDGE_TRIM)
    gasrtu_results_sc.save(GASTRU_RES_DIR)

    #train and test neighborhood-informed
    decoding_wn = wt_droso.train_and_test_wn(test_data, encode_genes)
    wt_droso_results_wn = TestResults(decoding_wn, 'wn_wt', wt_droso.means_wn, wt_droso.covs_wn, encode_genes,
                                      edge_trim=EDGE_TRIM)
    wt_droso_results_wn.save(DROSO_RES_DIR)