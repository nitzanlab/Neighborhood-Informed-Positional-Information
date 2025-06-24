from test_results_analysis.TestResults import *
from src._utils import *

def calculate_gastru_decoding_maps_one_gene_group(encode_genes, data_path:str):
    """
    """
    gastru = Gastruloids(data_path=data_path,training=True, edge_trim=20)
    gastu_test_data = format_gastru_like_droso(load_gastruloid_data(data_path))

    #train and test cell-independent
    decoding_sc = gastru.train_and_test_sc(gastu_test_data, encode_genes)
    gasrtu_results_sc = TestResults(decoding_sc, 'sc_gastru', gastru.means_sc, gastru.std_sc, encode_genes,
                                      edge_trim=EDGE_TRIM)
    gasrtu_results_sc.save(GASTRU_RES_DIR)

    # #train and test neighborhood-informed
    decoding_wn = gastru.train_and_test_wn(gastu_test_data, encode_genes)
    gastru_results_wn = TestResults(decoding_wn, 'wn_gastru', gastru.means_wn, gastru.covs_wn, encode_genes,
                                      edge_trim=EDGE_TRIM)
    gastru_results_wn.save(GASTRU_RES_DIR)