from test_results_analysis.TestResults import *
from src._utils import *
from data.NeuralTube import *
def calculate_neural_tube_decoding_maps_one_gene_group(encode_genes, data_path:str):
    """
    """
    neuraltube = NeuralTube(data_path=data_path,training=True, edge_trim=20)
    neuraltube_test_data = load_gastruloid_data(data_path)

    #train and test cell-independent
    decoding_sc = neuraltube.train_and_test_sc(neuraltube_test_data, encode_genes)
    gasrtu_results_sc = TestResults(decoding_sc, 'sc_gastru', neuraltube.means_sc, neuraltube.std_sc, encode_genes,
                                      edge_trim=EDGE_TRIM)
    gasrtu_results_sc.save(GASTRU_RES_DIR)

    # #train and test neighborhood-informed
    decoding_wn = neuraltube.train_and_test_wn(neuraltube_test_data, encode_genes)
    gastru_results_wn = TestResults(decoding_wn, 'wn_gastru', neuraltube.means_wn, neuraltube.covs_wn, encode_genes,
                                      edge_trim=EDGE_TRIM)
    gastru_results_wn.save(GASTRU_RES_DIR)