from data.WildTypeDrosophila import *



def wt_droso_load_preprocess_train_and_eval_neighborhood_informed(encode_genes):
    wt_droso = WildTypeDrosophilaData(training=True)
    test_data = load_wt_gap_test_data()
    # #train and test sc
    decoding_map_wn = wt_droso.train_and_test_wn(test_data, ENCODE_GENES)
    wt_droso_results = TestResults(decoding_map_wn, 'wn', wt_droso.means_wn, wt_droso.covs_wn, ENCODE_GENES)
    #wt_droso_results.plot_decoding_map('')
    wt_droso_results.save(DROSO_RES_DIR)


def calculate_wt_decoding_map(encode_genes):
    wt_droso = WildTypeDrosophilaData(None)
    test_data = load_wt_gap_test_data()
    decoding_map_sc = wt_droso.train_and_test_sc(test_data, encode_genes)
    decoding_map_wn = wt_droso.train_and_test_wn(test_data, encode_genes)

    pass