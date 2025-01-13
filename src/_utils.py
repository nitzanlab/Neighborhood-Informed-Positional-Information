from src._constants import *
from test_results_analysis.TestResults import *
from data._data import *
def set_style():
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'figure.titlesize': 26, 'figure.titleweight': 'bold', 'axes.titlesize': 22,
                         'axes.titleweight': "bold", 'axes.labelsize': 24,'axes.labelweight': 'bold',
                         "ytick.labelsize": 28, "xtick.labelsize": 28, 'legend.fontsize': 24,'font.family': 'Candara,  Arial'})
    plt.rcParams.update({'figure.figsize': (8, 6)})
    plt.rcParams.update({'savefig.dpi': 300})


def caclulate_positional_error_per_decoding_map_MAP_positions(decoding_type):
    wt_decoding_res = TestResults.from_pickle(DROSO_RES_DIR, decoding_type, GAP_GENES)
    wt_cov = wt_decoding_res.train_std
    decoding_map = wt_decoding_res.normalized_decoding_map
    MAP_positions = np.argmax(decoding_map, axis=1)
    mean_exp = wt_decoding_res.train_mean

    mean_exp_slopes = np.diff(mean_exp,axis=0)
    mean_exp_extended = np.vstack([mean_exp_slopes, mean_exp_slopes[-1]])

    num_genes = mean_exp_slopes.shape[1]
    num_embryos, num_positions = MAP_positions.shape[0], MAP_positions.shape[1] - 1
    MAP_pos_mean_exp_slopes = np.zeros((num_embryos, num_positions,num_genes))

    # Iterate over embryos to extract slopes corresponding to MAP positions
    for embryo in range(num_embryos):
        MAP_pos_mean_exp_slopes[embryo] = mean_exp_extended[MAP_positions[embryo, 1:]]

    wt_covs_MAP = wt_cov[MAP_positions][:,1:,:,:]
    num_positions = MAP_pos_mean_exp_slopes.shape[1]
    num_embryos = MAP_pos_mean_exp_slopes.shape[0]
    position_error = np.zeros((num_embryos,num_positions))
    for embryo in range(num_embryos):
        for pos in range(num_positions):
            position_error[embryo,pos] = MAP_pos_mean_exp_slopes[embryo,pos,:]@np.linalg.inv(wt_covs_MAP[embryo,pos,:,:])@MAP_pos_mean_exp_slopes[embryo,pos,:]
    pos_err = np.sqrt((1/(position_error)))/num_positions
    return pos_err

def calculate_positional_inf(sigma_x):
    return np.log2(1/(sigma_x*np.sqrt(2*np.pi*np.e)))


def get_all_gap_gene_subsets_list(save=False):
    all_decoding_gene_subsets = []
    for r in range(1, len(GAP_GENES) + 1):
        for subset_decode_genes in combinations(GAP_GENES, r):
            all_decoding_gene_subsets.append(list(subset_decode_genes))
    if save:
        filename = os.path.join(DROSO_RES_DIR, 'decoding_gene_subsets.pkl')
        with open(filename, 'wb') as file:
            pickle.dump(all_decoding_gene_subsets, file)
    return all_decoding_gene_subsets

def get_decoding_maps():
    sc_map = TestResults.from_pickle(DROSO_RES_DIR, 'sc_wt_sigmax', GAP_GENES)
    sc_map.normalized_decoding_map = sc_map.raw_test_results[:, 1:-1, 1:-1]
    wn_map = TestResults.from_pickle(DROSO_RES_DIR, 'wn_wt_sigmax', GAP_GENES)
    wn_map.normalized_decoding_map = wn_map.raw_test_results
    return sc_map, wn_map

def get_pr_expected_expression_profiles(decoding_map):
    mean_wt_pr_exp = get_mean_wt_pr_per_position(EDGE_TRIM)[1:-1, :]
    pr_predictions = decoding_map.expected_pr_exp(mean_wt_pr_exp)
    return pr_predictions

def compare_expression_correlations(mutant_type, wn_pred, sc_pred):
    mutant_pr_exp = get_mutant_pr_data(mutant_type)
    ground_truth = reformat_exp_data_to_arr(mutant_pr_exp, EDGE_TRIM)[:, 1:-1, :]

    results = {}
    all_results = {}

    for i, gene in enumerate(PAIR_RULE_GENES):
        # Step 1: Compute ground truth mean slope
        gt_gene = ground_truth[:, :, i]
        smooth_gt = gaussian_filter1d(gt_gene, sigma=20, axis=1)
        mean_gt_gene = np.mean(smooth_gt, axis=0)

        pred_wn_gene = wn_pred[:, :, i]  # Extract one feature
        smoothed_preds_wn = gaussian_filter1d(pred_wn_gene, sigma=20, axis=1)

        pred_sc_gene = sc_pred[:, :, i]  # Extract one feature
        smoothed_preds_sc = gaussian_filter1d(pred_sc_gene, sigma=20, axis=1)

        num_pred_embryos = pred_wn_gene.shape[0]

        corr_sc = np.zeros(num_pred_embryos)
        corr_wn = np.zeros(num_pred_embryos)
        for j in range(num_pred_embryos):
            corr_sc[j] = np.corrcoef(smoothed_preds_sc[j], mean_gt_gene)[0, 1]
            corr_wn[j] = np.corrcoef(smoothed_preds_wn[j], mean_gt_gene)[0, 1]

        num_gt_embryos = gt_gene.shape[0]
        corr_gt = np.zeros(num_gt_embryos)
        for k in range(num_gt_embryos):
            corr_gt[k] = np.corrcoef(gt_gene[k], mean_gt_gene)[0, 1]
        results[f'{gene}'] = {'mean_Correlation_sc': np.mean(corr_sc),
                              'std_Correlation_sc': np.std(corr_sc),
                              'mean_Correlation_wn': np.mean(corr_wn), 'std_Correlation_wn': np.std(corr_wn)}
        all_results[f'{gene} CI'] = corr_sc
        all_results[f'{gene} NI'] = corr_wn
        all_results[f'{gene} GT'] = corr_gt

    return results, all_results


def calculate_mutant_pr_pred(mutant_type):
    wn_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_wn', GAP_GENES)
    sc_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_sc', GAP_GENES)

    mean_wt_pr_exp = get_mean_wt_pr_per_position(EDGE_TRIM)
    mean_mutant_pr_exp = get_mutant_pr_data(mutant_type)
    full_mutant_pr_data = reformat_exp_data_to_arr(mean_mutant_pr_exp, EDGE_TRIM)
    mean_exp_mutant_pr_per_position = np.mean(full_mutant_pr_data, axis=0)
    mean_exp_mutant_per_positions_normed = (
                (mean_exp_mutant_pr_per_position - np.min(mean_exp_mutant_pr_per_position, axis=0))
                / (np.max(mean_exp_mutant_pr_per_position, axis=0) - np.min(
            mean_exp_mutant_pr_per_position, axis=0)))

    sc_pr_pred = sc_map.expected_pr_exp(mean_wt_pr_exp)
    wn_pr_pred = wn_map.expected_pr_exp(mean_wt_pr_exp[1:-1, :])
    return sc_pr_pred, wn_pr_pred




def calculate_mutant_pr_pred_err(mutant_type, sc_pr_pred, wn_pr_pred,mean_exp_mutant_pr_per_position):
    mean_mutant_pr_exp = get_mutant_pr_data(mutant_type)
    full_mutant_pr_data = reformat_exp_data_to_arr(mean_mutant_pr_exp, EDGE_TRIM)
    sc_pr_pred_err = pr_prediction_err(sc_pr_pred, mean_exp_mutant_pr_per_position)[:, 1:-1, :]
    wn_pr_pred_err = pr_prediction_err(wn_pr_pred, mean_exp_mutant_pr_per_position[1:-1, :])
    wt_pr_exp = get_wt_pr_data()

    std_pr_gene_per_pos = np.std(full_mutant_pr_data, axis=0)[1:-1, :]

    sc_weighted_pr_err = (sc_pr_pred_err) / (std_pr_gene_per_pos)  # [:,200:-100,:]
    wn_weighted_pr_err = (wn_pr_pred_err) / (std_pr_gene_per_pos)
    return sc_weighted_pr_err, wn_weighted_pr_err
