"""
Figure 1 visuals
"""
import matplotlib.pyplot as plt

from src._imports import *
from src._constants import *
from data._preprocessing import *
from data._data import *
from data.WildTypeDrosophila import *
from src._utils import *
#### panels b and c, gap gene and pair rule genes pairwise correlations:
def gene_expression_pairwise_correlation_plots():
    """
    This function plots the pairwise correlation per gene of all of the gap genes together and all of
    the pair rule genes together.
    """
    plot_gap_gene_pairwise_correlations()
    plot_pair_rule_gene_pairwise_correlations()


def plot_pair_rule_gene_pairwise_correlations():
    """
    This function plots the pairwise pair rule gene expression correlations over all of the positions
    """
    pr_corr = pairwise_correlation_pair_rule_genes()
    bins = np.linspace(0.975, 1, 100)
    for i in range(pr_corr.shape[0]):  # per gene
        plt.hist(
            pr_corr[i],
            bins=bins,  # Use density to calculate the density
            weights=np.ones_like(pr_corr[i]) / len(pr_corr[i]),  # Normalize by total count
            color=PAIR_RULE_COLORS[PAIR_RULE_GENES[i]],
            label=PAIR_RULE_GENES[i],
            alpha=0.5
        )

    plt.legend(fontsize=20, loc='upper left')
    plt.xlabel('pair-rule genes pairwise correlation')
    plt.ylabel('relative density')

    plt.ylim(0, .1)
    plt.xlim(0.97,1)
    plt.tight_layout()
    plt.show()
    pass

def plot_gap_gene_pairwise_correlations():
    """
    This function plots the pair wise correlation  in pair rule gene expression in the wild type embryo over positions
    """
    gg_corr = pairwise_correlations_gap_genes()
    bins = np.linspace(0.975, 1, 100)
    for i in range(gg_corr.shape[0]):  # per gene
        plt.hist(
            gg_corr[i],
            bins=bins,  # Use density to calculate the density
            weights=np.ones_like(gg_corr[i]) / len(gg_corr[i]),  # Normalize by total count
            color=GAP_GENE_COLORS[GAP_GENES[i]],
            label=GAP_GENES[i],
            alpha=0.5
        )
    plt.legend(fontsize=20, loc='upper left')
    plt.xlabel('gap gene pairwise correlation')
    plt.ylabel('relative density')
    plt.tight_layout()
    plt.xlim(0.97, 1)
    plt.ylim(0, .1)
    plt.tight_layout()
    plt.show()
def pairwise_correlations_gap_genes():
    """
    This function calculates the pairwise correlations in gap genes
    """
    wt_gap_gene = load_wt_gap_test_data()
    wt_arr = reformat_exp_data_to_arr(wt_gap_gene)[:,EDGE_TRIM:-EDGE_TRIM,:]
    gene_names = np.array(wt_gap_gene.columns)
    all_genes_corrs = []
    for i in range(wt_arr.shape[2]):
        correlations = []
        wt_arr_gene = wt_arr[:, :, i].reshape((wt_arr.shape[0], wt_arr.shape[1]))
        for j in range(wt_arr_gene.shape[1] - 1):
            corr = np.corrcoef(wt_arr_gene[:, j], wt_arr_gene[:, j + 1])[0, 1]  # Correlation between column i and i+1
            correlations.append(corr)
        all_genes_corrs.append(np.array(correlations))
    return np.array(all_genes_corrs)

def pairwise_correlation_pair_rule_genes():
    """
        This function measures the pairwise correlations in pair rules gene expression
    """
    wt_pr = get_wt_pr_data()
    wt_arr_pr = reformat_exp_data_to_arr(wt_pr)[:, EDGE_TRIM:-EDGE_TRIM, :]
    all_genes_corrs = []
    for i in range(wt_arr_pr.shape[2]):
        correlations = []
        wt_arr_gene = wt_arr_pr[:, :, i].reshape((wt_arr_pr.shape[0], wt_arr_pr.shape[1]))
        for j in range(wt_arr_gene.shape[1] - 1):
            corr = np.corrcoef(wt_arr_gene[:, j], wt_arr_gene[:, j + 1])[0, 1]  # Correlation between column i and i+1
            correlations.append(corr)
        all_genes_corrs.append(np.array(correlations))
    return np.array(all_genes_corrs)

### figure 1d, position information in bits panel
def plot_positional_information_in_bits_empirical_std(wn_stds, sc_stds):
    """
    This function plots the positional information in bits based on the empirical posterior position decoding distribution
    """
    positions = np.linspace(0.1, 0.9, sc_stds.shape[1])

    plt.plot(positions, np.log2(N) * np.ones_like(positions), color='gray', linestyle='--',
             label='Unique cell specification')
    plt.fill_between(x=positions, y1=np.log2(N + N_STD), y2=np.log2(N - N_STD), color='gray', alpha=0.5)

    mean_pos_inf_sc = np.mean(calculate_positional_inf(sc_stds), axis=0)
    mean_pos_inf_wn = np.mean(calculate_positional_inf(wn_stds), axis=0)
    std_pos_inf_sc = np.std(calculate_positional_inf(sc_stds), axis=0)
    std_pos_inf_wn = np.std(calculate_positional_inf(wn_stds), axis=0)
    plt.plot(positions, mean_pos_inf_sc, label=DECODER_NAMES['sc'], color=DECODER_TYPE_COLOR['sc'])
    plt.plot(positions, mean_pos_inf_wn, label=DECODER_NAMES['wn'], color=DECODER_TYPE_COLOR['wn'])
    plt.fill_between(positions, mean_pos_inf_sc + std_pos_inf_sc, mean_pos_inf_sc - std_pos_inf_sc,
                     color=DECODER_TYPE_COLOR['sc'], alpha=0.3)
    plt.fill_between(positions, mean_pos_inf_wn + std_pos_inf_wn,
                     mean_pos_inf_wn - std_pos_inf_wn, color=DECODER_TYPE_COLOR['wn'], alpha=0.3)

    plt.legend(loc='lower right')
    plt.xlabel(POSITION_X_LABEL)
    plt.ylabel('positional information in bits')
    plt.ylim(0, 10)
    plt.xlim(0.1, 0.9)
    plt.tight_layout()
    plt.show()


def plot_positional_information_in_bits():
    """
    This function measures and plots the positional information in bits based on the position estimate error
    """
    sigma_x_wn = caclulate_positional_error_per_decoding_map_MAP_positions(f'wn_wt')  # vs 'wn'
    sigma_x_sc = caclulate_positional_error_per_decoding_map_MAP_positions(f'sc_wt')


    mean_wn_stds = np.mean(sigma_x_wn, axis=0)
    std_wn_stds = np.std(sigma_x_wn, axis=0)

    mean_sc_stds = np.mean(sigma_x_sc, axis=0)[1:-1]
    std_sc_stds = np.std(sigma_x_sc, axis=0)[1:-1]

    positions = np.linspace(POSITIONS_START, POSITIONS_END, len(mean_sc_stds))

    plt.plot(positions, mean_sc_stds, label=DECODER_NAMES['sc'], color=DECODER_TYPE_COLOR['sc'])
    plt.plot(positions, mean_wn_stds, label=DECODER_NAMES['wn'], color=DECODER_TYPE_COLOR['wn'])

    plt.fill_between(positions, mean_sc_stds - std_sc_stds, mean_sc_stds + std_sc_stds, color=DECODER_TYPE_COLOR['sc'],
                     alpha=0.5)
    plt.fill_between(positions, mean_wn_stds - std_wn_stds, mean_wn_stds + std_wn_stds, color=DECODER_TYPE_COLOR['wn'],
                     alpha=0.5)
    plt.xlabel(POSITION_X_LABEL)
    plt.ylabel('posterior standard deviation')
    plt.ylim(0, 0.05)
    plt.tight_layout()
    plt.show()

    positions = np.linspace(0.1, 0.9, sigma_x_sc.shape[1])
    plt.plot(positions,np.log2(N)*np.ones_like(positions), color='gray', linestyle='--', label='Unique cell specification')
    plt.fill_between(x=positions,y1=np.log2(N+N_STD), y2=np.log2(N-N_STD),  color='gray', alpha=0.5)

    mean_pos_inf_sc = np.mean(calculate_positional_inf(sigma_x_sc),axis=0)
    mean_pos_inf_wn = np.mean(calculate_positional_inf(sigma_x_wn),axis=0)
    std_pos_inf_sc = np.std(calculate_positional_inf(sigma_x_sc),axis=0)
    std_pos_inf_wn = np.std(calculate_positional_inf(sigma_x_wn),axis=0)
    plt.plot(positions, mean_pos_inf_sc, label=DECODER_NAMES['sc'], color=DECODER_TYPE_COLOR['sc'])
    plt.plot(positions[1:-1], mean_pos_inf_wn, label=DECODER_NAMES['wn'], color=DECODER_TYPE_COLOR['wn'])
    plt.fill_between(positions, mean_pos_inf_sc+ std_pos_inf_sc, mean_pos_inf_sc- std_pos_inf_sc, color=DECODER_TYPE_COLOR['sc'],alpha=0.3)
    plt.fill_between(positions[1:-1], mean_pos_inf_wn + std_pos_inf_wn,
                     mean_pos_inf_wn - std_pos_inf_wn, color=DECODER_TYPE_COLOR['wn'], alpha=0.3)

    plt.legend(loc='lower right')
    plt.xlabel(POSITION_X_LABEL)
    plt.ylabel('positional information in bits')
    plt.ylim(0,10)
    plt.xlim(0.1,0.9)
    plt.tight_layout()
    plt.show()

def plot_decoding_maps_comparison_wt(genes:list, vmaxs, xlim):
    """
    This function plots the comparison between the decoding map for neighborhood-informed decoding
    and the cell-independent decoding in wild type embryos
    """
    wt_droso = WildTypeDrosophilaData(training=True)
    wt_droso.plot_gene_exp_over_positions(genes)
    plot_decoding_maps([genes], WT_DECODE_TYPES,vmaxs, xlim)

def plot_all_decoding_maps_comparison_wt():
    """
    This function plots the comparison between the decoding maps when decoding based on multiple subsets of gap gene expression
    separately
    """
    for gene_subset in [['Kr'],['Kr','Hb'],['Kr','Gt','Hb'],GAP_GENES]:
        print(gene_subset)
        plot_decoding_maps([gene_subset], WT_DECODE_TYPES,VMAXS_THREE_GENES, xlim=False)

def plot_decoding_maps(genes_subset_list, types, vmaxs, xlim):
    """
    This function plots the decoding maps of the position posterior distribution
    """
    for i,gene_subset in enumerate(genes_subset_list):
        for type in types:
            decoding_map_obj = TestResults.from_pickle(DROSO_RES_DIR, type, gene_subset)
            #decoding_map_obj.normalize_decoding_map()
            decoding_map_obj.plot_decoding_map(' '.join(gene_subset), vmax=vmaxs[i], xlim=xlim)


def plot_position_posterior_std_all_gene_subsets():
    """
    This function plots the position standard deviation when decoding based on each of the
    gap gene subsets.
    """
    all_gene_subsets = get_all_gap_gene_subsets_list()
    sc_mean_errs = []
    sc_std_errs = []
    wn_mean_errs = []
    wn_std_errs = []
    genes_as_label = []
    error_by_num_genes_sc = {1: [], 2: [], 3: [], 4: []}
    error_by_num_genes_wn = {1: [], 2: [], 3: [], 4: []}
    t_tests_pvals = {1: [], 2: [], 3: [], 4: []}
    mean_t_test_pvals = {1: [], 2: [], 3: [], 4: []}
    for gene_subset in all_gene_subsets:
        num_genes = len(gene_subset)
        genes_as_label.append('\n'.join(gene_subset))
        sc_map = TestResults.from_pickle(DROSO_RES_DIR, 'sc_wt', gene_subset)
        sc_map.normalized_decoding_map = sc_map.normalized_decoding_map[:, 1:-1, 1:-1]
        wn_map = TestResults.from_pickle(DROSO_RES_DIR, 'wn_wt', gene_subset)

        wn_weighted_err = wn_map.measure_weighted_dist_prediction_error()
        sc_weighted_err = sc_map.measure_weighted_dist_prediction_error()
        mean_wn_err = np.mean(wn_weighted_err)
        std_wn_err = np.std(wn_weighted_err)
        mean_sc_err = np.mean(sc_weighted_err)
        std_sc_err = np.std(sc_weighted_err)

        sc_mean_errs.append(mean_sc_err)
        sc_std_errs.append(std_sc_err)

        wn_mean_errs.append(mean_wn_err)
        wn_std_errs.append(std_wn_err)

        error_by_num_genes_sc[num_genes].append(mean_sc_err)
        error_by_num_genes_wn[num_genes].append(mean_wn_err)
        mean_t_per_embryo = 0
        for i in range(sc_weighted_err.shape[0]):
            t_stat, p_value = ttest_ind(sc_weighted_err.flatten(), wn_weighted_err.flatten(), equal_var=False)
            mean_t_per_embryo +=  1 if p_value<=0.05 else 1
        t_tests_pvals[len(gene_subset)].append(mean_t_per_embryo/sc_weighted_err.shape[0])

    for i in range(1,5):
        mean_t_test_pvals[i] = np.mean(t_tests_pvals[i])
    print(f' mean t test vals: {mean_t_test_pvals}')
    # plot_bar_comparison(wn_mean_errs, wn_std_errs, sc_mean_errs, sc_std_errs, genes_as_label, label_axes=True, x_label='number of genes used for decoding position',y_label='mean position uncertainty')
    error_by_num_genes_to_plot_vals(error_by_num_genes_wn, error_by_num_genes_sc, label_axes=True,
                                    x_label='number of genes used for\n position decoding',
                                    y_label='mean position uncertainty')

def get_ttest_error_by_nums(error_by_num_genes_wn, error_by_num_genes_sc):
    ttest_results = []
    for i in range(1,len(error_by_num_genes_wn)+1):
        # Perform independent t-test for the ith entries
        t_stat, p_value = ttest_ind(error_by_num_genes_wn[i], error_by_num_genes_sc[i], equal_var=False)
        ttest_results.append((t_stat, p_value))
    print(ttest_results)
    return ttest_results
def error_by_num_genes_to_plot_vals(error_by_num_genes_wn, error_by_num_genes_sc,label_axes=False, x_label=None, y_label=None ):
    """
    This function plots the errors as a bar plot where each bar is the mean over all of the gene subsets of the same size.
    """

    wn_mean_num_genes = []
    wn_std_num_genes = []
    sc_mean_num_genes = []
    sc_std_num_genes = []
    for key in error_by_num_genes_sc.keys():
        wn_mean_num_genes.append(np.mean(np.array(error_by_num_genes_wn[key])))
        wn_std_num_genes.append(np.std(np.array(error_by_num_genes_wn[key])))
        sc_mean_num_genes.append(np.mean(np.array(error_by_num_genes_sc[key])))
        sc_std_num_genes.append(np.std(np.array(error_by_num_genes_sc[key])))
    plot_bar_comparison(wn_mean_num_genes, wn_std_num_genes, sc_mean_num_genes, sc_std_num_genes, np.arange(1, 5),label_axes=label_axes, x_label=x_label, y_label=y_label)

def plot_bar_comparison(wn_mean_errs, wn_std_errs, sc_mean_errs, sc_std_errs, labels,data_type='',label_axes=False, x_label=None, y_label=None):
    """
    This is a generic function for plot bar plots
    """
    lower_wn_err = sc_std_errs
    lower_sc_err = sc_std_errs

    barWidth = 0.25
    # setting position of bar on x axis
    br1 = np.arange(len(labels))
    br2 = br1 + barWidth
    plt.bar(br2, sc_mean_errs, yerr=[lower_sc_err, sc_std_errs], label=DECODER_NAMES["sc"], color=DECODER_TYPE_COLOR['sc'],width=barWidth)
    plt.bar(br1, wn_mean_errs, yerr=[lower_wn_err, wn_std_errs], label=DECODER_NAMES["wn"],color=DECODER_TYPE_COLOR['wn'], width=barWidth)

    if label_axes:
        plt.xlabel(x_label)
        plt.ylabel(y_label)

    plt.xticks(br1 + (0.5) * barWidth, labels)
    #plt.title(data_type)
    plt.legend(loc="upper right")
    #plt.ylim(0,0.4)
    plt.tight_layout()
    plt.subplots_adjust(left=0.2)
    plt.show()


def plot_position_MAP_error_all_gene_subsets():
    """
    This function plots the MAP position decoding error for all of the gap gene subsets
    """
    all_gene_subsets = get_all_gap_gene_subsets_list()
    sc_mean_errs = []
    sc_std_errs = []
    wn_mean_errs = []
    wn_std_errs = []
    genes_as_label = []
    error_by_num_genes_sc = {1: [], 2: [], 3: [], 4: []}
    error_by_num_genes_wn = {1: [], 2: [], 3: [], 4: []}
    t_tests_pvals = {1: [], 2: [], 3: [], 4: []}
    mean_t_test_pvals = {1: [], 2: [], 3: [], 4: []}
    for gene_subset in all_gene_subsets:
        num_genes = len(gene_subset)
        genes_as_label.append('\n'.join(gene_subset))
        sc_map = TestResults.from_pickle(DROSO_RES_DIR, 'sc_wt', gene_subset)
        sc_map.normalized_decoding_map = sc_map.normalized_decoding_map[:, 1:-1, 1:-1]
        wn_map = TestResults.from_pickle(DROSO_RES_DIR, 'wn_wt', gene_subset)

        wn_weighted_err = wn_map.calculate_MAP_error()
        sc_weighted_err = sc_map.calculate_MAP_error()
        mean_wn_err = np.mean(wn_weighted_err)
        std_wn_err = np.std(wn_weighted_err)
        mean_sc_err = np.mean(sc_weighted_err)
        std_sc_err = np.std(sc_weighted_err)

        sc_mean_errs.append(mean_sc_err)
        sc_std_errs.append(std_sc_err)

        wn_mean_errs.append(mean_wn_err)
        wn_std_errs.append(std_wn_err)

        error_by_num_genes_sc[num_genes].append(mean_sc_err)
        error_by_num_genes_wn[num_genes].append(mean_wn_err)
        mean_t_per_embryo = 0
        for i in range(sc_weighted_err.shape[0]):
            t_stat, p_value = ttest_ind(sc_weighted_err.flatten(), wn_weighted_err.flatten(), equal_var=False)
            mean_t_per_embryo += 1 if p_value<=0.05 else 0
        t_tests_pvals[len(gene_subset)].append(mean_t_per_embryo/sc_weighted_err.shape[0])

    for i in range(1,5):
        mean_t_test_pvals[i] = np.mean(t_tests_pvals[i])
    print(f' mean t test vals: {mean_t_test_pvals}')

    error_by_num_genes_to_plot_vals(error_by_num_genes_wn, error_by_num_genes_sc, label_axes=True,
                                    x_label='number of genes used for\n position decoding',
                                    y_label='mean MAP position\n prediction error')


####################### WT pair rule prediction related functions ##########################################
def pair_rule_prediction_across_positions():
    """
    This function calculates and then plots the pair rule gene expression predictions across the AP axis
    """
    sc_map, wn_map = get_wt_decoding_maps()
    sc_pr_pred = get_pr_expected_expression_profiles(sc_map)
    wn_pr_pred = get_pr_expected_expression_profiles(wn_map)
    wt_pr_exp = get_wt_pr_data()
    wt_pr_arr = reformat_exp_data_to_arr(wt_pr_exp, EDGE_TRIM)
    plot_predicted_prs_across_positions(wn_pr_pred, sc_pr_pred, wt_pr_arr[:, 1:-1, :])


def plot_predicted_prs_across_positions(wn_pr_pred, sc_pr_pred, wt_pr_data):
    """
    This function plots the predicted pair rule expression across positions
    """
    mean_pr_pred = np.mean(wt_pr_data, axis=0)
    mean_pr_pred_normed = ((mean_pr_pred - np.min(mean_pr_pred, axis=0))
                                        / (np.max(mean_pr_pred, axis=0) - np.min(
                mean_pr_pred, axis=0)))
    std_pr_pred = np.std(wn_pr_pred, axis=0)


    mean_prediction_sc = np.mean(sc_pr_pred, axis=0)
    mean_prediction_wn = np.mean(wn_pr_pred, axis=0)
    std_prediction_sc = np.std(sc_pr_pred, axis=0)
    std_prediction_wn = np.std(wn_pr_pred, axis=0)
    x_pos = np.linspace(0.1,0.9, wn_pr_pred.shape[1])

    for i, gene in enumerate(PAIR_RULE_GENES):
        mean_sc_gene = mean_prediction_sc[:, i]
        std_sc_gene = std_prediction_sc[:, i]
        mean_wn_gene = mean_prediction_wn[:, i]
        std_wn_gene = std_prediction_wn[:, i]
        plt.plot(x_pos, mean_pr_pred_normed[:, i], label=f'{gene} mean expression', color='red', linewidth=2.5)
        plt.fill_between(x_pos, np.maximum(mean_pr_pred_normed[:,i]-std_pr_pred[:,i],0), mean_pr_pred_normed[:,i]+std_pr_pred[:,i], color='red', alpha=0.2)
        plt.plot(x_pos, mean_sc_gene, label=DECODER_NAMES['sc'], color=DECODER_TYPE_COLOR['sc'])
        plt.fill_between(x_pos, np.maximum(mean_sc_gene-std_sc_gene,0), mean_sc_gene+std_sc_gene, color=DECODER_TYPE_COLOR['sc'], alpha=0.2)
        plt.xlabel(POSITION_X_LABEL)
        plt.ylabel(EXP_Y_LABEL)
        plt.legend()
        plt.tight_layout()
        custom_yticks = np.arange(0, 1.2, 0.2)
        plt.yticks(custom_yticks)
        plt.show()
        plt.plot(x_pos, mean_pr_pred_normed[:, i], label=f'{gene} mean expression', color='red', linewidth=2.5)
        plt.fill_between(x_pos, np.maximum(mean_pr_pred_normed[:, i] - std_pr_pred[:, i], 0),
                         mean_pr_pred_normed[:, i] + std_pr_pred[:, i], color='red', alpha=0.2)
        plt.plot(x_pos,mean_wn_gene, label=DECODER_NAMES['wn'], color=DECODER_TYPE_COLOR['wn'])
        plt.fill_between(x_pos,np.maximum(mean_wn_gene-std_wn_gene,0), mean_wn_gene+std_wn_gene, alpha=0.4, color=DECODER_TYPE_COLOR['wn'])
        plt.xlabel(POSITION_X_LABEL)
        custom_yticks = np.arange(0, 1.2, 0.2)
        plt.yticks(custom_yticks)
        plt.ylabel(EXP_Y_LABEL)
        plt.legend()
        plt.tight_layout()
        plt.show()



def plot_all_pr_genes_mean_prediction_error():
    """
    This function plots all of the pair rule genes mean prediction errors
    """
    sc_map, wn_map = get_wt_decoding_maps()
    sc_pr_pred = get_pr_expected_expression_profiles(sc_map)
    wn_pr_pred = get_pr_expected_expression_profiles(wn_map)

    mean_wt_pr_exp = get_mean_wt_pr_per_position(EDGE_TRIM)[1:-1, :]
    sc_pr_pred_err = pr_prediction_err(sc_pr_pred, mean_wt_pr_exp)
    wn_pr_pred_err = pr_prediction_err(wn_pr_pred, mean_wt_pr_exp)

    wt_pr_exp = get_wt_pr_data()
    wt_pr_arr = reformat_exp_data_to_arr(wt_pr_exp, EDGE_TRIM)
    std_pr_gene_per_pos = np.std(wt_pr_arr, axis=0)[1:-1, :]

    sc_weighted_pr_err = (sc_pr_pred_err / std_pr_gene_per_pos)
    wn_weighted_pr_err = (wn_pr_pred_err / std_pr_gene_per_pos)

    mean_sc_err = np.mean(np.mean(sc_weighted_pr_err, axis=1), axis=0)
    std_sc_err = np.std(np.mean(sc_weighted_pr_err, axis=1), axis=0)
    mean_wn_err = np.mean(np.mean(wn_weighted_pr_err, axis=1), axis=0)
    std_wn_err = np.std(np.mean(wn_weighted_pr_err, axis=1), axis=0)
    print(sc_weighted_pr_err.shape[0])
    for j in range(len(PAIR_RULE_GENES)):
        sig_pvals = 0
        for i in range(sc_weighted_pr_err.shape[0]):
            ttest, pval  = ttest_ind(sc_weighted_pr_err[i,:,j], wn_weighted_pr_err[i,:,j])
            sig_pvals += 1 if pval<=0.05 else 0
        print(sig_pvals/sc_weighted_pr_err.shape[0])
    plot_bar_comparison(mean_wn_err, std_wn_err, mean_sc_err, std_sc_err, PAIR_RULE_GENES, label_axes=True,
                        x_label='pair-rule gene', y_label='mean prediction error')
   #
   #  sc_weighted_pr_err_std = sc_pr_pred_err / (std_pr_gene_per_pos)
    wn_weighted_pr_err_std = wn_pr_pred_err / (std_pr_gene_per_pos)
    #plot_error_pr_error_per_position(sc_pr_pred_err, wn_pr_pred_err)
    #plot_error_pr_error_per_position(sc_weighted_pr_err_std, wn_weighted_pr_err_std)


def plot_error_pr_error_per_position(sc_weighted_pr_err_std, wn_weighted_pr_err_std):
    """
    This function plots the pair rule error per position
    """
    width=0.4
    for j in range(sc_weighted_pr_err_std.shape[2]):
        #plt.figure(figsize=(12, 10))
        bin_size = 50
        num_bins = sc_weighted_pr_err_std.shape[1] // bin_size
        positions = np.arange(num_bins)
        wn_agg = [
            wn_weighted_pr_err_std[:, i * bin_size:(i + 1) * bin_size, j].flatten()
            for i in range(num_bins)
        ]

        sc_agg = [
            sc_weighted_pr_err_std[:, i * bin_size:(i + 1) * bin_size, j].flatten()
            for i in range(num_bins)
        ]

        box1= plt.boxplot(wn_agg, positions=positions - width / 2, widths=width, patch_artist=True,
                    boxprops=dict(facecolor=DECODER_TYPE_COLOR['wn']),showfliers=False)
        box2 = plt.boxplot(sc_agg, positions=positions + width / 2, widths=width, patch_artist=True,
                           boxprops=dict(facecolor=DECODER_TYPE_COLOR['sc']), showfliers=False)

        plt.xlabel(POSITION_X_LABEL)
        plt.ylabel(f"{PAIR_RULE_GENES[j]} expression prediction error")
        plt.xticks(positions[::2], ((np.arange(num_bins))*bin_size/800)[::2])
        plt.legend([box1["boxes"][0], box2["boxes"][0]],
                   [DECODER_NAMES['wn'], DECODER_NAMES['sc']],
                   loc="upper right")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        #plt.ylim(0,30)
        # Show the plot
        plt.tight_layout()
        plt.show()

####################### mutant background embryo analyses related functions ##############################
def plot_all_mutant_results(mutant_type):
    """
    This function plots all of the mutant results shown in the paper
    """
    plot_mutant_gap_gene(mutant_type)
    plot_decoding_maps_mutants([mutant_type], DECODING_TYPES)
    ##binned over positions
    # reconstructions
    #aggregated results
    #plot_position_uncertainty_decoding_comparison_mutants(mutant_type)
    plot_decoding_map_uncertainty_mutants_binned_positions(mutant_type)

def plot_posterior_standard_deviation_comparison(wn_stds, sc_stds):
    """
    This function plots the posterior standard deviation comparison when neighborhood-informed
    compared to cell-independent
    """
    labels = ['neighborhood\n informed', 'cell\n independent']
    data = [wn_stds.flatten(), sc_stds.flatten()]  # Flatten to combine all samples

    print(f'wn mean: {np.mean(data[0])} sc mean: {np.mean(data[1])}')
    # t_stat, p_value = ttest_ind(wn_stds, sc_stds, equal_var=False)
    # plt.scatter(np.arange(len(p_value)), p_value)
    # plt.title('pvals')
    #
    # plt.show()

    t_stat, p_value = ttest_ind(wn_stds.flatten(), sc_stds.flatten(), equal_var=False)

    print(f't_stat: {t_stat}, p_value {p_value}')
    # Create the figure and axes
    fig, ax = plt.subplots()

    # Create the violin plot
    parts = ax.violinplot(data, showmeans=True, showmedians=True, widths=0.6)

    # Customize colors
    colors = [DECODER_TYPE_COLOR['wn'], DECODER_TYPE_COLOR['sc']]
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])  # Set the fill color for each violin
        pc.set_edgecolor('black')   # Set the edge color
        pc.set_alpha(0.7)           # Adjust transparency

    # Customize mean and median lines
    for key in ['cmeans', 'cmedians']:
        if key in parts:
            parts[key].set_color('black')  # Make mean/median lines black
            parts[key].set_linewidth(2)

    # Add labels and title
    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels)
    ax.set_ylabel('Posterior Standard Deviation')

    # Add a legend
    legend_patches = [
        Patch(facecolor=DECODER_TYPE_COLOR['wn'], edgecolor='black', label='Neighborhood Informed'),
        Patch(facecolor=DECODER_TYPE_COLOR['sc'], edgecolor='black', label='Cell Independent')
    ]
    ax.legend(handles=legend_patches, loc='upper center')

    plt.tight_layout()
    # plt.ylim(0,0.06)
    # Show the plot
    plt.show()

def plot_posterior_standard_deviation_comparison_wt():
    """
    This function plots the posterior distribution standard deviation comparison in wild type embryos
    """
    wn_map = TestResults.from_pickle(DROSO_RES_DIR, 'wn_wt', GAP_GENES)
    sc_map = TestResults.from_pickle(DROSO_RES_DIR, 'sc_wt', GAP_GENES)
    wn_stds = get_std_per_position_decoding_map(wn_map.normalized_decoding_map)
    sc_stds = get_std_per_position_decoding_map(sc_map.normalized_decoding_map)[:, 1:-1]

    plot_positional_information_in_bits_empirical_std(wn_stds, sc_stds)

    mean_wn_stds = np.mean(wn_stds, axis=0)
    std_wn_stds = np.std(wn_stds, axis=0)

    mean_sc_stds = np.mean(sc_stds, axis=0)
    std_sc_stds = np.std(sc_stds, axis=0)

    positions = np.linspace(POSITIONS_START, POSITIONS_END, len(mean_sc_stds))

    plt.plot(positions, mean_sc_stds, label=DECODER_NAMES['sc'], color=DECODER_TYPE_COLOR['sc'])
    plt.plot(positions, mean_wn_stds, label=DECODER_NAMES['wn'], color=DECODER_TYPE_COLOR['wn'])

    plt.fill_between(positions, mean_sc_stds - std_sc_stds, mean_sc_stds + std_sc_stds, color=DECODER_TYPE_COLOR['sc'],
                     alpha=0.5)
    plt.fill_between(positions, mean_wn_stds - std_wn_stds, mean_wn_stds + std_wn_stds, color=DECODER_TYPE_COLOR['wn'],
                     alpha=0.5)
    plt.xlabel(POSITION_X_LABEL)
    plt.ylabel('posterior standard deviation')
    plt.ylim(0,0.05)
    plt.tight_layout()
    plt.show()

def plot_posterior_standard_deviation_comparison_mutants(mutant_type):
    """
    This function plots the position posterior distribution standard deviation comparison in mutant embryos
    """
    wn_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_wn', GAP_GENES)
    sc_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_sc', GAP_GENES)
    wn_stds = get_std_per_position_decoding_map(wn_map.normalized_decoding_map)
    sc_stds = get_std_per_position_decoding_map(sc_map.normalized_decoding_map)[:, 1:-1]

    mean_wn_stds = np.mean(wn_stds, axis=0)
    std_wn_stds = np.std(wn_stds, axis=0)

    mean_sc_stds = np.mean(sc_stds, axis=0)
    std_sc_stds = np.std(sc_stds, axis=0)

    positions = np.linspace(POSITIONS_START, POSITIONS_END, len(mean_sc_stds))

    plt.plot(positions, mean_sc_stds, label=DECODER_NAMES['sc'], color=DECODER_TYPE_COLOR['sc'])
    plt.plot(positions, mean_wn_stds, label=DECODER_NAMES['wn'], color=DECODER_TYPE_COLOR['wn'])

    plt.fill_between(positions, mean_sc_stds - std_sc_stds, mean_sc_stds + std_sc_stds, color=DECODER_TYPE_COLOR['sc'],
                     alpha=0.5)
    plt.fill_between(positions, mean_wn_stds - std_wn_stds, mean_wn_stds + std_wn_stds, color=DECODER_TYPE_COLOR['wn'],
                     alpha=0.5)
    plt.xlabel(POSITION_X_LABEL)
    plt.ylabel('posterior standard deviation')
    plt.ylim(0,0.05)
    plt.tight_layout()
    plt.show()

def plot_decoding_map_uncertainty_mutants_binned_positions(mutant_type):
    """
    This function plots the decoding map standard deviation in mutants, binned across positions
    """
    wn_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_wn', GAP_GENES)
    sc_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_sc', GAP_GENES)
    wn_stds = get_std_per_position_decoding_map(wn_map.normalized_decoding_map)
    sc_stds= get_std_per_position_decoding_map(sc_map.normalized_decoding_map)[:,1:-1]
    #plot_posterior_standard_deviation_comparison(wn_stds, sc_stds)
    #plot_positional_information_in_bits()

    width = 0.4
    bin_size = 50
    num_bins = wn_stds.shape[1] // bin_size
    positions = np.arange(num_bins)
    wn_agg = [
        wn_stds[:,i * bin_size:(i + 1) * bin_size].flatten()
        for i in range(num_bins)
    ]

    sc_agg = [
        sc_stds[:, i * bin_size:(i + 1) * bin_size].flatten()
        for i in range(num_bins)
    ]

    box1 = plt.boxplot(wn_agg, positions=positions - width / 2, widths=width, patch_artist=True,
                       boxprops=dict(facecolor=DECODER_TYPE_COLOR['wn']), showfliers=False)
    box2 = plt.boxplot(sc_agg, positions=positions + width / 2, widths=width, patch_artist=True,
                       boxprops=dict(facecolor=DECODER_TYPE_COLOR['sc']), showfliers=False)

    plt.xlabel(POSITION_X_LABEL)
    plt.ylabel("posterior standard deviation")
    plt.xticks(positions[::2], (np.arange(num_bins)*bin_size/800)[::2])
    plt.legend([box1["boxes"][0], box2["boxes"][0]],
               [DECODER_NAMES['wn'], DECODER_NAMES['sc']],
               loc="upper right")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    #plt.ylim(0, 30)
    # Show the plot
    plt.tight_layout()
    plt.show()


def plot_wt_gap_gene():
    wt_gap_exp = load_all_wt_droso_train_data()[0]
    wt_gap_exp_arr = reformat_exp_data_to_arr(wt_gap_exp, EDGE_TRIM)
    min_vals = np.min(wt_gap_exp_arr, axis=1, keepdims=True)  # Shape: (num_samples, 1, num_features)
    max_vals = np.max(wt_gap_exp_arr, axis=1, keepdims=True)  # Shape: (num_samples, 1, num_features)

    wt_gap_exp_arr_normed = (wt_gap_exp_arr - min_vals) / (max_vals - min_vals)
    mean_pr_exp_per_gene = np.mean(wt_gap_exp_arr_normed, axis=0)
    std_pr_exp_per_gene = np.std(wt_gap_exp_arr_normed, axis=0)
    positions = np.linspace(.1, 0.9, wt_gap_exp_arr_normed.shape[1])

    for i, gene in enumerate(wt_gap_exp.columns):
        mean_per_pos_one_gene = mean_pr_exp_per_gene[:, i]
        std_per_pos_one_gene = std_pr_exp_per_gene[:, i]
        plt.plot(positions, mean_per_pos_one_gene, color=GAP_GENE_COLORS[gene])
        plt.fill_between(positions, mean_per_pos_one_gene - std_per_pos_one_gene,
                         mean_per_pos_one_gene + std_per_pos_one_gene, alpha=0.5, label=gene,
                         color=GAP_GENE_COLORS[gene])

    plt.legend(loc='upper right')
    plt.xlabel('position (x/L)')
    plt.ylabel('gap gene expression')
    plt.tight_layout()
    plt.show()



def plot_mutant_gap_gene(mutant_type):
    """
    This function plot the mutant gap gene expression, given the specific mutant type
    """
    mutant_gap_exp = get_mutant_gap_data(mutant_type)
    mutant_gap_exp_arr = reformat_exp_data_to_arr(mutant_gap_exp)[:, EDGE_TRIM:-EDGE_TRIM, :]
    min_vals = np.min(mutant_gap_exp_arr, axis=1, keepdims=True)  # Shape: (num_samples, 1, num_features)
    max_vals = np.max(mutant_gap_exp_arr, axis=1, keepdims=True)  # Shape: (num_samples, 1, num_features)

    mutant_pr_arr_normed = (mutant_gap_exp_arr - min_vals) / (max_vals - min_vals)
    mean_pr_exp_per_gene = np.mean(mutant_pr_arr_normed, axis=0)
    std_pr_exp_per_gene = np.std(mutant_pr_arr_normed, axis=0)
    positions = np.linspace(.1, 0.9, mutant_pr_arr_normed.shape[1])

    for i, gene in enumerate(mutant_gap_exp.columns):
        mean_per_pos_one_gene = mean_pr_exp_per_gene[:, i]
        std_per_pos_one_gene = std_pr_exp_per_gene[:, i]
        plt.plot(positions, mean_per_pos_one_gene, color=GAP_GENE_COLORS[gene])
        plt.fill_between(positions, mean_per_pos_one_gene - std_per_pos_one_gene,
                         mean_per_pos_one_gene + std_per_pos_one_gene, alpha=0.5, label=gene,
                         color=GAP_GENE_COLORS[gene])

    plt.legend(loc='upper right')
    plt.xlabel('position (x/L)')
    plt.ylabel('gap gene expression')
    plt.tight_layout()
    plt.show()

def plot_decoding_maps_mutants(mutant_types, decoding_types):
    """
    This function plots the decoding maps of the mutant types given the specific mutant type
    :param mutant_types:
    :param decoding_types:
    :return:
    """
    for mutant_type in mutant_types:
        print(mutant_type)
        for decoding_type in decoding_types:
            decoding_map_obj = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_{decoding_type}', GAP_GENES)
            #decoding_map_obj.normalize_decoding_map()
            decoding_map_obj.plot_decoding_map(' '.join(GAP_GENES), vmax=0.05, is_wt=False)


def plot_position_uncertainty_decoding_comparison_mutants(mutant_type):
    """
    This function plots the position standard deviation comparison between mutants given a specific mutant type
    """
    wn_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_wn', GAP_GENES)
    sc_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_sc', GAP_GENES)
    results = plot_prediction_loss_weighted_by_expression_variance_mutant(mutant_type, sc_map, wn_map, label_axes=True, x_label='pair rule gene', y_label='mean prediction error')
    return results

def plot_one_mutation_expression_correlation():
    """
    This function plots the mutant expression in correlation
    """
    mean_exp_correlations_wn = {'Eve':[], 'Run':[], 'Prd':[]}
    mean_exp_correlations_sc = {'Eve': [], 'Run': [], 'Prd': []}
    for mutant_type in ONE_MUTATION_TYPES:
        mutant_results = plot_position_uncertainty_decoding_comparison_mutants(mutant_type)
        for pair_rule_gene in PAIR_RULE_GENES:
            mean_exp_correlations_wn[pair_rule_gene].append(mutant_results[pair_rule_gene]['mean_Correlation_wn'])
            mean_exp_correlations_sc[pair_rule_gene].append(mutant_results[pair_rule_gene]['mean_Correlation_sc'])
    means_over_mutants_wn = [np.mean(mean_exp_correlations_wn[pr_gene]) for pr_gene in PAIR_RULE_GENES]
    stds_over_mutants_wn = [np.std(mean_exp_correlations_wn[pr_gene]) for pr_gene in PAIR_RULE_GENES]

    means_over_mutants_sc = [np.mean(mean_exp_correlations_sc[pr_gene]) for pr_gene in
                             PAIR_RULE_GENES]
    stds_over_mutants_sc = [np.std(mean_exp_correlations_sc[pr_gene]) for pr_gene in
                            PAIR_RULE_GENES]

    plot_bar_comparison(means_over_mutants_wn,stds_over_mutants_wn,means_over_mutants_sc, stds_over_mutants_sc, PAIR_RULE_GENES, mutant_type,
                        label_axes=True, x_label='pair-rule gene', y_label='expression correlations')


def plot_prediction_loss_weighted_by_expression_variance_mutant(mutant_type, sc_map, wn_map, label_axes=False, x_label=None, y_label=None):
    mean_wt_pr_exp = get_mean_wt_pr_per_position(EDGE_TRIM)
    mean_mutant_pr_exp = get_mutant_pr_data(mutant_type)
    full_mutant_pr_data = reformat_exp_data_to_arr(mean_mutant_pr_exp, EDGE_TRIM)
    mean_exp_mutant_pr_per_position = np.mean(full_mutant_pr_data, axis=0)
    mean_exp_mutant_per_positions_normed = ((mean_exp_mutant_pr_per_position - np.min(mean_exp_mutant_pr_per_position, axis=0))
                                        / (np.max(mean_exp_mutant_pr_per_position, axis=0) - np.min(
                mean_exp_mutant_pr_per_position, axis=0)))

    sc_pr_pred = sc_map.expected_pr_exp(mean_wt_pr_exp)
    wn_pr_pred = wn_map.expected_pr_exp(mean_wt_pr_exp[1:-1, :])

    sc_pr_pred_err = pr_prediction_err(sc_pr_pred, mean_exp_mutant_per_positions_normed)[:, 1:-1, :]
    wn_pr_pred_err = pr_prediction_err(wn_pr_pred, mean_exp_mutant_per_positions_normed[1:-1, :])

    wt_pr_exp = get_wt_pr_data()
    wt_pr_arr = reformat_exp_data_to_arr(wt_pr_exp, EDGE_TRIM)
    std_pr_gene_per_pos = np.std(wt_pr_arr, axis=0)[1:-1, :]

    sc_weighted_pr_err = (sc_pr_pred_err / std_pr_gene_per_pos)
    wn_weighted_pr_err = (wn_pr_pred_err / std_pr_gene_per_pos)

    plot_bar_comparison(np.mean(wn_weighted_pr_err, axis=(1, 0)),
                        np.std(np.mean(wn_weighted_pr_err, axis=1), axis=0),
                        np.mean(sc_weighted_pr_err, axis=(1, 0)),
                        np.std(np.mean(sc_weighted_pr_err, axis=1), axis=0), PAIR_RULE_GENES, mutant_type,
                        label_axes=label_axes, x_label=x_label, y_label=y_label)



def plot_mutant_pr_predictions_and_errors(mutant_type, plot_binned_position_errors=False):
    """
    This function plots the mutant pair rule prediction and errors given the specific mutant type
    """
    sc_pr_pred, wn_pr_pred = calculate_mutant_pr_pred(mutant_type)

    mean_mutant_pr_exp = get_mutant_pr_data(mutant_type)
    full_mutant_pr_data = reformat_exp_data_to_arr(mean_mutant_pr_exp, EDGE_TRIM)
    mean_exp_mutant_pr_per_position = np.mean(full_mutant_pr_data, axis=0)

    sc_weighted_pr_err, wn_weighted_pr_err = calculate_mutant_pr_pred_err(mutant_type, sc_pr_pred, wn_pr_pred, mean_exp_mutant_pr_per_position)
    plot_predicted_prs_across_positions(wn_pr_pred, sc_pr_pred[:, 1:-1, :], full_mutant_pr_data[:, 1:-1, :])
    if plot_binned_position_errors:
        plot_error_pr_error_per_position(sc_weighted_pr_err, wn_weighted_pr_err)

def plot_mutant_pr_exp(mutant_type='osk'):
    """
    This function plots the mutant pair rule gene expression give the specific mutant type
    """
    mutant_pr = get_mutant_pr_data(mutant_type)
    mutant_pr_arr = reformat_exp_data_to_arr(mutant_pr)[:, EDGE_TRIM:-EDGE_TRIM, :]
    plot_pr_exp(mutant_pr_arr)

def plot_wt_pr_exp():
    """
    This function plots the wild type pair rule gene expression profiles
    """
    wt_pr = get_wt_pr_data()
    wt_pr_arr = reformat_exp_data_to_arr(wt_pr, EDGE_TRIM)
    plot_pr_exp(wt_pr_arr)

def plot_pr_exp(pr_arr):
    """
    This function plots pair rule gene expression given the pair rule array
    """
    min_vals = np.min(pr_arr, axis=1, keepdims=True)  # Shape: (num_samples, 1, num_features)
    max_vals = np.max(pr_arr, axis=1, keepdims=True)  # Shape: (num_samples, 1, num_features)

    mutant_pr_arr_normed = (pr_arr - min_vals) / (max_vals - min_vals)
    mean_pr_exp_per_gene = np.mean(mutant_pr_arr_normed, axis=0)
    std_pr_exp_per_gene = np.std(mutant_pr_arr_normed, axis=0)
    positions = np.linspace(POSITIONS_START, POSITIONS_END, std_pr_exp_per_gene.shape[0])
    for i, gene in enumerate(PAIR_RULE_GENES):
        mean_per_pos_one_gene = mean_pr_exp_per_gene[:, i]
        std_per_pos_one_gene = std_pr_exp_per_gene[:, i]
        plt.plot(positions, mean_per_pos_one_gene, color=PAIR_RULE_COLORS[gene])
        plt.fill_between(positions, mean_per_pos_one_gene - std_per_pos_one_gene,
                         mean_per_pos_one_gene + std_per_pos_one_gene, alpha=0.5, label=gene,
                         color=PAIR_RULE_COLORS[gene])
        plt.xlabel(POSITION_X_LABEL)
        plt.ylabel(EXP_Y_LABEL)
        plt.tight_layout()
        plt.legend()
        plt.show()