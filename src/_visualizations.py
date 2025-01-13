"""
Figure 1 visuals
"""
from src._imports import *
from src._constants import *
from data._preprocessing import *
from data._data import *
from src._utils import *
#### panels b and c, gap gene and pair rule genes pairwise correlations:
def gene_expression_pairwise_correlation_plots():
    plot_gap_gene_pairwise_correlations()
    plot_pair_rule_gene_pairwise_correlations()


def plot_pair_rule_gene_pairwise_correlations():
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

    plt.legend(fontsize=20, loc='upper center')
    plt.xlabel('pair-rule genes pairwise correlation')
    plt.ylabel('relative density')
    plt.tight_layout()
    plt.ylim(0, .1)
    plt.show()
    pass

def plot_gap_gene_pairwise_correlations():
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
    plt.legend(fontsize=20, loc='upper center')
    plt.xlabel('gap-gene pairwise correlation')
    plt.ylabel('relative density')
    plt.tight_layout()
    plt.ylim(0, .1)
    plt.show()
def pairwise_correlations_gap_genes():
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
def plot_positional_information_in_bits():
    sigma_x_wn = caclulate_positional_error_per_decoding_map_MAP_positions(f'wn_wt_sigmax')  # vs 'wn'
    sigma_x_sc = caclulate_positional_error_per_decoding_map_MAP_positions(f'sc_wt_sigmax')
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
    plt.xlabel('position (x/L)')
    plt.ylabel('positional information in bits')
    plt.ylim(0,10)
    plt.xlim(0.1,0.9)
    plt.tight_layout()
    plt.show()

def plot_decoding_with_one_gene_comparison(gene):
    pass




