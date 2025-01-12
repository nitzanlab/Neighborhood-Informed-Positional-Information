"""
Figure 1 visuals
"""
from src._imports import *
from src._constants import *
from data._preprocessing import *
from data._data import *
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
