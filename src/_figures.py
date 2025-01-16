from src._visualizations import *
def plot_figure1_panels():
    """
    This figure presents the positional information encoded in neighboring positions in wild type Drosophila embryos
    This functions plots panels b,c,d,e,f,h, and i. Panels a and g are schematic diagrams drawn in biorender.com
    """
    ##panels b and c: gene expression correlations
    gene_expression_pairwise_correlation_plots()

    #panel d: positional information in bits
    plot_positional_information_in_bits()

    #panel e: decoding neighborhood-informed and cell independent given Kr expression
    plot_decoding_maps_comparison_wt(ONE_GENE_EXAMPLE_GENE, VMAXS_ONE_GENE, xlim=True)

    #panel f: decoding neighborhood-informed and cell independent give Kr, Gt, and Hb expression
    plot_decoding_maps_comparison_wt(THREE_GENES_EXAMPLE_GENES, VMAXS_THREE_GENES, xlim=False)

    #aggregated results - mean MAP error over positions and over wild type embryos, mean per number of genes used
    #for decoding
    #panel h:
    plot_position_MAP_error_all_gene_subsets()
    #
    # #same aggregation , but looking at the standard deviation of the position posterior distribution
    # #neighborhood-informed vs cell-independent
    # #panel i:
    plot_position_posterior_std_all_gene_subsets()

def plot_figure2_and_related_supp_panels():
    """
    This figure presents pair-rule expression profile prediction given the position probabilistic decoder.
    We show both qualitatively, quantitivately across the AP axis, and aggregated summary results
    comparing neighborhood-informed decoding based vs cell-independent decoding based pair-rule expression profile predictions
    This function plots panels b,c,and d. Panel a is schematic and was drawn in biorender.com
    """
    ##pair rule prediction across AP axis - qualitative comparison
    #figure 2 panel b - Prd prediction and supplementary Run and Eve prediction
    pair_rule_prediction_across_positions()

    #figure 2 panel c -reconstruction error across AP axis, d - summary results - mean prediction error per pair-rule gene:
    plot_all_pr_genes_mean_prediction_error()

def plot_figure3_panels():
    """
    This figure presents position decoding and pair-rule expression prediction results
    for Osk mutant background embryos
    :return:
    """
    plot_all_mutant_results(mutant_type='osk')

    ##TODO:wild type expression of all gap genes for comparison

    ##mutant pair-rule prediction and errors:
    #plot_mutant_pr_predictions_and_errors(mutant_type='osk', plot_binned_position_errors=True)

    #plot_one_mutation_expression_correlation()




def plot_figure5_panels():
    for mutant_type in MUTANT_TYPES:
        plot_all_mutant_results(mutant_type)


def plot_figure6_panels():
    for mutant_type in MUTANT_TYPES:
        print(mutant_type)
        fig, axes = plt.subplots(1, 6, figsize=(
        24, 5))  # 1 row, 6 columns, wider width (sharex=False and sharey=False by default)

        handles = []  # Store handles for the legend
        labels = []
        mean_wt_pr_exp = get_mean_wt_pr_per_position(EDGE_TRIM)
        mutant_pr = get_mutant_pr_data(mutant_type)
        mutant_pr_arr = reformat_exp_data_to_arr(mutant_pr)[:, EDGE_TRIM:-EDGE_TRIM, :]
        min_vals = np.min(mutant_pr_arr, axis=1, keepdims=True)
        max_vals = np.max(mutant_pr_arr, axis=1, keepdims=True)

        # Normalize the data
        mutant_pr_arr_normed = (mutant_pr_arr - min_vals) / (max_vals - min_vals)
        mean_pr_exp_per_gene = np.mean(mutant_pr_arr_normed, axis=0)[1:-1, :]
        std_pr_exp_per_gene = np.std(mutant_pr_arr_normed, axis=0)[1:-1, :]

        wn_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_wn', GAP_GENES)
        sc_map = TestResults.from_pickle(DROSO_RES_DIR, f'{mutant_type}_sc', GAP_GENES)
        sc_pr_pred = sc_map.expected_pr_exp(mean_wt_pr_exp)[:, 1:-1, :]
        wn_pr_pred = wn_map.expected_pr_exp(mean_wt_pr_exp[1:-1, :])
        positions = np.linspace(0.1, 0.9, sc_pr_pred.shape[1])

        mean_wn_pred = np.mean(wn_pr_pred, axis=0)
        std_wn_pred = np.std(wn_pr_pred, axis=0)

        mean_sc_pred = np.mean(sc_pr_pred, axis=0)
        std_sc_pred = np.std(sc_pr_pred, axis=0)

        for i, gene in enumerate(PAIR_RULE_GENES[:3]):
            # WN plot
            ax_wn = axes[2 * i]
            ax_wn.fill_between(positions, mean_pr_exp_per_gene[:, i] - std_pr_exp_per_gene[:, i],
                               mean_pr_exp_per_gene[:, i] + std_pr_exp_per_gene[:, i], color='grey', alpha=0.3)
            ax_wn.fill_between(positions, mean_wn_pred[:, i] - std_wn_pred[:, i],
                               mean_wn_pred[:, i] + std_wn_pred[:, i], color=DECODER_TYPE_COLOR['wn'], alpha=0.3)
            line_wn = ax_wn.plot(positions, mean_wn_pred[:, i], color=DECODER_TYPE_COLOR['wn'])[0]
            line_pr = ax_wn.plot(positions, mean_pr_exp_per_gene[:, i], color='grey')[0]
            if i == 0:
                ax_wn.set_ylabel('Gene Expression')
                ax_wn.set_yticks(np.arange(0, 1.05, .2))  # Set yticks for leftmost plot
            else:
                ax_wn.set_yticks([])  # Remove yticks for the other plots

            # SC plot
            ax_sc = axes[2 * i + 1]
            ax_sc.fill_between(positions, mean_pr_exp_per_gene[:, i] - std_pr_exp_per_gene[:, i],
                               mean_pr_exp_per_gene[:, i] + std_pr_exp_per_gene[:, i], color='grey', alpha=0.3)
            ax_sc.fill_between(positions, mean_sc_pred[:, i] - std_sc_pred[:, i],
                               mean_sc_pred[:, i] + std_sc_pred[:, i], color=DECODER_TYPE_COLOR['sc'], alpha=0.3)
            line_sc = ax_sc.plot(positions, mean_sc_pred[:, i], color=DECODER_TYPE_COLOR['sc'])[0]
            ax_sc.plot(positions, mean_pr_exp_per_gene[:, i], color='grey', linewidth=2)  # Add grey line here
            ax_sc.set_yticks([])  # Remove yticks for SC plot

            # Add gene name in the bottom-right corner of each subplot (larger and bold)
            ax_sc.text(0.95, 0.05, f'{gene}', transform=ax_sc.transAxes, ha='right', va='bottom', fontsize=22,
                       fontweight='bold')
            ax_wn.text(0.95, 0.05, f'{gene}', transform=ax_wn.transAxes, ha='right', va='bottom', fontsize=22,
                       fontweight='bold')

            # Collect handles and labels for the legend
            if i == 0:
                handles.extend([line_pr, line_wn, line_sc])
                labels.extend([f'{gene} mean expression', DECODER_NAMES['wn'], DECODER_NAMES['sc']])

        # Add a single legend for all subplots outside of the plots
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.03), ncol=3)

        for ax in axes:
            ax.set_xlabel('Position (x/L)')

        # Adjust layout to make more space and remove excess white space
        plt.subplots_adjust(top=0.85, bottom=0.2, left=0.05, right=0.97, hspace=0.05, wspace=0.05)
        plt.show()
    # for mutant_type in MUTANT_TYPES:
    #     plot_mutant_pr_predictions_and_errors(mutant_type)

