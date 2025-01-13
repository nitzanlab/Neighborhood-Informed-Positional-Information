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
    plot_decoding_maps_comparison(ONE_GENE_EXAMPLE_GENE, VMAXS_ONE_GENE, xlim=True)

    #panel f: decoding neighborhood-informed and cell independent give Kr, Gt, and Hb expression
    plot_decoding_maps_comparison(THREE_GENES_EXAMPLE_GENES, VMAXS_THREE_GENES, xlim=False)

    #aggregated results - mean MAP error over positions and over wild type embryos, mean per number of genes used
    #for decoding
    #panel h:
    plot_position_MAP_error_all_gene_subsets()

    #same aggregation , but looking at the standard deviation of the position posterior distribution
    #neighborhood-informed vs cell-independent
    #panel i:
    plot_position_posterior_std_all_gene_subsets()

def plot_figure2_and_supp_related_panels():
    """
    This figure presents pair-rule expression profile prediction given the position probabilistic decoder.
    We show both qualitatively, quantitivately across the AP axis, and aggregated summary results
    comparing neighborhood-informed decoding based vs cell-independent decoding based pair-rule expression profile predictions
    This function plots panels b,c,and d. Panel a is schematic and was drawn in biorender.com
    """
    ##pair rule prediction across AP axis - qualitative comparison
    #figure 2 panel b - Prd prediction and supplementary Run and Eve prediction
    #pair_rule_prediction_across_positions()

    #figure 2 panel c - reconstruction error across AP axis:

    #figure 2 panel d - summary results - mean prediction error per pair-rule gene:
    plot_all_pr_genes_mean_prediction_error()

def plot_figure3_and_supp_related_panels():
    pass

