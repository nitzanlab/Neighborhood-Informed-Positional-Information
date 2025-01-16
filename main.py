from src._figures import *
from src._utils import *
from src._figures import *
from src._create_datasets import *

def calculate_and_save_all_decoding_maps():
    ##define directory paths in src._constants
    ##
    calculate_all_gene_subset_decoding_maps_WT()
    calculate_all_mutant_decoding_maps()



def plot_all_figures():
    #plot_figure1_panels()
    #plot_figure2_and_related_supp_panels()
    #plot_figure3_panels()
    #plot_figure5_panels()
    plot_figure6_panels()


if __name__ == '__main__':
    set_style()

    #calculate_all_mutant_decoding_maps()
    #plot_all_figures()
    #plot_posterior_standard_deviation_comparison_wt()
    #plot_all_decoding_maps_comparison_wt()
    #calculate_all_gene_subset_decoding_maps_WT()
    #calculate_all_mutant_decoding_maps()
    plot_all_figures()
    # for mutant_type in MUTANT_TYPES:
    #     plot_decoding_map_uncertainty_mutants_binned_positions(mutant_type)
    #     print(mutant_type)
    #     plot_mutant_pr_predictions_and_errors(mutant_type, plot_binned_position_errors=False)


