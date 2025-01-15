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
    plot_figure3_panels()
    #plot_figure5_panels()
    #plot_figure6_panels()


if __name__ == '__main__':
    set_style()
    plot_all_figures()
    # for mutant_type in ONE_MUTATION_TYPES:
    #     print(mutant_type)
    #     plot_mutant_pr_predictions_and_errors(mutant_type, plot_binned_position_errors=True)
    #plot_one_mutation_expression_correlation()
    # for mutant_type in MUTANT_TYPES:
    #     print(mutant_type)
    #     plot_decoding_map_uncertainty_mutants_binned_positions(mutant_type)
    #     #plot_all_mutant_results(mutant_type)
