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
    plot_figure1_panels()
    plot_figure2_and_related_supp_panels()
    plot_figure3_panels()
    plot_figure5_panels()
    plot_figure6_panels()

def reproduce_all_results():
    calculate_and_save_all_decoding_maps()
    plot_all_figures()

if __name__ == '__main__':
    set_style()
    reproduce_all_results()

