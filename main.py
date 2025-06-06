from src._figures import *
from src._utils import *
from src._figures import *
from src._create_datasets import *
from data.Gastruloids import *

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
    with open(BRA_RES_PATH, 'rb') as f:
        bra_dict = pickle.load(f)
    with open(CDX2_RES_PATH, 'rb') as f:
        cdx2_dict = pickle.load(f)
    with open(FOXC1_RES_PATH, 'rb') as f:
        foxc1_dict = pickle.load(f)
    plot_gastruloids_data(bra_dict,'BRA')
    plot_gastruloids_data(cdx2_dict, 'CDX2')
    plot_gastruloids_data(foxc1_dict, 'FOXC1')
    #set_style()
    #reproduce_all_results()


