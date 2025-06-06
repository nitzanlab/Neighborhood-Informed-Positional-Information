from src._figures import *
from src._utils import *
from src._figures import *
from src._create_datasets_droso import *
from data.Gastruloids import *
from src._create_datasets_gastruloids import *

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
    with open(BRA_10_PATH, 'rb') as f:
        bra_10dict = pickle.load(f)
    #format_gastru_like_droso(bra_10dict)
    calculate_gastru_decoding_maps_one_gene_group(['Bra','Sox2'], data_path=BRA_10_PATH)
    #plot_decoding_maps([['Bra','Sox2']], ['sc_gastru'], VMAXS_THREE_GENES, xlim=False, results_dir=GASTRU_RES_DIR)
    #calculate_all_gene_subset_decoding_maps_WT()

    # with open(CDX2_RES_PATH, 'rb') as f:
    #     cdx2_dict = pickle.load(f)
    # with open(FOXC1_RES_PATH, 'rb') as f:
    #     foxc1_dict = pickle.load(f)
    #plot_gastruloids_data(bra_dict,'BRA')
    #plot_gastruloids_data(bra_10dict, 'BRA_10')
    #plot_gastruloids_data(foxc1_dict, 'FOXC1')
    #set_style()
    #reproduce_all_results()


