from src._figures import *
from src._utils import *
from src._figures import *
from src._create_datasets_droso import *
from data.Gastruloids import *
from src._create_datasets_gastruloids import *
from src._create_datasets_neural_tube import *

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
    # with open(os.path.join(NEURAL_TUBE_WT_PATH,'expressions_h=5.pkl'), 'rb') as f:
    #     nt_wt_5 = pickle.load(f)
    plot_summarized_neural_tube_over_axis_over_timepoints(NEURAL_TUBE_SET_A_GENES)

    #calculate_neural_tube_decoding_maps_one_gene_group(NEURAL_TUBE_SET_A_GENES, data_path=nt_path)

    # with open(BRA_10_PATH, 'rb') as f:
    #     bra_10dict = pickle.load(f)
    # format_gastru_like_droso(bra_10dict)
    # genes = ['Cdx2', 'Sox2']
    #calculate_gastru_decoding_maps_one_gene_group(genes, data_path=CDX2_RES_PATH)
    # print('wn')
    # plot_decoding_maps([genes], ['wn_gastru'], VMAXS_THREE_GENES, xlim=False, results_dir=GASTRU_RES_DIR)
    # print('sc')
    # plot_decoding_maps([genes], ['sc_gastru'], VMAXS_THREE_GENES, xlim=False,
    #                    results_dir=GASTRU_RES_DIR)
    #calculate_gastru_decoding_maps_one_gene_group(['Cdx2','Sox2'], data_path=CDX2_RES_PATH)
    #plot_decoding_maps([['Cdx2','Sox2']], ['sc_gastru', 'wn_gastru'], VMAXS_THREE_GENES, xlim=False, results_dir=GASTRU_RES_DIR)
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


