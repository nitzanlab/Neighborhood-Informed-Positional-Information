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
    set_style()
    plot_positional_information_neural_tube()
    #plot_summarized_neural_tube_gene_combos_one_timepoint(NEURAL_TUBE_SET_A_GENES,'35')
    #plot_positional_information_gastruloids()
    #plot_summarized_neural_tube_over_axis_over_timepoints(NEURAL_TUBE_SET_A_GENES)
    # neural_tube_data = NeuralTube(NEURAL_TUBE_WT_PATH,'expressions_h=5.pkl')
    # #neural_tube_data.plot_gene_exp_over_positions(NEURAL_TUBE_SET_A_GENES)
    # neural_tube_data.plot_comparison_position_inf_GT(NEURAL_TUBE_SET_A_GENES, 'wt')
    #print('')
    #get_all_subsets_pos_error(to_plot=True)
    #get_pos_error_three_genes_with_sox('Bra','Foxc1')
    #create_cov_and_mean_one_gene_wn('Bra')
    #plot_all_gastruloid_plots()
    #create_cov_and_mean_one_gene_wn('Bra')
    #create_cov_and_mean_one_gene_wn('Foxc1')
    #create_cov_and_mean_one_gene_wn('Cdx2')
    # compare_position_error_all_datasets()
    # covs_all_gastru_wn, mean_all_gastru_wn = create_cov_and_mean_joint_datasets_wn()
    # covs_gastru_sc, means_gastru_sc = create_covariance_sc_joint_datasets()
    #
    # # with open(os.path.join(NEURAL_TUBE_WT_PATH,'expressions_h=5.pkl'), 'rb') as f:
    # #     nt_wt_5 = pickle.load(f)
    #plot_summarized_neural_tube_over_axis_over_timepoints(NEURAL_TUBE_SET_A_GENES)

    #calculate_neural_tube_decoding_maps_one_gene_group(NEURAL_TUBE_SET_A_GENES, data_path=nt_path)

    # with open(BRA_10_PATH, 'rb') as f:
    #     bra_10dict = pickle.load(f)
    # format_gastru_like_droso(bra_10dict)
    # genes = ['Cdx2', 'Sox2']
    # calculate_gastru_decoding_maps_one_gene_group(genes, data_path=CDX2_RES_PATH)
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
    # #plot_gastruloids_data(bra_dict,'BRA')
    # plot_gastruloids_data(foxc1_dict, 'foxc1')
    #plot_gastruloids_data(foxc1_dict, 'FOXC1')
    #set_style()
    #reproduce_all_results()


