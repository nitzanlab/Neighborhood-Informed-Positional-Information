from src._constants import *
from test_results_analysis.TestResults import *

def set_style():
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({'figure.titlesize': 26, 'figure.titleweight': 'bold', 'axes.titlesize': 22,
                         'axes.titleweight': "bold", 'axes.labelsize': 24,'axes.labelweight': 'bold',
                         "ytick.labelsize": 28, "xtick.labelsize": 28, 'legend.fontsize': 24,'font.family': 'Candara,  Arial'})
    plt.rcParams.update({'figure.figsize': (8, 6)})
    plt.rcParams.update({'savefig.dpi': 300})


def caclulate_positional_error_per_decoding_map_MAP_positions(decoding_type):
    wt_decoding_res = TestResults.from_pickle(DROSO_RES_DIR, decoding_type, GAP_GENES)
    wt_cov = wt_decoding_res.train_std
    decoding_map = wt_decoding_res.normalized_decoding_map
    MAP_positions = np.argmax(decoding_map, axis=1)
    mean_exp = wt_decoding_res.train_mean

    mean_exp_slopes = np.diff(mean_exp,axis=0)
    mean_exp_extended = np.vstack([mean_exp_slopes, mean_exp_slopes[-1]])

    num_genes = mean_exp_slopes.shape[1]
    num_embryos, num_positions = MAP_positions.shape[0], MAP_positions.shape[1] - 1
    MAP_pos_mean_exp_slopes = np.zeros((num_embryos, num_positions,num_genes))

    # Iterate over embryos to extract slopes corresponding to MAP positions
    for embryo in range(num_embryos):
        MAP_pos_mean_exp_slopes[embryo] = mean_exp_extended[MAP_positions[embryo, 1:]]

    wt_covs_MAP = wt_cov[MAP_positions][:,1:,:,:]
    num_positions = MAP_pos_mean_exp_slopes.shape[1]
    num_embryos = MAP_pos_mean_exp_slopes.shape[0]
    position_error = np.zeros((num_embryos,num_positions))
    for embryo in range(num_embryos):
        for pos in range(num_positions):
            position_error[embryo,pos] = MAP_pos_mean_exp_slopes[embryo,pos,:]@np.linalg.inv(wt_covs_MAP[embryo,pos,:,:])@MAP_pos_mean_exp_slopes[embryo,pos,:]
    pos_err = np.sqrt((1/(position_error)))/num_positions
    return pos_err

def calculate_positional_inf(sigma_x):
    return np.log2(1/(sigma_x*np.sqrt(2*np.pi*np.e)))