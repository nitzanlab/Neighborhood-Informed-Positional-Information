import matplotlib.pyplot as plt
import numpy as np

from src._imports import *


from src._constants import *
from data._preprocessing import *
from data.droso_data import *
from data.Data import *
from data.gastruloid_data import *

class Gastruloids(Data):
    def __init__(self, data_path, data=None, training=False, save_training=False, save_dir=None, load_dir=None, edge_trim=None):
        self.data_path = data_path
        self.meta_data = None  # includes orient, dist, age, genotype,..
        self.save_training = save_training
        self.save_dir = save_dir
        self.edge_trim = edge_trim
        super().__init__(data, 'Gastruloid')
        if training:
            self.preprocess(data, edge_trim)

        else:
            self.means_sc = None
            self.covs_wn = None
            self.means_wn = None
            self.covs_wn = None


    def preprocess(self, data=None, edge_trim=None): #preprcoess training data
        print("Preprocessing Gastruloid data")
        all_training_data = format_gastru_like_droso(load_gastruloid_data(self.data_path))
        self.meta_data = ''
        self.define_data_structures(all_training_data)


    def define_data_structures(self, normalized_data):
        gene_exp_data = normalized_data
        self.genes = {gene: i for i, gene in enumerate(gene_exp_data.columns)}
        training_arr = reshape_gene_data_to_arr(gene_exp_data, self.genes)
        #TODO handle nans
        self.train_data = np.nan_to_num(training_arr, nan=0.0)
        if self.edge_trim is not None:
            self.train_data = self.train_data[:,self.edge_trim:-self.edge_trim,:]


    def train_wn(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        train_data_sbst_genes = self.train_data[:,:,decoding_genes_idx]
        train_wn_data = self.reshape_data_for_wn(train_data_sbst_genes)
        self.learn_mean_wn(train_wn_data,decoding_genes_idx)
        self.learn_covariance_wn(train_wn_data, decoding_genes_idx)
        if self.save_training:
            self.save_dir('wn')

    def learn_mean_sc(self, decoding_genes_idx=np.arange(len(GAP_GENES))):
        self.means_sc = np.mean(self.train_data[:, : , decoding_genes_idx], axis=0)

    def learn_covariance_sc(self, decoding_genes_idx=np.arange(len(GAP_GENES))):
        training_data = self.train_data[:,:,decoding_genes_idx]
        self.std_sc = get_cov(training_data)

    def test_wn(self, test_data, decoding_genes):
        processed_test_data = self.prepare_test_data(test_data, decoding_genes)
        reshaped_test_wn = self.reshape_data_for_wn(processed_test_data)
        decoding_map = self.get_position_distribution(reshaped_test_wn, self.means_wn, self.covs_wn)
        return decoding_map

    def test_sc(self, test_data, decoding_genes):
        processed_test_data = self.prepare_test_data(test_data, decoding_genes)
        decoding_map = self.get_position_distribution(processed_test_data, self.means_sc, self.std_sc)
        return decoding_map
    def prepare_test_data(self, test_data, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        processed_test_data = reshape_gene_data_to_arr(test_data, decoding_genes)[:, :, decoding_genes_idx]
        if self.edge_trim is not None:
            processed_test_data = processed_test_data[:, self.edge_trim:-self.edge_trim, :]
        return processed_test_data

    def train_and_test_sc(self, test_data,decoding_genes):
        self.train_sc(decoding_genes)
        decoding_map = self.test_sc(test_data, decoding_genes)
        return decoding_map

    def train_and_test_wn(self, test_data, decoding_genes):
        self.train_wn(decoding_genes)
        decoding_map = self.test_wn(test_data, decoding_genes)
        return decoding_map


    def learn_mean_wn(self, train_data_wn=None):
        reshaped_wn_position_means = np.concatenate(
            (self.means_sc[:-2, :], self.means_sc[1:-1, :], self.means_sc[2:, :]), axis=1)
        self.means_wn = reshaped_wn_position_means

    def learn_covariance_wn(self, train_data_wn=None):
        #TODO depending on which genes were measured together
        self.covs_wn = get_cov(train_data_wn)

    def reshape_data_for_wn(self, data):
        reshaped_wn_data = np.concatenate(
            (data[:, :-2, :], data[:, 1:-1, :], data[:, 2:, :]),axis=2)
        return reshaped_wn_data

    def train_sc(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        self.learn_mean_sc(decoding_genes_idx)
        self.learn_covariance_sc(decoding_genes_idx)
        if self.save_training:
            self.save_dir('sc')

    def train_wn(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        train_data_wn = self.reshape_data_for_wn(self.train_data[:,:, decoding_genes_idx])
        self.learn_mean_wn(train_data_wn)#, decoding_genes_idx)
        self.learn_covariance_wn(train_data_wn)#, decoding_genes_idx)

    def plot_gene_exp_over_positions(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        data_gene_subset = self.train_data[:,:,decoding_genes_idx]
        mean_gene_exp_over_positions = np.mean(data_gene_subset, axis=0)[EDGE_TRIM:-EDGE_TRIM,:]
        std_gene_exp_over_positions = np.std(data_gene_subset, axis=0)[EDGE_TRIM:-EDGE_TRIM,:]
        positions = np.linspace(POSITIONS_START, POSITIONS_END, mean_gene_exp_over_positions.shape[0])
        for i, gene in enumerate(decoding_genes):
            mean_per_pos_one_gene = mean_gene_exp_over_positions[:,i]
            std_per_pos_one_gene = std_gene_exp_over_positions[:, i]
            plt.plot(positions, mean_per_pos_one_gene, color=GAP_GENE_COLORS[gene])
            plt.fill_between(positions, mean_per_pos_one_gene-std_per_pos_one_gene, mean_per_pos_one_gene+std_per_pos_one_gene,alpha=0.5, label=gene,  color=GAP_GENE_COLORS[gene])
        plt.xlabel(POSITION_X_LABEL)
        plt.ylabel(EXP_Y_LABEL)
        plt.legend()
        plt.tight_layout()
        plt.show()

    def calculate_positional_error_per_decoding_map_GT_positions(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        self.learn_mean_sc(decoding_genes_idx)
        self.learn_covariance_sc(decoding_genes_idx)
        self.learn_mean_wn()
        self.learn_covariance_wn()

    def calculate_position_inf_GT(self, decoding_type):
        if decoding_type == "sc":
            mean_exp = self.means_sc[1:-1, :]
            covs = self.std_sc[1:-1, :, :]
        elif decoding_type == "wn":
            mean_exp = self.means_wn
            covs = self.covs_wn
        else:
            print("Unknown decoding")
            return
        num_genes = mean_exp.shape[1]
        num_pos = mean_exp.shape[0]
        mean_exp_slopes = np.diff(mean_exp, axis=0)
        mean_exp_slopes = np.vstack([mean_exp_slopes, mean_exp_slopes[-1]])
        position_error = np.zeros(num_pos)
        for pos in range(num_pos):
            position_error[pos] = 1 / (mean_exp_slopes[pos, :] @ np.linalg.inv(
                covs[pos, :, :]) @ mean_exp_slopes[pos, :])
        # TODO add calculation and plot
        return position_error

    def plot_comparison_position_inf_GT(self, genes):
        self.calculate_positional_error_per_decoding_map_GT_positions(genes)
        position_error_sc = self.calculate_position_inf_GT('sc')
        position_error_wn = self.calculate_position_inf_GT('wn')
        plt.plot(np.linspace(0, 1, len(position_error_sc)), position_error_sc, label='sc')
        plt.plot(np.linspace(0, 1, len(position_error_wn)), position_error_wn, label='wn')
        plt.legend()
        plt.title('position information ground truth positions Neural Tube')
        plt.ylim(0, 100)
        plt.show()

    # def calculate_positional_error_per_decoding_map_GT_positions(self, decoding_genes):
    #     decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
    #     self.learn_mean_sc(decoding_genes_idx)
    #     self.learn_covariance_sc(decoding_genes_idx)
    #     self.learn_mean_wn()
    #     self.learn_covariance_wn()
def get_cov(training_data):
    num_positions = training_data.shape[1]
    num_features = training_data.shape[2]
    covs = np.zeros((num_positions, num_features, num_features))
    for pos in range(num_positions):
        covs[pos] = np.cov(training_data[:, pos, :], rowvar=False)
    return covs

def reshape_gene_data_to_arr(gene_exp_data, genes):
    reshaped_gene_data = []
    for gene in genes:
        reshaped_gene_data.append(np.vstack(gene_exp_data[gene]))
    reshaped_data = np.dstack(np.array(reshaped_gene_data))
    return reshaped_data

def get_sox2_exp():
    sox_all_exp = []
    for gene_name in gene_data_path_dict.keys():
        with open(gene_data_path_dict[gene_name], 'rb') as f:
            exp_dict = pickle.load(f)
        sox_exp = np.array(exp_dict['Sox2'])
        sox_exp[np.isnan(sox_exp)] = 0
        sox_all_exp.append(sox_exp)
    return np.concatenate(sox_all_exp)

def get_one_gene_exp_over_AP_axis(gene_name:str):
    with open(gene_data_path_dict[gene_name], 'rb') as f:
        exp_dict = pickle.load(f)
    exp_arr = np.array(exp_dict[gene_name])
    exp_arr[np.isnan(exp_arr)] = 0
    return exp_arr

def get_wn_exp_from_sc_exp(sc_exp_arr):
    wn_exp_arr = sliding_window_view(sc_exp_arr, window_shape=3, axis=1)
    return wn_exp_arr
def get_wn_exp(gene_name:str):
    sc_exp_arr = get_one_gene_exp_over_AP_axis(gene_name)
    wn_exp_arr = get_wn_exp_from_sc_exp(sc_exp_arr)
    return wn_exp_arr
def get_sox2_wn_exp():
    sc_exp_arr = get_sox2_exp()
    wn_exp_arr = get_wn_exp_from_sc_exp(sc_exp_arr)
    return wn_exp_arr

def get_sox2_wn_exp_in_one_dataset(second_gene_name:str):
    with open(gene_data_path_dict[second_gene_name], 'rb') as f:
        exp_dict = pickle.load(f)
    exp_arr = np.array(exp_dict['Sox2'])
    exp_arr[np.isnan(exp_arr)] = 0
    sox2_wn_exp = get_wn_exp_from_sc_exp(exp_arr)
    return sox2_wn_exp

def get_sox2_sc_exp_in_one_dataset(second_gene_name:str):
    with open(gene_data_path_dict[second_gene_name], 'rb') as f:
        exp_dict = pickle.load(f)
    exp_arr = np.array(exp_dict['Sox2'])
    exp_arr[np.isnan(exp_arr)] = 0
    return exp_arr


def get_joint_genes_wn_exp(second_gene_name:str):
    """
    The joint expression will be of sox2 and additional gene.
    This expression will then be used to calculate the covariance in gene expression
    across the ap axis of these two genes that were measured simulatenously in the embryo
    :param second_gene_name:
    :return:
    """
    second_gene_wn_arr = get_wn_exp(second_gene_name)
    joint_wn_exp = np.dstack(
        (get_sox2_wn_exp_in_one_dataset(second_gene_name), second_gene_wn_arr))
    return joint_wn_exp

def get_joint_genes_sc_exp(second_gene_name:str):
    """
    The joint expression SC of sox2 and an additional gene
    :param second_gene_name:
    :return:
    """
    second_gene_sc_arr = get_one_gene_exp_over_AP_axis(second_gene_name)
    joint_sc_exp = np.dstack((get_sox2_sc_exp_in_one_dataset(second_gene_name), second_gene_sc_arr))
    return joint_sc_exp

def plot_all_gasturloid_genes_exp_together():
    cdx2_exp = get_one_gene_exp_over_AP_axis('Cdx2')
    bra_exp = get_one_gene_exp_over_AP_axis('Bra')
    foxc1_exp = get_one_gene_exp_over_AP_axis('Foxc1')
    sox2_exp = get_sox2_exp()
    all_genes_exp_dict = {
        'Cdx2': cdx2_exp,
        'Bra' : bra_exp,
        'Foxc1' : foxc1_exp,
        'Sox2' : sox2_exp
    }
    for gene, expr in all_genes_exp_dict.items():
        mean_expr = np.nanmean(expr, axis=0)  # In case of NaNs
        std_expr = np.nanstd(expr, axis=0)
        x = np.linspace(0,1, expr.shape[1])  # x-axis positions (e.g., 0 to M-1)

        plt.plot(x, (mean_expr)/1000, label=gene,color=GASTRULOID_GENE_COLORS[gene], linewidth=2)
        plt.fill_between(x, (mean_expr - std_expr)/1000, (mean_expr + std_expr)/1000, color=GASTRULOID_GENE_COLORS[gene], alpha=0.5)

    plt.xlabel("x/L")
    plt.ylabel("I (a.u)")
    plt.title("Gastruloid Gene Expression")
    plt.legend()
    plt.xlim(0.1,0.9)
    plt.tight_layout()
    plt.show()


def plot_gastruloids_data(data_dict, dict_name):
    plt.figure(figsize=(12, 6))

    for key in data_dict:
        arr = np.array(data_dict[key])

        arr[np.isnan(arr)] = 0
        # shape should be (49, 200)
        mean_vals = np.mean(arr, axis=0)  # mean across the 49 rows
        std_vals = np.std(arr, axis=0)  # std across the 49 rows

        x = np.arange(arr.shape[1])  # 0 to 199 (200 points)

        plt.plot(x, mean_vals, label=key)
        plt.fill_between(x, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3)

    plt.xlabel('x/L')
    plt.ylabel('Mean ± Std over 49 samples')
    plt.title(f'Mean ± Std Dev over Positions for Each Key- {dict_name}')
    plt.legend()
    #plt.ylim(0,4000)
    #plt.xlim(0.4,0.8)
    plt.tight_layout()
    plt.show()



def create_cov_and_mean_one_gene_sc(gene_name):

    pass

def normalize_gastruloid_gene_expression(training_data, gene_name):
    """
    This function normalizes the given data as done in Petkova et al 2019.
    """
    min_mean_exp, max_mean_exp = min_and_max_mean_gene_expression(training_data, gene_name)
    data_decode = normalize_gene_exp(training_data, np.array(min_mean_exp), np.array(max_mean_exp), gene_name)
    return data_decode

def pos_error_one_gene(gene_name):
    if gene_name != 'Sox2':
        data_path = gene_data_path_dict[gene_name]
        gastru_data = format_gastru_like_droso(load_gastruloid_data(data_path))
    else:
        data_path_bra = gene_data_path_dict['Bra']
        gastru_data_bra =format_gastru_like_droso(load_gastruloid_data(data_path_bra))
    normalized_gastru_data = normalize_gastruloid_gene_expression(gastru_data, [gene_name])
    print(normalized_gastru_data.shape)
    gastru_arr = reshape_gene_data_to_arr(gastru_data,[gene_name])
    print(gastru_arr.shape)
    mean_sc = np.mean(gastru_arr, axis=0).flatten()
    mean_wn = np.vstack(
            (mean_sc[:-2], mean_sc[1:-1],mean_sc[2:])).T
    cov_sc = np.var(gastru_arr, axis = 0).flatten()
    neigh_arr = np.concatenate((gastru_arr[:,:-2,:],gastru_arr[:,1:-1,:], gastru_arr[:,2:,:]),axis=2)
    cov_wn = get_cov(neigh_arr)
    mean_sc = mean_sc[1:-1]
    cov_sc = cov_sc[1:-1]
    sc_pos_err = calculate_position_error_one_gene(cov_sc, mean_sc)
    wn_pos_err = calculate_position_error_gt_pos(cov_wn, mean_wn)
    window_size = 10
    sc_pos_err_smoothened = uniform_filter1d(sc_pos_err[42:], window_size)
    wn_pos_err_smoothened = uniform_filter1d(wn_pos_err[42:], window_size)
    x_pos = np.linspace(0,1,len(sc_pos_err_smoothened))
    plt.scatter(x_pos, sc_pos_err_smoothened/len(sc_pos_err_smoothened), color='blue',label='sc')
    plt.scatter(x_pos,wn_pos_err_smoothened/len(wn_pos_err_smoothened), color='orange',label='wn')
    #plt.scatter(np.linspace(0,1,len(wn_pos_err)), sc_pos_err/len(sc_pos_err), color='blue',label='sc')
    #plt.scatter(np.linspace(0,1,len(wn_pos_err)), wn_pos_err/len(sc_pos_err), color='orange',label='wn')
    plt.title(f'positional error gene:{gene_name} in gastruloids smoothened')
    plt.ylim(0, 1)
    plt.xlim(0, 1)
    plt.legend()
    plt.show()
    return wn_pos_err, sc_pos_err



def get_all_genes_full_wn_mean(sox2_wn_arr, bra_wn_arr, cdx2_wn_arr, foxc1_wn_arr):
    sox2_wn_mean = np.mean(sox2_wn_arr, axis=0)
    bra_wn_mean = np.mean(bra_wn_arr, axis=0)
    cdx2_wn_mean = np.mean(cdx2_wn_arr, axis=0)
    foxc1_wn_mean = np.mean(foxc1_wn_arr, axis=0)
    wn_mean = np.hstack((sox2_wn_mean, bra_wn_mean, cdx2_wn_mean, foxc1_wn_mean))
    return wn_mean

def get_pos_error_all_two_genes_combos(to_plot=True):
    all_two_gene_wn_err = []
    all_two_gene_sc_err = []
    gastru_genes = ['Sox2','Bra','Cdx2','Foxc1']
    for i in range(len(gastru_genes)):
        for j in range(i+1, len(gastru_genes)):
            print(f'pair: {gastru_genes[i]} {gastru_genes[j]}')
            gene1 = gastru_genes[i]
            gene2 = gastru_genes[j]
            if gene1 == 'Sox2':
                wn_pos_err, sc_pos_err = [], []
            #     #wn_pos_err, sc_pos_err = get_pos_error_two_dependent_genes(gene1, gene2)
            else:
                wn_pos_err, sc_pos_err = get_pos_error_two_independent_genes(gene1, gene2)
            all_two_gene_wn_err.append(wn_pos_err)
            all_two_gene_sc_err.append(sc_pos_err)
            if to_plot:
                gene_names = f'{gene1}_{gene2}'
                compare_position_error_all_datasets(wn_pos_err, sc_pos_err,gene_names)
    return all_two_gene_wn_err, all_two_gene_sc_err

def get_pos_error_all_three_gene_combos(to_plot=True):
    all_wn_pos_error_3genes, all_sc_pos_error_3genes = [],[]
    gene_set1 = ['Sox2','Bra','Foxc1']
    gene_set2 = ['Sox2', 'Bra','Cdx2']
    gene_set3 = ['Bra', 'Foxc1','Cdx2']
    gene_set4 = ['Sox2','Foxc1','Cdx2']
    #wn_pos_err_ind, sc_pos_err_ind = get_pos_error_three_independent_genes(gene_set3[0], gene_set3[1], gene_set3[2], to_plot)
    wn_pos_err_dep1 , sc_pos_err_dep1 = get_pos_error_three_genes_with_sox(gene_set1[1], gene_set1[2], to_plot)
    # wn_pos_err_dep2, sc_pos_err_dep2 = get_pos_error_three_genes_with_sox(gene_set2[1], gene_set2[2], to_plot)
    # wn_pos_err_dep3, sc_pos_err_dep3 = get_pos_error_three_genes_with_sox(gene_set4[1], gene_set4[2],to_plot)
    return all_wn_pos_error_3genes, all_sc_pos_error_3genes


def get_pos_error_three_genes_with_sox(gene2:str, gene3:str, to_plot=True):
    #the first gene is sox
    gene1 = 'Sox2'
    #embryosXpositionsX6 (3) per gene, first sox
    gene2_sox2_wn_exp = get_joint_genes_wn_exp(gene2)
    gene3_sox2_wn_exp = get_joint_genes_wn_exp(gene3)


    gene2_sox2_wn_cov = get_cov(gene2_sox2_wn_exp)[:, 3:, :3] #off diagonal covs
    gene3_sox2_wn_cov = get_cov(gene3_sox2_wn_exp)[:, 3:, :3]
    sox2_wn_exp = get_sox2_wn_exp()
    sox2_wn_mean = np.mean(sox2_wn_exp, axis=0)
    gene2_wn_mean = np.mean(gene2_sox2_wn_exp[:, :, 3:], axis=0) #the 3 last columns in axis 2 are of the non Sox gene
    gene3_wn_mean = np.mean(gene3_sox2_wn_exp[:, :, 3:], axis=0)
    wn_mean = np.hstack((sox2_wn_mean, gene2_wn_mean, gene3_wn_mean))


    #diagonal wn covs
    sox2_wn_cov = get_cov(sox2_wn_exp)
    gene2_wn_cov = get_cov(get_wn_exp(gene2))
    gene3_wn_cov = get_cov(get_wn_exp(gene3))

    batch_size = sox2_wn_cov.shape[0]  # 192
    block_size = sox2_wn_cov.shape[1]  # 3
    num_blocks = 3
    final_size = block_size * num_blocks
    full_wn_covs = np.zeros((batch_size, final_size, final_size))
    #the diagonal, the covariance in expression between neighboring positions, same gene
    for i in range(batch_size):
        for j, A in enumerate([sox2_wn_cov, gene2_wn_cov, gene3_wn_cov]):
            start = j * block_size
            end = (j + 1) * block_size
            full_wn_covs[i, start:end, start:end] = A[i]

    for j in range(batch_size):
        for k, B in enumerate([gene2_sox2_wn_cov, gene3_sox2_wn_cov], start=1):
            row_start = k * block_size
            row_end = (k + 1) * block_size
            col_start = 0 * block_size
            col_end = 1 * block_size

            # Lower block: [k,0]
            full_wn_covs[j, row_start:row_end, col_start:col_end] = B[j]

            # Symmetric upper block: [0,k] is B.T
            full_wn_covs[j, col_start:col_end, row_start:row_end] = B[j].T
    print(full_wn_covs.shape)
    full_wn_covs[np.isnan(full_wn_covs)] = 0

    # TODO remaining :  sc covs
    gene2_sc_exp = get_one_gene_exp_over_AP_axis(gene2)
    gene3_sc_exp = get_one_gene_exp_over_AP_axis(gene3)
    sox2_sc_exp = get_sox2_exp()

    sc_mean = np.hstack((np.mean(sox2_sc_exp[:, :, np.newaxis], axis=0), np.mean(gene2_sc_exp[:, :, np.newaxis], axis=0),
                         np.mean(gene3_sc_exp[:, :, np.newaxis],axis=0)))

    # sox2_arr = np.vstack([np.stack(cdx2_sox2['Sox2'].to_list()), np.stack(bra2_sox2['Sox2'].to_list()), np.stack(foxc1_sox2['Sox2'].to_list())])
    sox2_var = np.var(sox2_sc_exp, axis=0)
    gene2_var = np.var(gene2_sc_exp, axis=0)
    gene3_var = np.var(gene3_sc_exp, axis=0)

    sox_gene2_sc_exp = get_joint_genes_sc_exp(gene2)
    sox_gene3_sc_exp = get_joint_genes_sc_exp(gene3)
    sox2_gene2_cov = get_cov(sox_gene2_sc_exp)
    sox2_gene3_cov = get_cov(sox_gene3_sc_exp)

    full_covs_sc = np.zeros((sox2_sc_exp.shape[1], 3, 3))  # will be sox2,bra,cdx2,foxc1 order
    for pos in np.arange(full_covs_sc.shape[0]):
        full_covs_sc[pos, 0, 0] = sox2_var[pos]  # sox2
        full_covs_sc[pos, 0, 1] = sox2_gene2_cov[pos, 0, 1]  # cov sox2 bra
        full_covs_sc[pos, 0, 2] = sox2_gene3_cov[pos, 0, 1]
        full_covs_sc[pos, 1, 0] = sox2_gene2_cov[pos, 0, 1]
        full_covs_sc[pos, 1, 1] = gene2_var[pos]  # 0 cov bra2 , cdx2, and bra2 foxc1 , so 1,2 = 0 , 1,3 = 0
        full_covs_sc[pos, 2, 0] = sox2_gene3_cov[pos, 0, 1]
        full_covs_sc[pos, 2, 2] = gene3_var[pos]
    sc_pos_error = calculate_position_error_gt_pos(full_covs_sc, sc_mean)
    wn_pos_error = calculate_position_error_gt_pos(full_wn_covs, wn_mean)
    if to_plot:
        gene_names = f'{gene1}_{gene2}_{gene3}'
        compare_position_error_all_datasets(wn_pos_error, sc_pos_error, gene_names)
    return wn_pos_error, sc_pos_error




def get_pos_error_three_independent_genes(gene1:str, gene2:str, gene3:str, to_plot=True):
    """
    In the case that the genes are jointly independent.
    In this case, this means that the genes are Bra, Foxc1, and Cdx2 (not necessarily in that order)
    because these genes
    were measured separately """
    gene1_exp = get_one_gene_exp_over_AP_axis(gene1)
    gene2_exp = get_one_gene_exp_over_AP_axis(gene2)
    gene3_exp = get_one_gene_exp_over_AP_axis(gene3)
    gene1_wn_exp = get_wn_exp_from_sc_exp(gene1_exp)
    gene2_wn_exp = get_wn_exp_from_sc_exp(gene2_exp)
    gene3_wn_exp = get_wn_exp_from_sc_exp(gene3_exp)

    num_positions = gene1_exp.shape[1]
    cov_sc = np.zeros((num_positions, 3, 3))
    cov_sc[:, 0, 0] = np.var(gene1_exp, axis=0)
    cov_sc[:, 1, 1] = np.var(gene2_exp, axis=0)
    cov_sc[:, 2, 2] = np.var(gene3_exp, axis=0)

    mean_sc = np.hstack((np.mean(gene1_exp[:, :, np.newaxis], axis=0), np.mean(gene2_exp[:, :, np.newaxis], axis=0), np.mean(gene3_exp[:, :, np.newaxis], axis=0)))
    mean_wn = np.hstack((np.mean(gene1_wn_exp, axis=0), np.mean(gene1_wn_exp, axis=0), np.mean(gene3_wn_exp,axis=0)))
    cov_wn = np.zeros((num_positions-2, 9, 9))
    gene1_wn_cov = get_cov(gene1_wn_exp)
    gene2_wn_cov = get_cov(gene2_wn_exp)
    gene3_wn_cov = get_cov(gene3_wn_exp)
    cov_wn[:, 0:3, 0:3] = gene1_wn_cov
    cov_wn[:, 3:6, 3:6] = gene2_wn_cov
    cov_wn[:, 6:9, 6:9] = gene3_wn_cov
    wn_pos_err = calculate_position_error_gt_pos(cov_wn, mean_wn)
    sc_pos_err = calculate_position_error_gt_pos(cov_sc, mean_sc)
    if to_plot:
        gene_names = f'{gene1}_{gene2}_{gene3}'
        compare_position_error_all_datasets(wn_pos_err, sc_pos_err, gene_names)
    return wn_pos_err, sc_pos_err

def get_pos_error_two_independent_genes(gene1:str, gene2:str):
    #both are not sox2
    gene1_exp = get_one_gene_exp_over_AP_axis(gene1)
    gene2_exp = get_one_gene_exp_over_AP_axis(gene2)
    gene1_wn_exp = get_wn_exp_from_sc_exp(gene1_exp)
    gene2_wn_exp = get_wn_exp_from_sc_exp(gene2_exp)
    num_positions = gene1_exp.shape[1]
    cov_sc = np.zeros((num_positions, 2, 2))
    cov_sc[:, 0, 0] = np.var(gene1_exp, axis=0)
    cov_sc[:, 1, 1] = np.var(gene2_exp, axis=0)
    mean_sc = np.hstack((np.mean(gene1_exp[:,:,np.newaxis],axis=0), np.mean(gene2_exp[:,:,np.newaxis],axis=0)))
    mean_wn = np.hstack((np.mean(gene1_wn_exp, axis=0),np.mean(gene1_wn_exp,axis=0)))
    cov_wn = np.zeros((num_positions-2,6,6))
    gene1_wn_cov = get_cov(gene1_wn_exp)
    gene2_wn_cov = get_cov(gene2_wn_exp)
    cov_wn[:, 0:3, 0:3] = gene1_wn_cov
    cov_wn[:, 3:6, 3:6] = gene2_wn_cov
    wn_pos_err = calculate_position_error_gt_pos(cov_wn, mean_wn)
    sc_pos_err = calculate_position_error_gt_pos(cov_sc, mean_sc)
    return wn_pos_err, sc_pos_err


def get_pos_error_two_dependent_genes(gene1:str, gene2:str):
    # first gene is sox
    exp_arr_sc = get_joint_genes_sc_exp(gene2)
    exp_arr_wn = get_joint_genes_wn_exp(gene2)
    mean_sc = np.mean(exp_arr_sc, axis=0)
    mean_wn = np.mean(exp_arr_wn, axis=0)
    cov_sc = get_cov(exp_arr_sc)

    sox2_gene2_wn_exp = get_joint_genes_wn_exp(gene2)
    sox2_gene_2wn_cov = get_cov(sox2_gene2_wn_exp)
    sc_pos_err = calculate_position_error_gt_pos(cov_sc, mean_sc)
    wn_pos_err = calculate_position_error_gt_pos(sox2_gene_2wn_cov, mean_wn)
    return wn_pos_err, sc_pos_err

def get_all_subsets_pos_error(to_plot=True):
    ##one gene (4)
    #wn_pos_error_1gene, sc_pos_error_1gene = get_pos_error_all_one_gene(to_plot)
    ##two genes (6)
    #wn_pos_error_2genes, sc_pos_error_2genes = get_pos_error_all_two_genes_combos(to_plot)
    ###three genes (4)
    wn_pos_error_3genes, sc_pos_error_3genes = get_pos_error_all_three_gene_combos(to_plot)
    ### four genes (1)
    #wn_pos_4, sc_pos_4 = pos_err_all_four_genes(to_plot)

def get_pos_error_all_one_gene():
    all_wn_pos_error_1gene, all_sc_pos_error_1gene = [],[]
    for gene in gene_data_path_dict.keys():
        wn_pos_error_1gene, sc_pos_error_1gene = pos_error_one_gene(gene)
        all_wn_pos_error_1gene.append(wn_pos_error_1gene)
        all_sc_pos_error_1gene.append(sc_pos_error_1gene)
    return all_wn_pos_error_1gene, all_sc_pos_error_1gene

def create_cov_and_mean_joint_datasets_wn():
    sox2_wn_arr = get_sox2_wn_exp()
    bra_wn_arr = get_wn_exp('Bra')
    cdx2_wn_arr = get_wn_exp('Cdx2')
    foxc1_wn_arr = get_wn_exp('Foxc1')

    #TODO normalize
    #need to calculate the covariance with neighbors for every pair of genes
    #start by getting the expression of the pairs as arrays
    bra2_sox2_wn_exp = get_joint_genes_wn_exp('Bra')
    cdx2_sox2_wn_exp = get_joint_genes_wn_exp('Cdx2')
    foxc1_sox2_wn_exp = get_joint_genes_wn_exp('Foxc1')

    #TODO explain the slicing here
    bra2_sox2_wn_cov = get_cov(bra2_sox2_wn_exp)[:, 3:, :3]
    cdx2_sox2_wn_cov = get_cov(cdx2_sox2_wn_exp)[:, 3:, :3]
    foxc1_sox2_wn_cov = get_cov(foxc1_sox2_wn_exp)[:,3:, :3]

    wn_mean = get_all_genes_full_wn_mean(sox2_wn_arr, bra_wn_arr, cdx2_wn_arr, foxc1_wn_arr)
    #diagonal covs
    sox2_wn_cov = get_cov(sox2_wn_arr)
    bra_wn_cov = get_cov(bra_wn_arr)
    cdx2_wn_cov = get_cov(cdx2_wn_arr)
    foxc1_wn_cov = get_cov(foxc1_wn_arr)

    batch_size = sox2_wn_cov.shape[0]  # 190
    block_size = sox2_wn_cov.shape[1]  # 3
    num_blocks = 4
    final_size = block_size * num_blocks
    full_wn_covs = np.zeros((batch_size, final_size, final_size))
    #the diagonal, the covariance in expression between neighboring positions, same gene
    for i in range(batch_size):
        for j, A in enumerate([sox2_wn_cov, bra_wn_cov, cdx2_wn_cov, foxc1_wn_cov]):
            start = j * block_size
            end = (j + 1) * block_size
            full_wn_covs[i, start:end, start:end] = A[i]

    for j in range(batch_size):
        for k, B in enumerate([bra2_sox2_wn_cov, cdx2_sox2_wn_cov, foxc1_sox2_wn_cov], start=1):
            row_start = k * block_size
            row_end = (k + 1) * block_size
            col_start = 0 * block_size
            col_end = 1 * block_size

            # Lower block: [k,0]
            full_wn_covs[j, row_start:row_end, col_start:col_end] = B[j]

            # Symmetric upper block: [0,k] is B.T
            full_wn_covs[j, col_start:col_end, row_start:row_end] = B[j].T
    print(full_wn_covs.shape)
    full_wn_covs[np.isnan(full_wn_covs)] = 0
    return full_wn_covs, wn_mean
# def create_cov_and_mean_joint_datasets_wn():
#     sox2_wn_arr = get_sox2_wn_exp()
#     bra_wn_arr = get_wn_exp('Bra')
#     cdx2_wn_arr = get_wn_exp('Cdx2')
#     foxc1_wn_arr = get_wn_exp('Foxc1')
#
#     #for covariance between sox2 and each other gene
#     cdx2_sox2 = format_gastru_like_droso(load_gastruloid_data(CDX2_RES_PATH))
#     bra2_sox2 = format_gastru_like_droso(load_gastruloid_data(BRA_RES_PATH))
#     foxc1_sox2 = format_gastru_like_droso(load_gastruloid_data(FOXC1_RES_PATH))
#     #TODO normalize
#
#     #need to calculate the covariance with neighbors for every pair of genes
#     #start by getting the expression of the pairs as arrays
#     bra2_sox2_wn_exp = np.dstack((sliding_window_view(np.stack(bra2_sox2['Sox2'].to_list()),window_shape=3, axis=1),bra_wn_arr))
#     cdx2_sox2_wn_exp = np.dstack((sliding_window_view(np.stack(cdx2_sox2['Sox2'].tolist()), window_shape=3, axis=1), cdx2_wn_arr))
#     foxc1_sox2_wn_exp = np.dstack(
#         (sliding_window_view(np.stack(foxc1_sox2['Sox2'].tolist()), window_shape=3, axis=1), foxc1_wn_arr))
#
#     bra2_sox2_wn_cov = get_cov(bra2_sox2_wn_exp)[:, 3:, :3]
#     cdx2_sox2_wn_cov = get_cov(cdx2_sox2_wn_exp)[:, 3:, :3]
#     foxc1_sox2_wn_cov = get_cov(foxc1_sox2_wn_exp)[:,3:, :3]
#     block_size = 3
#
#     wn_mean = get_all_genes_full_wn_mean(sox2_wn_arr, bra_wn_arr, cdx2_wn_arr, foxc1_wn_arr)
#     #diagonal covs
#     sox2_wn_cov = get_cov(sox2_wn_arr)
#     bra_wn_cov = get_cov(bra_wn_arr)
#     cdx2_wn_cov = get_cov(cdx2_wn_arr)
#     foxc1_wn_cov = get_cov(foxc1_wn_arr)
#
#     batch_size = sox2_wn_cov.shape[0]  # 190
#     block_size = sox2_wn_cov.shape[1]  # 3
#     num_blocks = 4
#     final_size = block_size * num_blocks
#     full_wn_covs = np.zeros((batch_size, final_size, final_size))
#     #the diagonal
#     for i in range(batch_size):
#         for j, A in enumerate([sox2_wn_cov, bra_wn_cov, cdx2_wn_cov, foxc1_wn_cov]):
#             start = j * block_size
#             end = (j + 1) * block_size
#             full_wn_covs[i, start:end, start:end] = A[i]
#
#     for j in range(batch_size):
#         for k, B in enumerate([bra2_sox2_wn_cov, cdx2_sox2_wn_cov, foxc1_sox2_wn_cov], start=1):
#             row_start = k * block_size
#             row_end = (k + 1) * block_size
#             col_start = 0 * block_size
#             col_end = 1 * block_size
#
#             # Lower block: [k,0]
#             full_wn_covs[j, row_start:row_end, col_start:col_end] = B[j]
#
#             # Symmetric upper block: [0,k] is B.T
#             full_wn_covs[j, col_start:col_end, row_start:row_end] = B[j].T
#     # for k in range(batch_size):
#     #     for m, B in enumerate([])
#     print(full_wn_covs.shape)
#     return full_wn_covs, wn_mean




def create_covariance_sc_joint_datasets():
    cdx2_sox2 = format_gastru_like_droso(load_gastruloid_data(CDX2_RES_PATH))
    bra2_sox2 = format_gastru_like_droso(load_gastruloid_data(BRA_RES_PATH))
    foxc1_sox2 = format_gastru_like_droso(load_gastruloid_data(FOXC1_RES_PATH))

    cdx2_arr = get_one_gene_exp_over_AP_axis('Cdx2')
    bra2_arr = get_one_gene_exp_over_AP_axis('Bra')
    foxc1_arr = get_one_gene_exp_over_AP_axis('Foxc1')
    sox2_arr = get_sox2_exp()

    # sox2_arr = np.vstack([np.stack(cdx2_sox2['Sox2'].to_list()), np.stack(bra2_sox2['Sox2'].to_list()), np.stack(foxc1_sox2['Sox2'].to_list())])
    sox2_var = np.var(sox2_arr, axis=0)
    bra2_var = np.var(bra2_arr, axis=0)
    cdx2_var = np.var(cdx2_arr, axis=0)
    foxc1_var = np.var(foxc1_arr, axis=0)

    bra2_sox2_cov = get_cov(reshape_gene_data_to_arr(bra2_sox2,  ['Bra', 'Sox2']))
    cdx2_sox2_cov = get_cov(reshape_gene_data_to_arr(cdx2_sox2,  ['Cdx2', 'Sox2']))
    foxc1_sox2_cov = get_cov(reshape_gene_data_to_arr(foxc1_sox2, ['Foxc1', 'Sox2']))

    full_covs = np.zeros((bra2_arr.shape[1],4,4)) #will be sox2,bra,cdx2,foxc1 order
    for pos in np.arange(full_covs.shape[0]):
        full_covs[pos,0,0] = sox2_var[pos] #sox2
        full_covs[pos,0, 1] = bra2_sox2_cov[pos,0,1] #cov sox2 bra
        full_covs[pos,0, 2] = cdx2_sox2_cov[pos,0,1]
        full_covs[pos,0, 3] = foxc1_sox2_cov[pos,0,1]
        full_covs[pos, 1, 0] = bra2_sox2_cov[pos,0,1]
        full_covs[pos,1,1] = bra2_var[pos] #0 cov bra2 , cdx2, and bra2 foxc1 , so 1,2 = 0 , 1,3 = 0
        full_covs[pos, 2,0] = cdx2_sox2_cov[pos,0,1]
        full_covs[pos,2, 2 ] = cdx2_var[pos]
        full_covs[pos,3,0] = foxc1_sox2_cov[pos,0,1]
        full_covs[pos, 3, 3] = foxc1_var[pos]
    means_sc = np.vstack((np.mean(sox2_arr, axis=0), np.mean(bra2_arr, axis=0), np.mean(cdx2_arr, axis=0),np.mean(foxc1_arr, axis=0))).T
    full_covs[np.isnan(full_covs)] = 0
    return full_covs, means_sc

def calculate_position_error_gt_pos(full_covs, full_means):
    num_pos = full_means.shape[0]
    mean_exp_slopes = np.abs(np.diff(full_means, axis=0))
    mean_exp_slopes = np.vstack([mean_exp_slopes, mean_exp_slopes[-1]])
    position_error = np.zeros(num_pos)
    for pos in range(num_pos):
        position_error[pos] = 1 / (mean_exp_slopes[pos, :] @ np.linalg.inv(
            full_covs[pos, :, :]) @ mean_exp_slopes[pos, :])
    return np.sqrt(position_error)

def calculate_position_error_one_gene(full_covs, means):
    num_pos = means.shape[0]
    mean_exp_slopes = np.diff(means, axis=0)
    mean_exp_slopes = np.append(mean_exp_slopes, mean_exp_slopes[-1])
    position_error = (1/(np.abs(mean_exp_slopes)))*np.sqrt(full_covs)
    return position_error

def plot_positional_information_gastruloids():
    covs_all_gastru_wn, mean_all_gastru_wn = create_cov_and_mean_joint_datasets_wn()
    covs_gastru_sc, means_gastru_sc = create_covariance_sc_joint_datasets()
    wn_pos_error = calculate_position_error_gt_pos(covs_all_gastru_wn,mean_all_gastru_wn)
    sc_pos_error = calculate_position_error_gt_pos(covs_gastru_sc, means_gastru_sc)
    i_sc = np.log2(GASTRULOID_L/((np.sqrt(2*np.pi))*sc_pos_error))[1:-1]
    i_wn = np.log2(GASTRULOID_L/((np.sqrt(2*np.pi))*wn_pos_error))
    i_unique = np.log2(GASTRULOID_L/((np.sqrt(2*np.pi))))*np.ones_like(i_sc)
    x_pos = np.linspace(0,1,len(i_sc))
    plt.plot(x_pos, i_sc, color='blue', label=DECODER_NAMES['sc'])
    plt.plot(x_pos, i_wn, color='orange', label=DECODER_NAMES['wn'])
    plt.plot(x_pos, i_unique, color='black', label='Unique cell specification', linestyle='--')
    plt.legend()
    plt.xlim(0.1,0.9)
    plt.xlabel('position (x/L)')
    plt.ylabel('positional information in bits')
    plt.tight_layout()
    plt.show()

def pos_err_all_four_genes():
    covs_all_gastru_wn, mean_all_gastru_wn = create_cov_and_mean_joint_datasets_wn()
    covs_gastru_sc, means_gastru_sc = create_covariance_sc_joint_datasets()
    wn_pos_error = calculate_position_error_gt_pos(covs_all_gastru_wn, mean_all_gastru_wn)
    sc_pos_error = calculate_position_error_gt_pos(covs_gastru_sc, means_gastru_sc)
    return wn_pos_error, sc_pos_error

def compare_position_error_all_datasets(wn_pos_error, sc_pos_error, genes):
    plt.plot(np.linspace(0, 1, len(wn_pos_error)), sc_pos_error[1:-1], label='sc')
    plt.plot(np.linspace(0, 1, len(wn_pos_error)), wn_pos_error, label='wn')
    plt.legend()
    plt.title(f'pos inf gt \n positions Gastruloids {genes}')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.show()

    print('')

