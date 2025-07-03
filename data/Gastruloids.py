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



def plot_gastruloids_data(data_dict, dict_name):
    plt.figure(figsize=(12, 6))

    for key in data_dict:
        arr = np.array(data_dict[key])  # shape should be (49, 200)
        mean_vals = np.mean(arr, axis=0)  # mean across the 49 rows
        std_vals = np.std(arr, axis=0)  # std across the 49 rows

        x = np.arange(arr.shape[1])  # 0 to 199 (200 points)

        plt.plot(x, mean_vals, label=key)
        plt.fill_between(x, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3)

    plt.xlabel('x/L')
    plt.ylabel('Mean ± Std over 49 samples')
    plt.title(f'Mean ± Std Dev over Positions for Each Key- {dict_name}')
    plt.legend()
    plt.ylim(0,4000)
    plt.tight_layout()
    plt.show()

def create_cov_and_mean_joint_datasets_wn():
    cdx2_sox2 = format_gastru_like_droso(load_gastruloid_data(CDX2_RES_PATH))
    bra2_sox2 = format_gastru_like_droso(load_gastruloid_data(BRA_RES_PATH))
    foxc1_sox2 = format_gastru_like_droso(load_gastruloid_data(FOXC1_RES_PATH))

    cdx2_arr = cdx2_sox2['Cdx2'].to_list()
    cdx2_arr = np.stack(cdx2_arr)  # list of arrays

    bra2_arr = bra2_sox2['Bra'].to_list()
    bra2_arr = np.stack(bra2_arr)

    foxc1_arr = foxc1_sox2['Foxc1'].to_list()
    foxc1_arr = np.stack(foxc1_arr)

    sox2_arr = np.vstack([np.stack(cdx2_sox2['Sox2'].to_list()), np.stack(bra2_sox2['Sox2'].to_list()),
                          np.stack(foxc1_sox2['Sox2'].to_list())])
    sox2_wn_arr = sliding_window_view(sox2_arr, window_shape=3, axis=1)
    bra_wn_arr = sliding_window_view(bra2_arr, window_shape=3, axis=1)
    cdx2_wn_arr = sliding_window_view(cdx2_arr, window_shape=3, axis=1)
    foxc1_arr = sliding_window_view(foxc1_arr, window_shape=3, axis=1)


def create_covariance_sc_joint_datasets():
    cdx2_sox2 = format_gastru_like_droso(load_gastruloid_data(CDX2_RES_PATH))
    bra2_sox2 = format_gastru_like_droso(load_gastruloid_data(BRA_RES_PATH))
    foxc1_sox2 = format_gastru_like_droso(load_gastruloid_data(FOXC1_RES_PATH))

    cdx2_arr = cdx2_sox2['Cdx2'].to_list()
    cdx2_arr = np.stack(cdx2_arr)  # list of arrays

    bra2_arr = bra2_sox2['Bra'].to_list()
    bra2_arr = np.stack(bra2_arr)

    foxc1_arr = foxc1_sox2['Foxc1'].to_list()
    foxc1_arr = np.stack(foxc1_arr)

    sox2_arr = np.vstack([np.stack(cdx2_sox2['Sox2'].to_list()), np.stack(bra2_sox2['Sox2'].to_list()), np.stack(foxc1_sox2['Sox2'].to_list())])

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
    means_wn = np.concatenate(
        (means_sc[:-2, :],means_sc[1:-1, :],means_sc[2:, :]), axis=1)
    return full_covs[42:,:,:], means_sc[42:,:], means_wn[42:,:]

def calculate_position_error_full_exp_profiles_sc(full_covs, means_sc):
    num_pos = means_sc.shape[0]
    mean_exp_slopes = np.diff(means_sc, axis=0)
    mean_exp_slopes = np.vstack([mean_exp_slopes, mean_exp_slopes[-1]])
    position_error = np.zeros(num_pos)
    for pos in range(num_pos):
        position_error[pos] = 1 / (mean_exp_slopes[pos, :] @ np.linalg.inv(
            full_covs[pos, :, :]) @ mean_exp_slopes[pos, :])
    return position_error




def create_covariance_wn_joint_datasets():
    pass