
from src._constants import *
from data._preprocessing import *
from data._data import *
from data.Data import *

class WildTypeDrosophilaData(Data):
    def __init__(self, data=None, training=False, save_training=False, save_dir=None, load_dir=None, edge_trim=None):
        self.meta_data = None  # includes orient, dist, age, genotype,..
        self.save_training = save_training
        self.save_dir = save_dir
        self.edge_trim = edge_trim
        super().__init__(data, 'WT Drosophila')
        if training:
            self.preprocess(data, edge_trim)

        else:
            self.means_sc = None
            self.covs_wn = None
            self.means_wn = None
            self.covs_wn = None


    def preprocess(self, data=None, edge_trim=None): #preprcoess training data
        print("Preprocessing wild type data")
        ##Here one can decide what is the data used for training
        ## either 254 WT embryos measured between 38-48 minutes into n.c 14
        ### or 38 WT embryos measured between 40-44 minutes into n.c 14.
        ### in both cases, the test data is define as the latter
        ###what is currently uncommented is when using the 38 WT embryos to learn the decode
        ### and predict position estimate

        #all_training_data, training_meta_data = load_all_wt_droso_train_data()
        #self.meta_data = training_meta_data

        #all_training_data,_ = load_all_wt_droso_train_data()
        all_training_data = load_wt_gap_test_data()
        self.meta_data = ''
        self.define_data_structures(all_training_data)


    def define_data_structures(self, normalized_data):
        gene_exp_data = normalized_data[GAP_GENES]
        self.genes = {gene: i for i, gene in enumerate(gene_exp_data.columns)}
        self.train_data = reshape_gene_data_to_arr(gene_exp_data)
        if self.edge_trim is not None:
            self.train_data = self.train_data[:,self.edge_trim:-self.edge_trim,:]


    def train_wn(self, decoding_genes=GAP_GENES):
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

    def test_wn(self, test_data, decoding_genes=GAP_GENES):
        processed_test_data = self.prepare_test_data(test_data, decoding_genes)
        reshaped_test_wn = self.reshape_data_for_wn(processed_test_data)
        decoding_map = self.get_position_distribution(reshaped_test_wn, self.means_wn, self.covs_wn)
        return decoding_map

    def test_sc(self, test_data, decoding_genes=GAP_GENES):
        processed_test_data = self.prepare_test_data(test_data, decoding_genes)
        decoding_map = self.get_position_distribution(processed_test_data, self.means_sc, self.std_sc)
        return decoding_map
    def prepare_test_data(self, test_data, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        processed_test_data = reshape_gene_data_to_arr(test_data)[:, :, decoding_genes_idx]
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
        self.means_wn = np.mean(train_data_wn, axis=0)

    def learn_covariance_wn(self, train_data_wn=None):
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

def get_cov(training_data):
    num_positions = training_data.shape[1]
    num_features = training_data.shape[2]
    covs = np.zeros((num_positions, num_features, num_features))
    for pos in range(num_positions):
        covs[pos] = np.cov(training_data[:, pos, :], rowvar=False)
    return covs

def reshape_gene_data_to_arr(gene_exp_data):
    reshaped_gene_data = []
    for gene in GAP_GENES:
        reshaped_gene_data.append(np.vstack(gene_exp_data[gene]))
    reshaped_data = np.dstack(np.array(reshaped_gene_data))
    return reshaped_data

