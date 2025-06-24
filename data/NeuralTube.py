from src._imports import *


from src._constants import *
from data._preprocessing import *
from data.droso_data import *
from data.Data import *
from data.gastruloid_data import *

class NeuralTube(Data):
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
        all_training_data = load_gastruloid_data(self.data_path)
        self.train_data = all_training_data
        self.meta_data = ''
        self.define_data_structures(all_training_data)


    def define_data_structures(self, normalized_data):
        gene_exp_data = normalized_data
        self.genes = {gene: i for i, gene in enumerate(gene_exp_data.keys())}
        #training_arr = reshape_gene_data_to_arr(gene_exp_data, self.genes)
        #TODO handle nans n
        # need to trim and turn to array without nans per gene
        # self.train_data = np.nan_to_num(training_arr, nan=0.0)
        # if self.edge_trim is not None:
        #     self.train_data = self.train_data[:,self.edge_trim:-self.edge_trim,:]


    def train_wn(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        train_data_sbst_genes = self.train_data[:,:,decoding_genes_idx]
        train_wn_data = self.reshape_data_for_wn(train_data_sbst_genes)
        self.learn_mean_wn(train_wn_data,decoding_genes_idx)
        self.learn_covariance_wn(train_wn_data, decoding_genes_idx)
        if self.save_training:
            self.save_dir('wn')

    def learn_mean_sc(self, decoding_genes_idx=np.arange(len(GAP_GENES))):
        gene_means = []
        for gene in self.genes.keys():
            if self.genes[gene] in decoding_genes_idx:
                gene_data = self.train_data[gene]
                gene_data_arr = np.vstack(gene_data)
                gene_means.append(gene_data_arr.mean(axis=0))
        self.means_sc = np.array(gene_means).T

    def learn_covariance_sc(self, decoding_genes_idx=np.arange(len(GAP_GENES))):
        """
        The genes are measured on separate embryos, so we conduct the harshest assumption - that the gene expression
        of the genes is independent
        :param decoding_genes_idx:
        :return:
        """
        covs = np.zeros((self.means_sc.shape[0], 2, 2))
        i=0
        for gene in self.genes.keys():
            if self.genes[gene] in decoding_genes_idx:
                gene_data = self.train_data[gene]
                gene_data_arr = np.vstack(gene_data)
                covs[:,i,i] = np.var(gene_data_arr, axis=0)
                i+=1
        self.std_sc = covs

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

def reshape_gene_data_to_arr(gene_exp_data, genes):
    reshaped_gene_data = []
    for gene in genes:
        reshaped_gene_data.append(np.vstack(gene_exp_data[gene]))
    reshaped_data = np.dstack(np.array(reshaped_gene_data))
    return reshaped_data



def plot_neuraltube_data(data_dict, dict_name):
    plt.figure(figsize=(12, 6))

    for key in data_dict:
        arr = np.array(data_dict[key])
        mean_vals = np.mean(arr, axis=0)
        std_vals = np.std(arr, axis=0)

        x = np.arange(arr.shape[1])

        plt.plot(x, mean_vals, label=key)
        plt.fill_between(x, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3)

    plt.xlabel('x/L')
    plt.ylabel('Mean ± Std over samples')
    plt.title(f'Mean ± Std Dev over Positions for Each Key- {dict_name}')
    plt.legend()
    plt.ylim(0,4000)
    plt.tight_layout()
    plt.show()