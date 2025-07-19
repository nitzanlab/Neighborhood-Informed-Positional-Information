import matplotlib.pyplot as plt

from src._imports import *


from src._constants import *
from data._preprocessing import *
from data.droso_data import *
from data.Data import *
from data.gastruloid_data import *

class NeuralTube(Data):
    def __init__(self, data_path,data_dir=None, data_timepoint=None, data=None, training=False, save_training=False, save_dir=None, load_dir=None, edge_trim=None):
        self.data_dir = data_dir
        self.data_timepoint = data_timepoint
        self.data_path = data_path
        with open(self.data_path, 'rb') as f:
            nt_data = pickle.load(f)
        self.data = nt_data
        self.meta_data = None  # includes orient, dist, age, genotype,..
        self.save_training = save_training
        self.save_dir = save_dir
        self.edge_trim = edge_trim
        super().__init__(data, 'Gastruloid')
        #if training:
        self.preprocess(data, edge_trim)

        # else:
        #     self.means_sc = None
        #     self.covs_wn = None
        #     self.means_wn = None
        #     self.covs_wn = None


    def preprocess(self, data=None, edge_trim=None): #preprcoess training data
        print("Preprocessing Neural Tube data")
        self.define_data_structures()


    def define_data_structures(self):
        self.genes = {gene: i for i, gene in enumerate(self.data.keys())}
        #self.data_arr = reshape_gene_data_to_arr(self.data, self.genes)
        #TODO handle nans n
        # need to trim and turn to array without nans per gene
        # self.train_data = np.nan_to_num(training_arr, nan=0.0)
        # if self.edge_trim is not None:
        #     self.train_data = self.train_data[:,self.edge_trim:-self.edge_trim,:]


    def train_wn(self, decoding_genes):
        decoding_genes_idx = self.get_decode_genes_idx(decoding_genes)
        self.learn_mean_wn(decoding_genes_idx)
        self.learn_covariance_wn(train_wn_data, decoding_genes_idx)
        if self.save_training:
            self.save_dir('wn')

    def learn_mean_sc(self, genes):
        gene_means = []
        for gene in genes:
            if gene in NEURAL_TUBE_SET_A_GENES:
                gene_data = self.data[gene]
                gene_data_arr = np.vstack(gene_data)
                gene_means.append(gene_data_arr.mean(axis=0))
        self.means_sc = np.array(gene_means).T

    def learn_covariance_sc(self, decoding_genes):
        """
        The genes are measured on separate embryos, so we conduct the harshest assumption - that the gene expression
        of the genes is independent
        :param decoding_genes_idx:
        :return:
        """
        covs = np.zeros((self.means_sc.shape[0], 2, 2))
        i=0
        for gene in decoding_genes:
            if gene in NEURAL_TUBE_SET_A_GENES:
                gene_data = self.data[gene]
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
        reshaped_wn_position_means = np.concatenate(
           (self.means_sc[:-2, :], self.means_sc[1:-1, :], self.means_sc[2:, :]), axis=1)
        self.means_wn = reshaped_wn_position_means

    def learn_covariance_wn(self, train_data_wn=None):
        num_positions = self.means_sc.shape[0]-2
        covs_wn = np.zeros((num_positions, len(NEURAL_TUBE_SET_A_GENES)*3, len(NEURAL_TUBE_SET_A_GENES)*3))
        for i,gene in enumerate(NEURAL_TUBE_SET_A_GENES):
            gene_data = self.data[gene]
            gene_data_arr = np.vstack(gene_data)
            gene_data_arr_wn = sliding_window_view(gene_data_arr, window_shape=3, axis=1)
            for pos in range(num_positions):
                covs_wn_pos = np.cov(gene_data_arr_wn[:, pos, :], rowvar=False)
                covs_wn[pos, i, i] = covs_wn_pos[0,0]
                covs_wn[pos, i, i+2] = covs_wn_pos[0, 1]
                covs_wn[pos, i, i + 4] = covs_wn_pos[0, 2]

                covs_wn[pos, i+2, i] = covs_wn_pos[1, 0]
                covs_wn[pos, i+2, i + 2] = covs_wn_pos[1, 1]
                covs_wn[pos, i+2, i + 4] = covs_wn_pos[1, 2]

                covs_wn[pos, i+4, i ] = covs_wn_pos[2, 0]
                covs_wn[pos, i+4, i + 2] = covs_wn_pos[2, 1]
                covs_wn[pos, i+4, i + 4] = covs_wn_pos[2, 2]
        self.covs_wn = covs_wn


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
        positions = np.linspace(0,1, NEURAL_TUBE_BINS)
        for gene in decoding_genes:
            gene_exp = self.data[gene]
            mean_per_pos_one_gene = np.mean(gene_exp, axis=0)
            std_per_pos_one_gene = np.std(gene_exp,axis=0)
            plt.plot(positions,mean_per_pos_one_gene, color=NEURAL_TUBE_COLORS[gene])
            plt.fill_between(positions, mean_per_pos_one_gene - std_per_pos_one_gene,
                             mean_per_pos_one_gene + std_per_pos_one_gene, alpha=0.5, label=gene,
                             color=NEURAL_TUBE_COLORS[gene])
        plt.xlabel(POSITION_X_LABEL)
        plt.ylabel(EXP_Y_LABEL)
        plt.legend()
        plt.tight_layout()
        plt.show()

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
        self.learn_mean_sc(decoding_genes)
        self.learn_covariance_sc(decoding_genes)
        self.learn_mean_wn()
        self.learn_covariance_wn()


    def calculate_position_inf_GT(self, decoding_type):
        if decoding_type == "sc":
            mean_exp = self.means_sc[1:-1,:]
            covs = self.std_sc[1:-1,:,:]
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
            position_error[pos] = 1/(mean_exp_slopes[pos, :] @ np.linalg.inv(
                covs[pos, :, :]) @ mean_exp_slopes[pos, :])
        #TODO add calculation and plot
        return position_error

    def plot_comparison_position_inf_GT(self, genes, title):
        self.calculate_positional_error_per_decoding_map_GT_positions(genes)
        position_error_sc = self.calculate_position_inf_GT('sc')
        position_error_wn = self.calculate_position_inf_GT('wn')
        plt.plot(np.linspace(0,1, len(position_error_sc)) , position_error_sc, label='sc')
        plt.plot(np.linspace(0,1,len(position_error_wn)), position_error_wn, label='wn')
        plt.legend()
        plt.title(f'position information ground truth positions Neural Tube time {title}')
        plt.ylim(0,100)
        plt.show()

def plot_all_time_points(title):
    for time in NEURAL_TUBE_TIMES:
        if title == 'hypo':
            gene_exp_file = f'mutants_h={time}.pkl'
            data_path = NEURAL_TUBE_HYPO_PATH
        else:
            gene_exp_file = f'expressions_h={time}.pkl'
            data_path = NEURAL_TUBE_WT_PATH
        nt_path = os.path.join(data_path, gene_exp_file)
        neuraltube = NeuralTube(data_path=nt_path, training=True, edge_trim=20)
        neuraltube.plot_comparison_position_inf_GT(NEURAL_TUBE_SET_A_GENES, title=f'{time} {title}')

def plot_summarized_neural_tube_over_axis_over_timepoints(genes):
    sc_errors_wt = []
    wn_errors_wt = []
    sc_errors_hypo = []
    wn_errors_hypo = []
    for time in NEURAL_TUBE_TIMES:
        nt_path_hypo = os.path.join(NEURAL_TUBE_HYPO_PATH, f'mutants_h={time}.pkl')
        neuraltube_hypo = NeuralTube(data_path=nt_path_hypo, training=True, edge_trim=20)
        neuraltube_hypo.calculate_positional_error_per_decoding_map_GT_positions(genes)
        position_error_sc_hypo = neuraltube_hypo.calculate_position_inf_GT('sc')
        position_error_wn_hypo = neuraltube_hypo.calculate_position_inf_GT('wn')
        sc_errors_hypo.append(position_error_sc_hypo)
        wn_errors_hypo.append(position_error_wn_hypo)

        nt_path_wt = os.path.join(NEURAL_TUBE_WT_PATH, f'expressions_h={time}.pkl')
        neuraltube_wt = NeuralTube(data_path=nt_path_wt, training=True, edge_trim=20)
        neuraltube_wt.calculate_positional_error_per_decoding_map_GT_positions(genes)
        position_errors_sc_wt = neuraltube_wt.calculate_position_inf_GT('sc')
        position_error_wn_wt = neuraltube_wt.calculate_position_inf_GT('wn')
        sc_errors_wt.append(position_errors_sc_wt)
        wn_errors_wt.append(position_error_wn_wt)
    all_errors = np.concatenate(sc_errors_hypo + wn_errors_hypo + sc_errors_wt + wn_errors_wt)
    ymin, ymax = np.min(all_errors), np.max(all_errors)
    ymin, ymax = max(0, ymin - 0.05 * abs(ymin)), ymax + 0.05 * abs(ymax)
    ymax = 200
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # WT
    bps_sc_wt = axs[0].boxplot(sc_errors_wt, positions=np.array(range(len(NEURAL_TUBE_TIMES))) - 0.15, widths=0.3,
                   patch_artist=True, boxprops=dict(facecolor=DECODER_TYPE_COLOR['sc']))
    bps_wn_wt = axs[0].boxplot(wn_errors_wt, positions=np.array(range(len(NEURAL_TUBE_TIMES))) + 0.15, widths=0.3,
                   patch_artist=True, boxprops=dict(facecolor=DECODER_TYPE_COLOR['wn']))
    axs[0].set_title('WT')
    axs[0].set_xticks(range(len(NEURAL_TUBE_TIMES)))
    axs[0].set_xticklabels(NEURAL_TUBE_TIMES)
    axs[0].set_ylim([ymin, ymax])
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Positional Error')
    axs[0].legend([plt.Rectangle((0, 0), 1, 1, facecolor=DECODER_TYPE_COLOR['sc']), plt.Rectangle((0, 0), 1, 1, facecolor=DECODER_TYPE_COLOR['wn'])],
                  [DECODER_NAMES['sc'], DECODER_NAMES['wn']], loc='upper right')

    # Hypo
    bps_sc_hypo = axs[1].boxplot(sc_errors_hypo, positions=np.array(range(len(NEURAL_TUBE_TIMES))) - 0.15, widths=0.3,
                   patch_artist=True, boxprops=dict(facecolor=DECODER_TYPE_COLOR['sc']))
    bps_wn_hypo = axs[1].boxplot(wn_errors_hypo, positions=np.array(range(len(NEURAL_TUBE_TIMES))) + 0.15, widths=0.3,
                   patch_artist=True, boxprops=dict(facecolor=DECODER_TYPE_COLOR['wn']))
    axs[1].set_title('Hypo')
    axs[1].set_xticks(range(len(NEURAL_TUBE_TIMES)))
    axs[1].set_xticklabels(NEURAL_TUBE_TIMES)
    axs[1].set_ylim([ymin, ymax])
    axs[1].set_xlabel('Time')
    axs[1].legend([plt.Rectangle((0, 0), 1, 1, facecolor=DECODER_TYPE_COLOR['sc']), plt.Rectangle((0, 0), 1, 1, facecolor=DECODER_TYPE_COLOR['wn'])],
                  [DECODER_NAMES['sc'], DECODER_NAMES['wn']], loc='upper right')
    for median_line in bps_sc_wt['medians']:
        median_line.set(color='black', linewidth=2)
    for median_line in bps_wn_wt['medians']:
        median_line.set(color='black', linewidth=2)
    for median_line in bps_sc_hypo['medians']:
        median_line.set(color='black', linewidth=2)
    for median_line in bps_wn_hypo['medians']:
        median_line.set(color='black', linewidth=2)
    plt.tight_layout()
    plt.show()


def plot_positional_information_neural_tube():
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