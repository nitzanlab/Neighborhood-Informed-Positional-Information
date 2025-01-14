from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.ticker import MaxNLocator
import os
class TestResults(ABC):
    def __init__(self,  test_results, test_type, train_mean, train_std, decoding_genes, save=False, savedir='', edge_trim=False):
        self.raw_test_results = test_results #in the shape of num embryo samples X num absolute positions X num positions predicted
        self.test_type = test_type
        self.train_mean = train_mean
        self.train_std = train_std
        self.decoding_genes = decoding_genes
        self.edge_trim = edge_trim
        self.normalized_decoding_map = None
        self.normalize_decoding_map()
        if save:
            self.save(savedir)

    def normalize_decoding_map(self):
        predicted_positions_sums = np.sum(self.raw_test_results, axis=1, keepdims=True)
        self.normalized_decoding_map = self.raw_test_results/predicted_positions_sums


    def measure_weighted_dist_prediction_error(self):
        """
        This function meeasures the weighted distance prediction error:
        Sigma sqrt(p(x)*(x-x_gt)^2 ), the standard deviation
        - x_gt ground truth and x is the predicted position
        this measure is only relevant for the wild type embryos since in this case
        there is a ground truth
        the function returns the mean error per position - actually the variance
        """
        num_embryos = self.raw_test_results.shape[0]
        num_positions = self.raw_test_results.shape[2]
        all_embryos_error = []
        #pos_col = np.arange(1, num_positions + 1).reshape(-1, 1)
        pos_col = np.linspace(0.1,0.9, num_positions).reshape(-1,1)
        pos_columns = np.tile(pos_col, (1, num_positions))
        pos_minus_gt = np.square((pos_columns - pos_columns.T))
        for i in range(num_embryos):
            one_embryo_error = np.sum(self.raw_test_results[i] * pos_minus_gt, axis=0)
            all_embryos_error.append(one_embryo_error)
        return np.sqrt(np.array(all_embryos_error))



    def plot_decoding_map(self,title='', vmax=10, is_wt=True, xlim=None):
    #normalize the results first?

        plt.imshow(np.mean(np.log(self.raw_test_results+(1e-6*np.ones_like(self.raw_test_results))), axis=0), cmap='Greys',
                   origin='lower',vmax=vmax, vmin=-12)

        #plt.title(f'{self.test_type} {title}')
        #plt.xlabel('position along "AP axis"')
        #plt.ylabel('distribution over positions, implied positions')
        colorbar = plt.colorbar()  # Add a colorbar to show the scale
        #plt.savefig(save_dir + f'{self.test_type}_{title}.png')
        colorbar.set_label(''
                           ''
                           'mean log position probability', rotation=270, labelpad=30)

        colorbar.locator = MaxNLocator(integer=True)
        colorbar.update_ticks()

        tick_values = np.linspace(0.1,0.9,9)

        xtick_positions = np.linspace(0, self.raw_test_results.shape[1] - 1, len(tick_values), dtype=int)
        if is_wt:
            plt.plot(range(self.raw_test_results.shape[1]), range(self.raw_test_results.shape[1]), color='red', linestyle='--', linewidth=3,label='ground truth')
            plt.legend()
        plt.xticks(ticks=xtick_positions, labels=[f"{x:.1f}" for x in tick_values], rotation=45)
        plt.yticks(ticks=xtick_positions, labels=[f"{y:.1f}" for y in tick_values])

        if xlim:
            plt.xlim([250, 600])
            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            plt.ylim([250, 600])
        #plt.title(f'{self.test_type} {title}')


        plt.xlabel(r'actual position (x/L)')
        plt.ylabel(r'predicted position ($\hat{x}$/L)')
        #plt.tight_layout()
        plt.tight_layout(rect=[0, 0, 0.95, 1])
        plt.show()

    def calculate_MAP_position(self):
        map_per_position_sc = np.argmax(self.normalized_decoding_map, axis=1)
        return map_per_position_sc

    def calculate_MAP_error(self):
        absolute_locations = np.arange(self.normalized_decoding_map.shape[1])  # the number of rows
        predicted_locations = np.argmax(self.normalized_decoding_map, axis=1)
        prediction_distance = np.abs(predicted_locations - absolute_locations)/self.normalized_decoding_map.shape[1]
        return prediction_distance

    def save(self, savedir):
        filename = os.path.join(savedir, f"{self.test_type}_{'_'.join(self.decoding_genes)}_cell_res.pkl")
        with open(filename, 'wb') as file:
            pickle.dump(self, file)

    @classmethod
    def from_pickle(cls, savedir, test_type, decode_genes, cell_res=False):
        if cell_res:
            filename = os.path.join(savedir, f"{test_type}_{'_'.join(decode_genes)}_cell_res.pkl")
        else:
            filename = os.path.join(savedir, f"{test_type}_{'_'.join(decode_genes)}.pkl")
        with open(filename, 'rb') as file:
            return pickle.load(file)

    #TODO handle that only drosophila data has pair rule genes for reconstruction analysis

    def predicted_pr_map_across_positions(self, wt_mean_pr_exp):
        max_indices = np.argmax(self.normalized_decoding_map, axis=1)
        map_position_predictions = np.zeros_like(self.normalized_decoding_map)
        for i in range(map_position_predictions.shape[0]):
            map_position_predictions[i, max_indices[i], np.arange(self.normalized_decoding_map.shape[2])] = 1
        map_positions_for_mult = np.transpose(map_position_predictions, axes=(
        0, 2, 1))  # need the rows for each embryo to be the distribution for that position
        return map_positions_for_mult @ wt_mean_pr_exp

    def expected_pr_exp(self, wt_mean_pr_exp):
        expected_pr_exp_arr = []
        for i in range(self.normalized_decoding_map.shape[0]):
            per_embryo = self.normalized_decoding_map[i].T@wt_mean_pr_exp
            expected_pr_exp_arr.append(per_embryo)
        return np.array(expected_pr_exp_arr)


def pr_prediction_err(predicted_pr_exp, mean_pr_exp_per_pos):
    per_pos_error_all_pr_genes = np.mean(np.abs(predicted_pr_exp - mean_pr_exp_per_pos), axis=0)
    return per_pos_error_all_pr_genes

