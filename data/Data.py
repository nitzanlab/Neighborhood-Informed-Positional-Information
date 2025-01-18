from src._constants import *
from src._imports import *
class Data(ABC):
    def __init__(self, data, data_type, data_subtype=None):
        self.train_data = None
        self.data_type = data_type
        self.data_subtype = data_subtype
        self.means_sc = None
        self.std_sc = None
        self.means_wn = None
        self.covs_wn = None
        self.test_results_wn = None
        self.test_results_sc = None
        self.genes = None

    @abstractmethod
    def preprocess(self, data):
        pass

    def get_decode_genes_idx(self, decoding_genes):
        """
        This function returns the corresponding decoding genes to the gene indices given
        """
        decoding_genes_idx = [self.genes.get(decoding_gene) for decoding_gene in decoding_genes]
        return decoding_genes_idx

    @abstractmethod
    def train_wn(self, decode_genes=None):
        """
        This function trains the neighborhood-informed position decoder
        """
        pass

    @abstractmethod
    def reshape_data_for_wn(self):
        """
        This function reshapes the data for the neighborhood-informed position decoder
        """
        pass

    @abstractmethod
    def learn_mean_wn(self, train_data_wn=None, decoding_genes_idx=None):
        """
        This function learns the mean expression of each gap gene across the AP axis
        for the neighborhood-informed decoder
        """
        pass

    @abstractmethod
    def learn_covariance_wn(self, train_data_wn=None, decoding_genes_idx=None):
        """
        This function learns the covariance in expression of the gap genes across position for
        the neighborhood- informed decoding.
        """
        pass
    @abstractmethod
    def train_sc(self, decoding_genes):
        """
        This function trains the cell-independent decoder
        """
        pass

    @abstractmethod
    def learn_mean_sc(self, decoding_genes_idx=None):
        """
        This function learns the mean expression profiles over the positions
        used for the cell-independent decoder
        """
        pass
    @abstractmethod
    def learn_covariance_sc(self, decoding_genes_idx=None):
        """
        This function learns the covariance of the gap gene expression over positions
        for the cell-independent decoder
        """
        pass

    @abstractmethod
    def test_wn(self, test_data, decoding_genes=None):
        """
        This function predicts the position prediction over positions given the test data
        """
        self.preprocess_data_for_testing_with_neighbors(test_data)
        self.get_position_distribution(test_data)

    def test_sc(self, test_data, decoding_genes=None):
        """
        This function predicts the position decoder over positions given the test data
        """
        self.get_position_distribution(test_data, self.means_sc, self.std_sc)

    def get_log_likelihood(self,x_s, means, covs):
        """
        This function measures the log likelihood distribution of p(g(x)|\hat{x})
        """
        all_pdfs = []
        for mean, cov in zip(means, covs):
            if isinstance(cov, (int, float)):
                cov += GAMMA
            elif np.linalg.matrix_rank(cov) < min(cov.shape):  # check if the matrix is singular
                cov = cov + GAMMA * np.eye(cov.shape[0])
            all_pdfs.append(multivariate_normal.logpdf(x_s, mean=mean, cov=cov))
        all_pdfs = np.array(all_pdfs)
        return all_pdfs

    def get_posterior(self, log_likelihood):
        """
        This function measures the posterior distribution over positions P(\hat{x}|g(x))
        """
        p_x_given_gia = (np.exp(log_likelihood))
        num_positions = p_x_given_gia.shape[0]
        if p_x_given_gia.ndim==1:
            num_positions = len(p_x_given_gia)
        if p_x_given_gia.ndim==2:
            num_positions =  p_x_given_gia.shape[1]
        p_post_unnormalized = p_x_given_gia * (1 / num_positions)
        # after multiplied by the uniform distirbution prior
        col_sums = p_post_unnormalized.sum(axis=0)
        normed_posterior = p_post_unnormalized / col_sums
        # *1/z({g_i}
        return normed_posterior

    def get_position_distribution(self, test_data, means, covs):
        """
        This function measures the posterior distribution for each of the embryos over each positions
        """
        all_maps = []
        num_samples = test_data.shape[0]
        for sample_idx in range(num_samples):
            log_likelihood = self.get_log_likelihood(test_data[sample_idx], means, covs)
            p_x_given_gia = self.get_posterior(log_likelihood)
            all_maps.append(p_x_given_gia)
        return np.array(all_maps)




