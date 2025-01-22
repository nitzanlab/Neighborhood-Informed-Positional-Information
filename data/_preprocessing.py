from src._imports import *
from src._constants import *



def load_data(addr, non_genes_list, genes):
    """
    This function loads the Drosophila data
    :param addr: a string representing the address in the directory
    :return: panda dataframe of the data
    """
    mat = loadmat(addr)  # load mat-file
    mdata = mat['data']  # variable in mat file
    mdtype = mdata.dtype  # dtypes of structures are "unsized objects"
    ndata = {n: mdata[n] for n in mdtype.names}
    columns = [n[0] for n in ndata.items()]  # if v.size == ndata['numIntervals']]
    vals = [ndata[c][0] for c in columns] #each ndata[c] is 1Xnum embryos
    all_vals = [np.concatenate(vals[i]) for i in range(len(columns))]
    all_vals = [np.concatenate(all_vals[i]) if col in non_genes_list else all_vals[i] for i, col in
                enumerate(columns)]
    for i in range(len(non_genes_list), len(all_vals)):
        where_nans = np.isnan(all_vals[i])
        all_vals[i][where_nans] = 0
    df_dict_non_genes = {c: all_vals[i] for i, c in enumerate(non_genes_list)}
    df_dift_genes = {c: (all_vals[i + len(non_genes_list)]).tolist() for i, c in enumerate(genes)}
    df1 = pd.DataFrame(df_dict_non_genes)
    df2 = pd.DataFrame(df_dift_genes)
    df = pd.concat([df1, df2], axis=1)
    return df


def get_normalized_decode_data(file_names, non_genes_list, genes, low_time=LOW_TIME, upper_time=UPPER_TIME):
    """
    This function normalizes the decode data
    """
    loaded_data = load_data(file_names[0], non_genes_list, genes)
    data_time_filtered = time_filter_data(loaded_data, genes, low_time,upper_time)
    min_mean_exp, max_mean_exp = min_and_max_mean_gene_expression(data_time_filtered,genes)
    data_decode = normalize_gene_exp(data_time_filtered, np.array(min_mean_exp), np.array(max_mean_exp), genes)
    return data_decode

def time_filter_data(data_to_filter,genes, lower_thresh, upper_thresh):
    """
    receives data and filters so that only cells in a certain time frame are kept 38-48 min
    :param data:
    :return:
    """
    data = data_to_filter.loc[(data_to_filter['age'] >= lower_thresh) & (data_to_filter['age'] < upper_thresh)]
    return data[genes]


def min_and_max_mean_gene_expression(data, genes):
    """
    THis function calculates the mean gene expression of each location
    and then returns the minimum and maximum expression values per gene
    shoudl normalize data to 0-1
    :param data:
    :param num_genes:
    :param genes:
    :return:
    """
    mean_exps = get_mean(data, len(genes), genes)
    return np.min(mean_exps, axis=0), np.max(mean_exps, axis=0)  # have embyros Xnum_genes, want per column min and max value

def get_mean(data, num_genes, genes):
    """
    This function returns the mean expression across the positions
    """
    gene_means_across_locations = pd.DataFrame(0, index=np.arange(
        len(data[genes[0]].iloc[0])),
                                               columns=genes)  # means[0] is num of genes on number of positions , we want num of positions shape
    for i in range(len(data)):  # loop over the embryos
        for j in range(
                num_genes):  # for each gene, add the expression of that gene by the i'th embryo
            gene_means_across_locations[genes[j]] = \
                gene_means_across_locations[genes[j]] + (data[genes[j]].iloc[i])
    gene_means_across_locations = gene_means_across_locations / len(data)
    return gene_means_across_locations

def normalize_gene_exp(data_to_normalize, min_exp, max_exp, genes):
    """
    This function normalizes the expression of each gene by minimum and maximum
    mean expression over embryos of each gene
    :param data_decode:
    :return:
    """
    data = pd.DataFrame(columns=genes)
    for i, gene in enumerate(genes):
        data[gene] = ((np.array(np.vstack(data_to_normalize[gene])) - min_exp[i]) / (max_exp[i] - min_exp[i])).tolist()
    return data