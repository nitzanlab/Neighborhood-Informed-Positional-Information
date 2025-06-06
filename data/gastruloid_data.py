

from data._preprocessing import *

def load_gastruloid_data(gastuloid_data_path:str):
    """
    This function loads the 38 wild type embryos gap gene expression profiles
    :return:
    """
    with open(gastuloid_data_path, 'rb') as f:
        data_dict = pickle.load(f)
    return data_dict

def format_gastru_like_droso(gastru_dict)->pd.DataFrame:
    num_samples = len(next(iter(gastru_dict.values())))

    # Initialize a list of dictionaries, one per sample
    rows = []
    for sample_idx in range(num_samples):
        sample_row = {gene: gastru_dict[gene][sample_idx] for gene in gastru_dict}
        rows.append(sample_row)

    # Convert list of dicts to DataFrame
    gastru_df = pd.DataFrame(rows)
    return gastru_df