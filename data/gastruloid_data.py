

from data._preprocessing import *

def load_gastruloid_data(gastuloid_data_path:str):
    """
    :return:
    """
    with open(gastuloid_data_path, 'rb') as f:
        data_dict = pickle.load(f)
    return data_dict

def format_gastru_like_droso(data_dict)->pd.DataFrame:
    num_samples = len(next(iter(data_dict.values())))

    # Initialize a list of dictionaries, one per sample
    rows = []
    for sample_idx in range(num_samples):
        sample_row = {gene: data_dict[gene][sample_idx] for gene in data_dict}
        rows.append(sample_row)

    # Convert list of dicts to DataFrame
    gastru_df = pd.DataFrame(rows)
    return gastru_df