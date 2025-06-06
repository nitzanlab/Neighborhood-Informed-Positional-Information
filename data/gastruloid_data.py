

from data._preprocessing import *

def load_gastruloid_data(gastuloid_data_path:str):
    """
    This function loads the 38 wild type embryos gap gene expression profiles
    :return:
    """
    with open(gastuloid_data_path, 'rb') as f:
        data_dict = pickle.load(f)
    return data_dict