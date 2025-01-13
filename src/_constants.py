



"""
Directories:
"""
HOME = "C:/Users/micha/thesis/code/data/position_decoding/drosophila_gregor_optimal_decoding/Data/Gap/"
SAVE_HOME = "C:/Users/Michal/Documents/masters/thesis_res/neighboring_dist_results/"
LAB = "/cs/labs/mornitzan/michal.erez14/thesis/drosophilia_optimalcoding/Data/Data/Gap/"
WHERE = HOME
SAVE_ADDR = SAVE_HOME

TRAINING_DATA_DIR = r'C:\Users\micha\thesis\code\data\position_decoding\drosophila_gregor_optimal_decoding\Data\wt_training_data'
DROSO_RES_DIR = r'C:\Users\micha\thesis\code\data\position_decoding\drosophila_gregor_optimal_decoding\Results\all_results'



MUTANT_PAIR_RULE_FILES_PATHS = r'C:\Users\micha\thesis\code\data\position_decoding\drosophila_gregor_optimal_decoding\Data\Pair_rule\pair_rule_data_raw_dorsal_'
MUTANT_GAP_GENES_FILE_PATHS =r'C:\Users\micha\thesis\code\data\position_decoding\drosophila_gregor_optimal_decoding\Data\Gap\gap_data_raw_dorsal_wt_'

"""
Constants:
"""
LOW_TIME = 40
UPPER_TIME = 44

PAIR_RULE_UPPER_TIME = 55
PAIR_RULE_LOWER_TIME = 45

GAP_GENE_MUTANT_LOWER_TIME = 38
GAP_GENE_MUTANT_UPPER_TIME = 49

EDGE_TRIM = 100

INFORMATION_GAP_THRESHOLD = 1.75

N = 90 #Finding the Last Bits, McGough et. al
N_STD = 4

DECODING_TYPES = ['wn','sc']

NON_GENES_LIST = ['index', 'orient', 'dist', 'age', 'genotype']
GAP_GENES = ['Kni', 'Kr', 'Gt', 'Hb']
ALL_PAIR_RULE_GENES = ['Nuc', 'Eve', 'Run', 'Prd']

WT_PAIR_RULE_NON_GENE_LIST = ['index', 'orient', 'dist', 'age']
WT_PAIR_RULE_GENES = ['Hoechst', 'Eve', 'Run', 'Prd']

PAIR_RULE_GENES =  ['Eve', 'Run', 'Prd']
PAIR_RULE_GENES_NON_GENE_LIST =  NON_GENES_LIST
WT_EMBRYO_FILE_NAMES = [WHERE + "gap_data_raw_dorsal_wt.mat",
                        WHERE + "gap_data_raw_dorsal_wt_bcdE1.mat",
                        WHERE + "gap_data_raw_dorsal_wt_bcd_nos_tsl.mat",
                        WHERE + "gap_data_raw_dorsal_wt_bcd_only_germline_clones.mat",
                        WHERE + "gap_data_raw_dorsal_wt_bcd_osk.mat",
                        WHERE + "gap_data_raw_dorsal_wt_bcd_tsl.mat",
                        WHERE + "gap_data_raw_dorsal_wt_etsl.mat",
                        WHERE + "gap_data_raw_dorsal_wt_osk.mat",
                        WHERE + "gap_data_raw_dorsal_wt_time_series.mat"]

GAP_GENE_COLORS = {
    'Kr': '#8c564b',  # Brown (neutral and earthy)
    'Kni': '#2ca02c', # Green (bright and vibrant)
    'Hb': '#e377c2',  # Pink (far from orange and purple hues)
    'Gt': '#bcbd22'   }# Olive Green (muted yellow-green) }

DECODER_TYPE_COLOR = {'wn':'#ff7f0e','sc':'#0056b3'}
DECODER_NAMES = {'wn':'Neighborhood informed', 'sc': 'Cell independent'}
PAIR_RULE_COLORS = {'Eve':'#008080','Prd':'#800080', 'Run':'#FF8C00' }

MUTANT_PAIR_RULE_DATA = ['bcd_only_germline_clones', 'bcd_osk', 'bcd_tsl', 'bcdE1', 'etsl', 'osk','wt_bcd_nos_tsl']

WT_PAIR_RULE_DATA = 'wt_time_series'

MUTANT_TYPES  = ['etsl', 'bcdE1','osk','bcd_osk','bcd_only_germline_clones', 'bcd_tsl' ]


GAMMA = 1e-4

##plotting constants:
POSITION_X_LABEL = 'position (x/L)'
EXP_Y_LABEL = 'gene expression'
POSITIONS_START = 0.1
POSITIONS_END = 0.9

ONE_GENE_EXAMPLE_GENE = 'Kr'
THREE_GENES_EXAMPLE_GENS = ['Kr', 'Gt', 'Hb']