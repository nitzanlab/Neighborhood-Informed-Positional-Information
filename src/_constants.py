
import os

"""
Directories:
"""

HOME_DIR = r'C:\Users\micha\Neighborhood_Informed_PD'
DROSO_RES_DIR = os.path.join(HOME_DIR, 'Droso_Res')
DROSO_DATA_DIR_PATH = os.path.join(HOME_DIR, 'Drosophila_Data')
WHERE = os.path.join(DROSO_DATA_DIR_PATH,'Gap')

MUTANT_PAIR_RULE_DIR_PATH = os.path.join(DROSO_DATA_DIR_PATH, 'Pair_rule\pair_rule_data_raw_dorsal_')
MUTANT_GAP_GENE_DIR_PATH = os.path.join(DROSO_DATA_DIR_PATH, 'Gap\gap_data_raw_dorsal_wt_')

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

WT_DECODE_TYPES = ['wn_wt','sc_wt']

NON_GENES_LIST = ['index', 'orient', 'dist', 'age', 'genotype']
GAP_GENES = ['Kni', 'Kr', 'Gt', 'Hb']
ALL_PAIR_RULE_GENES = ['Nuc', 'Eve', 'Run', 'Prd']

WT_PAIR_RULE_NON_GENE_LIST = ['index', 'orient', 'dist', 'age']
WT_PAIR_RULE_GENES = ['Hoechst', 'Eve', 'Run', 'Prd']

PAIR_RULE_GENES =  ['Eve', 'Run', 'Prd']
PAIR_RULE_GENES_NON_GENE_LIST =  NON_GENES_LIST
WT_EMBRYO_FILE_NAMES = [os.path.join(WHERE,"gap_data_raw_dorsal_wt.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_bcdE1.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_bcd_nos_tsl.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_bcd_only_germline_clones.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_bcd_osk.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_bcd_tsl.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_etsl.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_osk.mat"),
                        os.path.join(WHERE, "gap_data_raw_dorsal_wt_time_series.mat")]

GAP_GENE_COLORS = {
    'Kr': '#8c564b',  # Brown (neutral and earthy)
    'Kni': '#2ca02c', # Green (bright and vibrant)
    'Hb': '#e377c2',  # Pink (far from orange and purple hues)
    'Gt': '#bcbd22'   }# Olive Green (muted yellow-green) }

DECODER_TYPE_COLOR = {'wn':'#ff7f0e','sc':'#0056b3'}
DECODER_NAMES = {'wn':'Neighborhood-informed', 'sc': 'Cell-independent'}
PAIR_RULE_COLORS = {'Eve':'#008080','Prd':'#800080', 'Run':'#FF8C00' }

MUTANT_PAIR_RULE_DATA = ['bcd_only_germline_clones', 'bcd_osk', 'bcd_tsl', 'bcdE1', 'etsl', 'osk','wt_bcd_nos_tsl']

WT_PAIR_RULE_DATA = 'wt_time_series'

MUTANT_TYPES  = ['etsl', 'bcdE1','osk','bcd_osk','bcd_only_germline_clones', 'bcd_tsl' ]

ONE_MUTATION_TYPES = ['etsl', 'bcdE1', 'osk']
DOUBLE_MUTATION_TYPES = ['bcd_osk','bcd_only_germline_clones', 'bcd_tsl']

GAMMA = 1e-4

##plotting constants:
POSITION_X_LABEL = 'position (x/L)'
EXP_Y_LABEL = 'gene expression'
POSITIONS_START = 0.1
POSITIONS_END = 0.9

ONE_GENE_EXAMPLE_GENE = ['Kr']
VMAXS_ONE_GENE = [-4,-4]

THREE_GENES_EXAMPLE_GENES = ['Kr', 'Gt', 'Hb']
VMAXS_THREE_GENES = [0.03, 0.03]
DECODING_TYPES = ['wn', 'sc']