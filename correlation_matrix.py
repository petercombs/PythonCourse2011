import numpy as np
from sequence_tools import read_fasta

gag_seq_file = '../data/HIV1_CON_2004_GAG_PRO.fasta'

gag_data_full = read_fasta(gag_seq_file)

gag_data = [gag_data[name] for name in gag_data_full]

filter_strains(gag_data)

consensus_sequence = get_consensus(gag_data)

x = generate_binary_matrix(gag_data, consensus_sequence)
# x is boolean 2D array, where 
# each column is a residue,
# each row is a strain of HIV,
# and 1 means that it is the same as the consensus sequence

corr_matrix = np.corrcoef(x)

