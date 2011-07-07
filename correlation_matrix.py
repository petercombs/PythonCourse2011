import numpy as np
from sequence_tools import read_fasta
from collections import Counter

def filter_strains(data):
    """ Filters out strains that are sufficiently different from the rest
    
    On my first look-through of this data, some of the aligned proteins are
    different lengths.  Since we don't know whether they were clipped at the
    beginning, the end, or what, we'll just take strains that are the most
    common length.

    In principle, you could have other filter conditions in here too!
    """
    # This is a python data type that counts things, and can tell you the most
    # common element
    length_counter = Counter()

    # Count up the lengths of the strains
    for strain in data:
        length_counter[len(strain)] += 1

    # A Counter's most_common method returns a list of tuples, with the first
    # element as the item, and the second element as the count of that item,
    # with the list sorted in descending order by counts
    most_common = length_counter.most_common()[0][0]

    lengths = [len(strain) for strain in data]
    # We need a copy of the data so we don't modify the list while looping over
    # it, which leads to odd behavior.
    good_data = []

    for sequence in data_copy:
        if len(sequence) == most_common:
            good_data.append(sequence)

    # Must not forget the return statement
    return good_data


gag_data_full = read_fasta(gag_seq_file)

gag_data = [gag_data[name] for name in gag_data_full]

# Remove strains that aren't similar enough to each other
gag_data = filter_strains(gag_data)

consensus_sequence = get_consensus(gag_data)

x = generate_binary_matrix(gag_data, consensus_sequence)
# x is boolean 2D array, where 
# each column is a residue,
# each row is a strain of HIV,
# and 1 means that it is the same as the consensus sequence

corr_matrix = np.corrcoef(x)

