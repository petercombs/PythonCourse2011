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

    for sequence in data:
        if len(sequence) == most_common:
            good_data.append(sequence)

    # Must not forget the return statement
    return good_data

def get_consensus(strains):
    """ Get the consensus sequence of aligned strains

    This assumes that all strains are the same length, and defines consensus
    sequence as the single most common residue at each position. If two or more
    residues are equally common, it chooses arbitrarily (but not necessarily
    randomly) 
    
    """

    # Set up a list of counters equal to the length of the sequence
    residue_counters = [Counter() for residue in strains[0]]

    for strain in strains:
        # Loop over each strain
        for index, residue in enumerate(strain):
            # Loop over each residue and count how many times that residue has
            # appeared at that position
            residue_counters[index][residue] += 1

    # Use a list comprehension to get the most common residue at each position
    consensus_list = [counter.most_common()[0][0] 
                      for counter in residue_counters]

    for i, counter in enumerate(residue_counters):
        if len(counter) == 1:
            print "No variation at residue: ", i

    # Efficiently convert a list into a string
    consensus = ''.join(consensus_list)

    return consensus

def generate_binary_matrix(data, consensus):
    """
    Generates a binary array x_i(s), where:
        * Each row corresponds to a particular strain
        * Each column corresponds to a particular site
        * The element is 1 if that strain at that site is indentical to the
           consensus sequence at that site

    """

    x = np.zeros( (len(data), len(consensus)), dtype=bool)
    for s, strain in enumerate(data):
        for i, site in enumerate(strain):
            x[s,i] = (site == consensus[i])

    return x

def clean_matrix(correlation_matrix):
    pass

def remove_phylogeny(binary_matrix):
    pass

def determine_sectors(correlation_matrix):
    """ 
    Determines the sectors of the protein

    Returns a list of lists, where each list contains the residue numbers
    (zero-indexed) of the components of each sector
    """
    pass

gag_seq_file = '../data/HIV1_FLT_2009_GAG_PRO.fasta'
gag_data_full = read_fasta(gag_seq_file)

# We don't actually need the names of the sequences, and a list is more
# convenient for what we're doing than a dictionary
gag_data = [gag_data_full[name] for name in gag_data_full]

# Remove strains that aren't similar enough to each other
gag_data = filter_strains(gag_data)

# Find the most common residue at each nucleotide location
consensus_sequence = get_consensus(gag_data)

x = generate_binary_matrix(gag_data, consensus_sequence)
# x is boolean 2D array, where 
# each column is a residue,
# each row is a strain of HIV,
# and 1 means that it is the same as the consensus sequence

x = remove_phylogeny(x)

corr_matrix = np.corrcoef(x)

corr_matrix = clean_matrix(corr_matrix)

sectors = determine_sectors(corr_matrix)

