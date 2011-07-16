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

    # If there's no variation at a residue, then it will mess with the corrcoef
    # function later. Accumulate these in a list to return as well
    novars = []
    for i, counter in enumerate(residue_counters):
        if len(counter) == 1:
            novars.append(i)

    # Efficiently convert a list into a string
    consensus = ''.join(consensus_list)

    return consensus, novars

def strip_positions(data, consensus, novar):
    data = [strip_single_position(seq, novar) for seq in data]
    consensus = strip_single_position(consensus, novar)

    return data, consensus

def strip_single_position(string, novar):
    string = np.array(list(string))
    take = set(range(len(string)))
    take.difference_update(novar)
    take = np.array(sorted(list(take)))
    return "".join(string[take])
    return "".join([char for idx, char in enumerate(string) if idx not in novar])





def generate_binary_matrix(data, consensus):
    """
    Generates a binary array x_i(s), where:
        * Each column corresponds to a particular strain
        * Each row corresponds to a particular site
        * The element is 1 if that strain at that site is indentical to the
           consensus sequence at that site

    """

    x = np.zeros( (len(consensus), len(data)), dtype=bool)
    for s, strain in enumerate(data):
        for i, site in enumerate(strain):
            x[i,s] = (site == consensus[i])

    return x

def find_cutoff(alignment):
    eigs = []

    # We don't want our permutations to mess with the original alignment, which
    # it would do if we don't make a copy of it.  Remember the interlude from
    # day 2?
    alignment = alignment.copy()

    nresidues, nstrains = np.shape(alignment)
    for i in range(1000):
        # Shuffle the residues at each position
        for r in range(nresidues):
            alignment[r, :] = np.random.permutation(alignment[r, :])

        # Calculate the correlation coefficient
        corr = np.corrcoef(alignment)
        
        # Add the eigenvalues to the running list of eigenvalues
        eigs.extend(np.linalg.eigvalsh(corr))

        # Poor-man's Progress bar
        if i%10 == 9: 
            print '.',
        if i%100 == 99: 
            print ''
    print
    return eigs

def clean_matrix(correlation_matrix, lambda_cutoff):
    """ Uses RMT to clean the correlation matrix

    Every eigenvector with an eigenvalue greater than the cutoff is used to
    generate a new correlation matrix
    """
    eigvals, vecs = np.linalg.eigh(correlation_matrix)
    clean = np.zeros_like(correlation_matrix)
    for k, eigval in enumerate(eigvals):
        if eigval > lambda_cutoff and eigval != max(eigvals):
            clean += eigval * np.outer(vecs[:,k], vecs[:,k])
    return clean

def remove_phylogeny(binary_matrix):
    return binary_matrix

def determine_sectors(correlation_matrix, lambda_cutoff):
    """ 
    Determines the sectors of the protein

    Returns a list of lists, where each list contains the residue numbers
    (zero-indexed) of the components of each sector
    """
    eigvals, vecs = np.linalg.eigh(correlation_matrix)

    n_residues = n_vectors = len(eigvals)

    loadings = [Counter() for i in range(n_residues)]

    for r in range(n_residues):
        for k in range(n_vectors):
            loading = np.dot(correlation_matrix[:,r], vecs[:,k])
            if eigvals[k] > lambda_cutoff and loading > .1:
                loadings[r][k] = loading

    print loadings
    best = [(l.most_common()[0][0] if len(l) else None)
            for l in loadings]
    print best


gag_seq_file = '../data/HIV1_ALL_2009_GAG_PRO.fasta'
gag_data_full = read_fasta(gag_seq_file)

# We don't actually need the names of the sequences, and a list is more
# convenient for what we're doing than a dictionary
gag_data = [gag_data_full[name] for name in gag_data_full 
            if 'B' in name.split('.')[0]]

# Remove strains that aren't similar enough to each other
gag_data = filter_strains(gag_data)

# Find the most common residue at each nucleotide location
consensus_sequence, novars = get_consensus(gag_data)

gag_data, consensus_sequence = remove_nonvarying(gag_data, consensus_sequence,
                                                 novars)

x = generate_binary_matrix(gag_data, consensus_sequence)
# x is boolean 2D array, where 
# each row is a residue,
# each column is a strain of HIV,
# and 1 means that it is the same as the consensus sequence


x = remove_phylogeny(x)

corr_matrix = np.corrcoef(x)

#eigs = find_cutoff(x)
#lambda_cutoff = max(eigs)
lambda_cutoff = 3.45

corr_matrix_clean = clean_matrix(corr_matrix, lambda_cutoff)

sectors = determine_sectors(corr_matrix, lambda_cutoff)

sec1 = [1, 88, 2, 94, 3, 97, 4, 99, 5, 100, 6, 108, 8, 118, 9, 120, 11, 122, 12,
        123, 14, 128, 16, 129, 19, 131, 20, 133, 21, 134, 24, 135, 27, 136, 29,
        138, 32, 139, 33, 141, 35, 142, 36, 143, 38, 144, 39, 145, 41, 148, 45,
        149, 48, 150, 50, 151, 51, 152, 52, 153, 57, 154, 60, 155, 63, 156, 73,
        158, 77, 160, 79, 251, 83, 276, 86, 279, 87, 433]

sec2 = [23, 446, 37, 450, 178, 452, 379, 455, 381, 457, 386, 459, 391, 461, 392,
        462, 393, 463, 394, 464, 395, 489, 399, 400, 402, 405, 406, 407, 408,
        412, 413, 414, 416, 417, 419, 420, 423, 430, 431, 432, 434, 435, 437,
        438, 439, 440, 442, 443, 444, 445]

sec3 = [53, 305, 140, 306, 163, 310, 167, 316, 169, 317, 170, 319, 171, 323,
        172, 326, 174, 338, 175, 344, 179, 345, 180, 346, 181, 182, 185, 186,
        187, 189, 191, 198, 199, 212, 221, 225, 229, 233, 240, 243, 245, 249,
        257, 260, 263, 265, 269, 284, 288, 291, 295, 347, 363, 364, 365, 366,
        367]

sec4 = [166, 197, 211, 222, 236, 237, 308, 318, 354, 396]

sec5 = 17, 31, 47, 137, 161, 261, 275, 278, 290, 298, 324, 334, 337, 343]

secQ = [18, 30, 54, 62, 69, 90, 125, 130, 146, 147, 159, 173, 176, 200, 218,
        219, 223, 224, 228, 230, 234, 242, 248, 252, 255, 256, 264, 267, 268,
        273, 280, 281, 286, 312, 341, 357, 362, 374, 375, 376, 401, 403]
