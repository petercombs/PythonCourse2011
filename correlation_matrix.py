import numpy as np
import sys
from Bio import SeqIO, AlignIO
from collections import Counter
from scipy import stats
from argparse import ArgumentParser

def filter_and_generate_binary(data):
    """Removes malformed strains and all sites that have non-varying residues or
    gaps.
    
    The return value x is a 2D boolean array, where:
        *each row corresponds to a residue
        *each column corresponds to a different strain
        *a 1 means that the given strain at that residue is identical to the
            consensus sequence
    """
    # Remove strains that aren't similar enough to each other
    data = filter_strains(data)

    # Find the most common residue at each nucleotide location
    consensus_sequence, novars = get_consensus(data)

    # For scientific and pratical reasons, we should exclude all sites that have
    # are a gap in the consensus strain. These won't be as useful for, say,
    # identifying locations for antibodies, and if we leave them in, we get way
    # too many sectors.
    gaps = [idx for idx, res in enumerate(consensus_sequence) if res is '-']
    novars.extend(gaps)


    data, consensus_sequence = strip_positions(data, consensus_sequence,
                                                     novars)

    x = generate_binary_matrix(data, consensus_sequence)

    return x


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
    length_counter = Counter([len(strain) for strain in data])


    # A Counter's most_common method returns a list of tuples, with the first
    # element as the item, and the second element as the count of that item,
    # with the list sorted in descending order by counts
    most_common = length_counter.most_common()[0][0]

    # Collect only strains that are the right length into the return variable
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
    """Remove positions given in novar from all of the strains as well as the
    consensus strain.
    """
    data = [strip_single_position(seq, novar) for seq in data]
    consensus = strip_single_position(consensus, novar)

    return data, consensus

def strip_single_position(string, novar):
    "Remove positions given in novar from a single string"
    novar = set(novar)  
    # Sets allow constant-time test of membership, as opposed to lists which
    # depend on the length of the list

    return "".join([char for i, char in enumerate(string) 
                    if i not in novar])

def generate_binary_matrix(data, consensus):
    """ Generates a binary array x_i(s), where:
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
    # it might do if we don't make a copy of it.  Remember the interlude from
    # day 2?
    alignment = alignment.copy()

    nresidues, nstrains = np.shape(alignment)
    max_corr = 0
    global allcorrs
    allcorrs = np.empty(100*nresidues**2)

    for i in range(100):
        # Shuffle the residues at each position
        for r in range(nresidues):
            alignment[r, :] = np.random.permutation(alignment[r, :])

        # Calculate the correlation coefficient
        corr = np.corrcoef(alignment, bias=1)
        
        # Add the eigenvalues to the running list of eigenvalues
        eigs.extend(np.linalg.eigvalsh(corr))

        allcorrs[i*nresidues**2:(i+1)*nresidues**2] = \
            abs((corr*~np.identity(nresidues,dtype=bool)).ravel())

        # Poor-man's Progress bar
        if i%10 == 9: 
            print '.',
        if i%100 == 99: 
            print ''
    return eigs

def clean_matrix(correlation_matrix, lambda_cutoff):
    """ Uses RMT to clean the correlation matrix

    Every eigenvector with an eigenvalue greater than the cutoff is used to
    generate a new correlation matrix
    """
    # Calculate the eigenvalues and eigenvectors
    # the h at the end means it is Hermitian (for real values: it's symmetric)
    eigvals, vecs = np.linalg.eigh(correlation_matrix)

    # The "clean" matrix starts out as just zeros
    clean = np.zeros_like(correlation_matrix)
    for k, eigval in enumerate(eigvals):
        if eigval > lambda_cutoff and eigval != max(eigvals):
            # For each eigenvalue larger than the cutoff, compute the outer
            # product, and add it to the clean matrix. This is equation S5
            clean += eigval * np.outer(vecs[:,k], vecs[:,k])
    return clean

def clean_phylogeny(binary_matrix):
    """ Cleans the binary matrix by removing the contribution of phylogeny
    
    This is section S4.

    """

    eigvals, eigvecs = np.linalg.eigh(np.corrcoef(binary_matrix))

    # "u_i^1 are the components of the eigenvector corresponding to the largest
    # eigenvalue"
    u1 = eigvecs[:,eigvals.argmax()]

    num_residues, num_strains = np.shape(binary_matrix)

    # Equation S11
    M = np.array([sum(u1[i] * binary_matrix[i,s] 
                      for i in range(num_residues))
                  for s in range(num_strains)])
    
    # Alpha could be a 1D array, but this is more convenient since we will need
    # to force it into a 2D array eventually
    alpha = np.zeros( (num_residues, num_strains) )
    beta = np.zeros( num_residues )
    epsilon = np.zeros( (num_residues, num_strains) )

    for i in range(num_residues):
        # "The value of the parameters alpha_i and beta_i are estimated through
        # a least square regression..."
        slope, intercept, rval, pval, stderr = stats.linregress(M, x[i,:])
        alpha[i,:] = intercept
        beta[i] = slope
    
    # Equation S10:
    # x_i(s) = alpha_i + beta_i M(s) + epsilon_i(s)
    epsilon = x - alpha - np.outer(beta, M)

    # Equation S12:
    # y_i(s) = alpha_i + epsilon_i(s)
    return alpha + epsilon

def remove_distinct_evo(binary_matrix):
    """ Removes evolutionarily distinct sequences
    
    Calculates a correlation matrix for the strains as they relate to each
    other, and then removes those that are significantly different
    """
    gamma = np.cov(binary_matrix.T)
    eigvals, vecs = np.linalg.eigh(gamma)
    vecs = vecs.T

    # Here there be dragons
    # Using the projections along eigenvector 2 and the cutoff of -.1 was
    # empirically determined. Your mileage may vary
    proj1 = [np.dot(gamma[i], vecs[-1]) for i in range(len(eigvals))]
    proj2 = [np.dot(gamma[i], vecs[-2]) for i in range(len(eigvals))]
    return [pos for pos, proj in enumerate(proj2) if proj > -.1]

def determine_sectors(correlation_matrix, lambda_cutoff):
    """ Determines the sectors of the protein

    See sections S6 and S7 of the supplementals.

    This function returns both the strongest eigenvalue at a given residue and
    the a list of counters with the projection of each residue onto significant
    eigenvectors
    """

    eigvals, vecs = np.linalg.eigh(correlation_matrix)

    n_residues = n_vectors = len(eigvals)

    loadings = [Counter() for i in range(n_residues)]

    # removes the autocorrelations, which should typically be much higher than
    # the inter-residue correlations
    # This works by multiplying by the inverse of the identity matrix
    othercorrs = abs(correlation_matrix 
                     * (~ np.identity(n_residues, dtype=bool)))

    for r in range(n_residues):
        if max(othercorrs[r]) < 0.15: 
            # "We chose to exclude from sectors those residues that did not
            # show any correlation higher than 0.1 (in absolute magnitude) with
            # any other sector residue"
            # I thought their threshold was a little low, since in a given
            # random matrix of about ~500x500, you would expect to see over 1000
            # correlations higher than that by chance, if you believe their P <
            # 5x10^-3 statistic
            continue

        for k in range(n_vectors):
            if eigvals[k] > lambda_cutoff: 
                # The loading is simply the projection of the correlation values for
                # each residue onto a given eigenvector
                loading = np.dot(correlation_matrix[:,r], vecs[:,k])
                loadings[r][k] = abs(loading)

    best = [(l.most_common()[0][0] if (len(l) > 0) else None)
            for l in loadings]

    return best, loadings

def imshow_with_boxes(corr_matrix, list_of_sectors, **kwargs):
    allsecs = np.hstack(list_of_sectors)
    mpl.pcolor(corr_matrix[np.meshgrid(allsecs, allsecs)],**kwargs)

    start = 0
    for sec in list_of_sectors:
        mpl.plot([start, start+len(sec), start+len(sec), start, start],
             [start, start, start+len(sec), start+len(sec), start], 
             'r-')
        start += len(sec)


#############################################################
### Main Program Start
#############################################################

if __name__ == "__main__":
    parser = ArgumentParser(description="Identify sectors of coevolving amino acids" )

    parser.add_argument('infile', nargs='?', 
                        help='The file containing alignments of the protein'
                        'currently supports: FASTA')

    args = parser.parse_args()

    filetype = os.path.splitext(args.infile)[-1]

    if (filetype is 'fasta') or (filetype is 'fa') or (filetype is 'fsa'):
        seq_data_full = [seq for seq in SeqIO.parse(args.infile, 'fasta')]
    else:
        print "The extension '%s' is not currently supported." % filetype



    # We don't actually need the names of the sequences, and a list is more
    # convenient for what we're doing than a dictionary
    seq_data = [seqrec.seq for seqrec in seq_data_full 
                if 'B' in seqrec.name.split('.')[0]]

    print "First Pass Filtering"
    x = filter_and_generate_binary(seq_data)

    distinct_strains = remove_distinct_evo(x)
    seq_data2 = [strain for idx, strain in enumerate(seq_data) if idx in
                 distinct_strains]

    rows, cols = np.shape(x)
    print "Found %d locations in %d strains" % (rows, cols)

    print "Second Pass filtering"

    # The claim in S4 is that the method works best if we remove the contribution
    # from evolutionarily distinct sequences, so we'll have to re-run once we've
    # taken those out.


    x = filter_and_generate_binary(seq_data2)
    x = clean_phylogeny(x)

    rows, cols = np.shape(x)
    print "Found %d locations in %d strains" % (rows, cols)


    print "Building matrix"

    corr_matrix = np.corrcoef(x, bias=1)


    # It actually takes a while for this to run, so I'm going to leave it commented
    # out, and just put in the result that I happen to know it will give us.

    #print "Finding Cutoff"
    #eigs = find_cutoff(x)
    #lambda_cutoff = max(eigs)
    lambda_cutoff = 3.45

    print "Cleaning matrix"
    corr_matrix_clean = clean_matrix(corr_matrix, lambda_cutoff)

    best, loadings = determine_sectors(corr_matrix_clean, lambda_cutoff)

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

    sec5 = [17, 31, 47, 137, 161, 261, 275, 278, 290, 298, 324, 334, 337, 343]

    secQ = [18, 30, 54, 62, 69, 90, 125, 130, 146, 147, 159, 173, 176, 200, 218,
            219, 223, 224, 228, 230, 234, 242, 248, 252, 255, 256, 264, 267, 268,
            273, 280, 281, 286, 312, 341, 357, 362, 374, 375, 376, 401, 403]

    # Takes the hard-coded sectors and combines them into a single list
    secs_them = sec1 + sec2 + sec3 + sec4 + sec5
    # Their coordinates are 1-based, whereas in Python everything is 0-based
    secs_them = np.array(secs_them) - 1

    l,v = np.linalg.eigh(corr_matrix_clean)
    v = v.T

    from matplotlib import pyplot as mpl

    n = 5

    best = np.array(best)



    for i in range(1,n+1):
        for j in range(1,n+1):
            proji = np.dot(corr_matrix_clean.T, v[-i])        
            projj = np.dot(corr_matrix_clean.T, v[-j]) 
            mpl.subplot(n,n,(i-1)*n+j)
            mpl.plot(proji, projj,'o',label='Other', markerfacecolor=(1,1,1,0))
            for sector in sorted(list(set(best))):
                if sector is None:
                    continue
                try:
                    sel = np.where(best == sector)[0]
                    mpl.plot(proji[sel], projj[sel], 'o', label=str(500-sector),)
                except (IndexError, TypeError) as error:
                    print "Skipping sector", sector, "after error", error
                    pass
            if i == 1:
                mpl.title(str(j))
            if j == 1:
                mpl.ylabel(str(i))
            if i == j == n:
                mpl.legend(numpoints=1)

    mpl.show()

    bestl = np.array([(ldg.most_common()[0][1] if len(ldg) else None) for ldg in loadings])

    secs = []
    for sec in sorted(list(set(best)), reverse=True):
        if sec is None:
            continue
        sites = np.where(best == sec)[0]
        ls = bestl[sites]
        secs.append(sites[ls.argsort()])


    secs_me = np.hstack(secs)

    me_only = set(secs_me).difference(secs_them)
    them_only = set(secs_me).symmetric_difference(secs_them) - me_only

    print "I did find   %3d that Dahirel et al didn't" % len(me_only)
    print "Did not find %3d that Dahirel et al did" % len(them_only)

    mpl.figure()

    imshow_with_boxes(corr_matrix_clean, secs)
