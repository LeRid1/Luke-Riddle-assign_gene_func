def scoring_function(a, b, substitution_matrix):

    #chosen gap penalty
    if a == '-' or b == '-':
        return -11  
    
    try:
        return substitution_matrix[a, b]
    except KeyError:
        return -1000
    

def fetch_spike_protein(accession_id, Entrez_email):
    from Bio import Entrez, SeqIO
    Entrez.email = Entrez_email

    # Fetch annotated GenBank record (NOT FASTA)
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession_id,
        idtype="acc",
        rettype="gb",
        retmode="text"
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Search for CDS with gene="S" or product containing "spike"
    for feature in record.features:
        if feature.type == "CDS":
            gene = feature.qualifiers.get("gene", [""])
            product = feature.qualifiers.get("product", [""])

            gene_name = gene[0].lower()
            product_name = product[0].lower()

            if gene_name == "s" or "spike" in product_name:
                cds_seq = feature.extract(record.seq)
                protein = cds_seq.translate(to_stop=True)
                return str(protein)

    raise ValueError("Spike protein (S gene) not found in GenBank record.")



def global_alignment(seq1, seq2, scoring_function, substitution_matrix):

    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n, m = len(seq1), len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    score[0][0] = 0
    back  = [[None] * (m + 1) for _ in range(n + 1)]

    gap_penalty = scoring_function('-', '-', substitution_matrix)

    # initialization
    for i in range(1, n + 1):
        score[i][0] = i*-11
        back[i][0] = "up"

    for j in range(1, m + 1):
        score[0][j] = j*-11
        back[0][j] = "left"

    # recurrence
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diag = score[i - 1][j - 1] + scoring_function(seq1[i-1], seq2[j-1], substitution_matrix)
            #print(diag)
            up = score[i - 1][j] + gap_penalty
            #print(up)
            left = score[i][j - 1] + gap_penalty
            #print(left)

            best = max(diag, up, left)
            #print(best)
            score[i][j] = best

            if best == diag:
                back[i][j] = "diag"
            elif best == up:
                back[i][j] = "up"
            else:
                back[i][j] = "left"

    # traceback
    aligned1 = []
    aligned2 = []

    i, j = n, m

    while i > 0 or j > 0:

        if i == 0:
            direction = "left"
        elif j == 0:
            direction = "up"
        else:
            direction = back[i][j]
        if direction == None:
            raise RuntimeError(f"Traceback stuck at i={i}, j={j}, direction={direction}")


        #print (direction)

        if direction == "diag":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == "up":
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1

        elif direction == "left":
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return "".join(aligned1), "".join(aligned2), score[n][m]


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
