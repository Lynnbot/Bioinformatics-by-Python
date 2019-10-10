# Yuelin Chen  F003kds

import collections, itertools, numpy as np, sys

from Bio.SubsMat import MatrixInfo
from Bio import SeqIO

def symmetrize(subst):
    d = {}
    for ((a,b),v) in subst.items():
        d[a,b] = v
        d[b,a] = v
    return d

blosum62 = symmetrize(MatrixInfo.blosum62)
blosum45 = symmetrize(MatrixInfo.blosum45)
exact = collections.defaultdict(lambda:0, (((a,a),1) for a in 'ACDEFGHIJKLMNPQRSTVWY'))


#%%
# Global NW
def nw(a, b, subst, gap):
    """Perform Needleman-Wunsch global alignment. 
    Args:
      a (string): first sequence
      b (string): second sequence
      subst (dictionary: (char,char)->int): substitution matrix
      gap (negative int): gap penalty
    Returns:
      (int, string, string): alignment score, gapped version of a, gapped version of b
    Example:
      >>> nw('SPADKTNVKAAWGKVGAHAGEYGAEALERMFLS', 'TPEEKSAVTALWGKVNVDEVGGEALGRLLVV', blosum62, -5)
      (62.0, 'SPADKTNVKAAWGKVGAHAGEYGAEALERMFLS', 'TPEEKSAVTALWGKV--NVDEVGGEALGRLLVV')
    """
        
    # TODO your code here
    # Create the scoring matrix
    M = np.zeros((len(a)+1, len(b)+1))
    for i in range(len(a)+1):
        M[i,0] = i*gap             #1st column
    for j in range(len(b)+1):
        M[0,j] = j*gap             # 1st row           
    for i in range(1,len(a)+1):
        for j in range(1,len(b)+1):
            dia = M[i-1,j-1]+subst[a[i-1], b[j-1]]
            ver = M[i, j-1]+gap
            hor = M[i-1,j]+gap
            M[i,j] = max(dia, ver, hor)            
    # Initialize
    score = M[i,j]
    ali_a=''
    ali_b=''
    # Trace back
    while i and j:
        if M[i,j] == M[i-1,j-1]+subst[a[i-1],b[j-1]]:
            ali_a += a[i-1]
            ali_b += b[j-1]
            i -= 1
            j -= 1
        elif M[i,j] == M[i, j-1]+gap: #vertical
            ali_a += '-'
            ali_b += b[j-1]
            j -= 1
        else:
            ali_a += a[i-1]
            ali_b += '-'
            i -= 1
    
    return score, ali_a[::-1], ali_b[::-1]
    
# test
nw('SPADKTNVKAAWGKVGAHAGEYGAEALERMFLS', 'TPEEKSAVTALWGKVNVDEVGGEALGRLLVV', blosum62, -5)


#%%
def sw(a, b, subst, gap):
    """Perform Smith-Waterman local alignment. 
    Args:
      a (string): first sequence
      b (string): second sequence
      subst (dictionary: (char,char)->int): substitution matrix
      gap (negative int): gap penalty
    Returns:
      (int, string, string): alignment score, gapped version of a, gapped version of b
    Example:
      >>> sw('SPADKTNVKAAWGKVGAHAGEYGAEALERMFLS', 'SADQISTVQASFDKVKGDPVGILYAVFKA', blosum62, -5)
      (29.0, 'ADK-TNVKAAWGKV-GAHAG', 'ADQISTVQASFDKVKGDPVG')
    """

    # TODO your code here
    # Create the scoring matrix
    M = np.zeros((len(a)+1, len(b)+1))
    for i in range(len(a)+1):
        M[i,0] = i*gap             #1st column
    for j in range(len(b)+1):
        M[0,j] = j*gap             # 1st row           
    for i in range(1,len(a)+1):
        for j in range(1,len(b)+1):
            dia = M[i-1,j-1]+subst[a[i-1], b[j-1]]
            ver = M[i, j-1]+gap
            hor = M[i-1,j]+gap
            M[i,j] = max(0, dia, ver, hor)       
    score = np.amax(M)
    # Initialize
    ali_a=''; ali_b=''
    # Trace back
    i,j = np.where(M == np.amax(M))
    i = int(i); j=int(j)

    while M[i,j]>0:
        if M[i,j] == M[i-1,j-1]+subst[a[i-1],b[j-1]]:
            ali_a += a[i-1]
            ali_b += b[j-1]
            i -= 1
            j -= 1
        elif M[i,j] == M[i, j-1]+gap: #vertical
            ali_a += '-'
            ali_b += b[j-1]
            j -= 1
        else:
            ali_a += a[i-1]
            ali_b += '-'
            i -= 1
    return score, ali_a[::-1], ali_b[::-1]


sw('SPADKTNVKAAWGKVGAHAGEYGAEALERMFLS', 'SADQISTVQASFDKVKGDPVGILYAVFKA', blosum62, -5)


#%%
def score_tup(seq1, seq2, subst, gap=-4):
    # Assumption: len(seq1) == len(seq2)
    s=0
    for i in range(len(seq1)):
       s += subst[seq1[i], seq2[i]]
    return s

#score_tup('NTT', 'NTP', blosum45,-4)


def make_lookups(query, w, subst, min_word_subst_score):
    """Generate lookup words (and their positions) for a query string.
    Args:
      query (string): sequence to generate queries words for
      w (int): length of lookup words
      subst (dictionary: (char,char)->int): substitution matrix
      min_word_subst_score: only return words with substitution score to original in query at least this value (could be equal)
    Returns:
      dictionary: string -> [int]: for each lookup word (of acceptable score), a list of positions in the query where it appears
    Example:
      >>> make_lookups('VNINV', 3, blosum62, 13)
      {'INI': [0, 2], 'VNI': [0], 'VNV': [0, 2], 'NIN': [1], 'NLN': [1], 'NMN': [1], 'NVN': [1], 'INV': [2]}
    """

    # TODO your code here
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    tup = ''
    tups = []
    for tup in itertools.product(AA, repeat=w):     # 'AAA'-'ZZZ'
        m = ''
        for n in tup:
            m += n
        tups.append(m)
    
    match = {}
    for i in range(len(query)-w+1):
        kmer = query[i:i+w]                  # kmer: snippet of the query seq
        for tup in tups:                     # tups: all possible AA seq
            score = score_tup(tup, kmer, subst=subst)
            if score >= min_word_subst_score:
                if tup in match:
                    match[tup].append(i)
                else:
                    match[tup] = [i]
    return match
    
lookups = make_lookups('VNINV', 3, blosum62, 13)

#%%
def find_hits(lookups, b):
    """Find hits in a database sequence for the lookups.
    Args:
      lookups (dictionary from make_lookups): the words and their locations in the query sequence
      b (string): the database sequence
    Returns:
      [(int,int)]: list of (position in query, position in database) for all matches to the lookups
    Example:
      >>> find_hits({'NIN': [1], 'INV': [2], 'INI': [0,2]}, 'RNINIYNINV')
      [(1, 1), (1, 6), (2, 7), (0, 2), (2, 2)]
    """

    # TODO your code here
    
    w = len(list(lookups.keys())[1])
    match_pos = []                                         
    for hit in lookups.keys():
        for i in range(len(b)-w+1):
            btup = b[i:i+w]
            if hit == btup:
                q_index = lookups[hit]
                for index in q_index:
                    match_pos.append((index, i))# SMALL BUG
    return match_pos

hits = find_hits({'NIN': [1], 'INV': [2], 'INI': [0,2]}, 'RNINIYNINV')
    
#%%
def extend_ungapped(hits, a, b, w, subst, min_single_subst_score=0):
    """Extend initial hits matching substrings of a and b, to the left and to the right, 
    as long as the corresponding substitution scores are good enough.
    Args:
      hits (list from find_hits): the initial matches
      a (string): the query sequence
      b (string): the database sequence
      w (int): length of lookup words, same as used for lookups
      subst (dictionary: (char,char)->int): substitution matrix
      min_single_subst_score: only allow substitutions of at least this value (can be equal)
    Returns:
      [(int,int,int,string,string)]: list of (substitution score, position in a, position in b, substring of a, substring of b) for all matches, possibly extended from input hits
    Example:
      >>> extend_ungapped([(2,3)], 'QKAAATF', 'MIRAAASE', 3, blosum62)
      [(15, 1, 2, 'KAAAT', 'RAAAS')]
    """

    # TODO your code here
    results = []
    for hit in hits:
        a_index = hit[0]
        b_index = hit[1]
        
        # The idea here is to truncate sequence a and b, so that they have the same length 
        # and the hits align at the same time.
        if a_index == b_index:
            a_short = a[:min(len(a), len(b))]
            b_short = b[:min(len(a), len(b))]
        
        if a_index != b_index:
            a_short = a[a_index - min(a_index, b_index):a_index+ w+ min(len(a)-a_index+w, len(b)-b_index+w)+1]
            b_short = b[b_index - min(a_index, b_index):b_index+ w+ min(len(a)-a_index+w, len(b)-b_index+w)+1]
            a_new_index = b_new_index = min(a_index, b_index)
        
        for i in range(min(len(a_short)-a_new_index+w, a_new_index)):
            single_subst_score = score_tup(a_short[a_new_index-i:a_new_index+i], b_short[b_new_index-i:b_new_index+i], subst)
            if single_subst_score > min_single_subst_score: 
                break
            
        a_substr = a_short[a_new_index-i:a_new_index+w+i]
        b_substr = b_short[b_new_index-i:b_new_index+w+i]
        score = score_tup(a_substr, b_substr, subst)
        result = (score, a_index-i, b_index-i, a_substr, b_substr)
    results += result
    return results

extend_ungapped([(2,3)], 'QKAAATF', 'MIRAAASE', 3, blosum62)
        
#%%
def extend_gapped(a, b, i, j, direction, subst, gap, min_gapped_score=0):
    """Extend in direction an alignment of a and b, starting at a[i]/b[j],
    as long as the total score of the extension doesn't drop below the threshold
    Returns the best such extension found over the whole search.

    Args:
      a,b (string): the two sequences
      i,j (int): where the alignment should start, at i in a and j in b
      direction (+1/-1): which way to extend, increasing i and/or j (+1) or decreasing them (-1)
      subst (dictionary: (char,char)->int): substitution matrix
      min_gapped_score: don't extend alignments whose score falls below this
      gap (negative int): gap penalty
    Returns:
      (number, string, string): score and substrings of (a,b) starting at (i,j) and proceeding in direction      
    Example:
      >>> extend_gapped('AAATHISALIGNS','ATHISALIGNS',0,0,1,blosum62,-4)
      (45, 'AAATHISALIGNS', 'A--THISALIGNS')
      >>> extend_gapped('THISALIGNSAAA','THISALIGNSA',12,10,-1,blosum62,-4)
      (45, 'THISALIGNSAAA', 'THISALIGNS--A')
    """

    # TODO your code here
    if direction == -1:
        a = a[i::-1]
        b = b[j::-1]
    else:
        a=a[i:]
        b=b[j:]
    

    cell = {'score': subst(a[i], b[i]), 'pos': (i,j)}
    visit_cell = cell
    
    while visit_cell:
        pass
    



def blast(a, b, w, subst, gap, min_word_subst_score, min_gapped_score):
    """Ungapped BLAST-based alignments.
    Args:
      a (string): first sequence
      b (string): second sequence
      w (int): length of lookup words (typically 3-5)
      subst (dictionary: (char,char)->int): substitution matrix
      gap (negative int): gap penalty
      min_word_subst_score: only return words with substitution score to original in query at least this value (could be equal)
      min_gapped_score: don't extend alignments whose score falls below this
    Returns:
        [(score,a',b'): a list of alignments of substrings of a and b and the alignment score
    Example:
        blast('SPADKTNVKAAWGKVGAHAGEYGAEALERMFLS', 'TPEEKSAVTALWGKVNVDEVGGEALGRLLVV', 3, blosum62, -5, 12, 0)
        => {(64, 'SPADKTNVKAAWGKVGAHAGEYGAEALERMFL', 'TPEEKSAVTALWGKV--NVDEVGGEALGRLLV')}
    """
    lookups = make_lookups(a, w, subst, min_word_subst_score)
    hits = find_hits(lookups, b)
    ungapped_hsps = extend_ungapped(hits, a, b, w, subst)
    gapped_hsps = set() # eliminate duplicates arrived at from different hsps
    for (score0,i0,j0,a0,b0) in ungapped_hsps:
        mi = i0+len(a0)//2; mj = j0+len(b0)//2 # start in middle of hsp
        (score_fw,a_fw,b_fw) = extend_gapped(a,b,mi,mj,1,subst,gap,min_gapped_score)
        (score_bw,a_bw,b_bw) = extend_gapped(a,b,mi-1,mj-1,-1,subst,gap,min_gapped_score)
        gapped_hsps.add((score_bw+score_fw, a_bw+a_fw, b_bw+b_fw))
    return gapped_hsps



'''
# simplistic command-line driver
# call with arguments
#   nw <multi-fasta-file> <idx1> <idx2> <subst> <gap
#     e.g. nw globin_fragments.fa 0 1 blosum62 -4
#   blast <multi-fasta-file> <idx1> <idx2> <subst> <gap> <w> <min_word_subst_score> <min_gapped_score>
#     e.g. blast globin_fragments.fa 0 1 blosum62 -4 3 12 0
if __name__ == '__main__':
    commands = {'nw':nw, 'sw':sw, 'blast':blast}
    command_name = sys.argv[1]
    if command_name not in commands: raise Exception('bad command')
    seqs = [str(r.seq) for r in SeqIO.parse(sys.argv[2], format='fasta')]
    index1 = int(sys.argv[3])
    if index1 < 0 or index1 >= len(seqs): raise Exception('bad index1')
    seq1 = seqs[index1]
    index2 = int(sys.argv[4])
    if index2 < 0 or index2 >= len(seqs): raise Exception('bad index2')
    seq2 = seqs[index2]
    matrices = {'blosum62':blosum62, 'blosum45':blosum45, 'exact':exact}
    if sys.argv[5] not in matrices: raise Exception('bad substitution matrix name')
    subst = matrices[sys.argv[5]]
    if command_name == 'nw' or command_name == 'sw':
        gap = int(sys.argv[6])
        (score, aligned1, aligned2) = commands[command_name](seq1, seq2, subst, gap)
        print('score', score, 'between', index1, 'and', index2)
        print(aligned1)
        print(aligned2)
        sys.exit(0)
    if command_name == 'blast':
        gap = int(sys.argv[6])
        w = int(sys.argv[7])
        min_word_subst_score = int(sys.argv[8])
        min_gapped_score = int(sys.argv[9])
        for (score, aligned1, aligned2) in sorted(blast(seq1, seq2, w, subst, gap, min_word_subst_score, min_gapped_score), reverse=True)[:10]:
            print
            print('score', score, 'between', index1, 'and', index2)
            print(aligned1)
            print(aligned2)
        sys.exit(0)
'''