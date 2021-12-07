#from Bio.Seq import Seq
#from Bio import SeqIO
#from Bio.Alphabet import generic_dna
from Bio import pairwise2


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def aln(seq1, seq2):
    return pairwise2.align.localms(seq1, seq2, 1, 0, -1.1, -.5)[0]
    
def printaln(alignment):
    print (pairwise2.format_alignment(*alignment))

def askn(smth,n):
    try:    out = smth[n]
    except: out = None
    return out

def get_diffs(wt,seq):
    """Outputs the differences between the first sequence
    (usually the wild type) and the second oneM"""
    assert len(seq)==len(wt)
    while wt[-1]=="-":
#         print "wild type ends with '-'!"
        wt = wt[:-1] 
        seq=seq[:-1]
    out = {"InPlace":[],"inserts":[],"deletions":[]}
    i = 0
    for l1,l2 in zip(wt,seq):
        if   l1 == "-" and l2 != "-":
            out["inserts"] += [[i,l2]]
        elif l1 != "-" and l2 == "-":
            i+=1
            out["deletions"] += [i]
        elif l1 != l2:
            i+=1
            out["InPlace"] += [[i,l1, l2]]
        else:
            i+=1
    return out

# def compare_to(wt,mut,trimmed = False, show = False):
#     myaln = aln(wt, mut)
#     try:
#         (s1, s2, sc, begin, end) = myaln[0]
#     except:
#         (s1, s2, sc, begin, end) = myaln
#     if trimmed:
#         if show: printseqs(s1[begin:end],s2[begin:end])
#         return   get_diffs(s1[begin:end],s2[begin:end])
#     else:
#         if show: printseqs(s1,s2)
#         return   get_diffs(s1,s2)

def compare_to(wt,mut,trimmed = False, show = False):
    myaln = aln(wt, mut)
    try:
        (s1, s2, sc, bb, ee) = myaln[0]
    except:
        (s1, s2, sc, bb, ee) = myaln
    if trimmed:
        #find the first base in the wild type aligned string
        begin = min(map(lambda let: str(s1)[  : ].index(let), "atcg")) 
    	#find the last base in the wild type aligned string
        end   = min(map(lambda let: str(s1)[::-1].index(let), "atcg"))
        if show: printseqs(s1[begin:-end],s2[begin:-end])
            
        return   get_diffs(s1[begin:-end],s2[begin:-end])
    else:
        if show:
            printseqs(s1,s2)
        return   get_diffs(s1,s2)

def printseqs(seq1, seq2, onlySecond=False):            
    assert len(seq1)==len(seq2)
    if not onlySecond:
        print (seq1)
    s2 = ""
    for l1,l2 in zip(seq1,seq2):
        if l1==l2: s2 += " "
        else: s2 += l2
    print (s2) #, " ", len(s2.split())


#def printaln(alignment):
#    print (pairwise2.format_alignment(*alignment))

# class Read:
#     '''
#     Read class contains the information about a read
#     Attributes: 
#     seq, score1, begin1, end1, score2, begin2, end2, second_inv
#     '''
#     
#     def __init__(self, row):
#         self.seq    =   Seq(row[0],generic_dna)
#         self.score1 = float(row[1])
#         self.begin1 =   int(row[2])
#         self.end1   =   int(row[3])
#         self.score2 = float(row[4])
#         self.begin2 =   int(row[5])
#         self.end2   =   int(row[6])
#         self.second_inv=int(row[7])