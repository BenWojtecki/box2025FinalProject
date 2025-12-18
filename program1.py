from readfa import readfq
from xopen import xopen
import argparse
import numpy as np
import time

SIZEMATRIX = 50000
COMP_TABLE = str.maketrans({'A':'T','T':'A','C':'G','G':'C', 'N':'N', 'M': 'K', 'K' : 'M', 'R':'Y', 'Y':'R', 'W': 'W', 'S' : 'S', 'V' : 'B', 'B':'V', 'D':'H', 'H':'D'})
BASES = ['A', 'C', 'G', 'T']

def get_sequences(fileName):
    sequences = []
    with xopen(fileName) as fasta:
        for name, seq, _ in readfq(fasta):
            sequences.append((name, seq, len(seq)))
    return sequences

def allKmers(seq, k):
    if len(seq) < k:
        return set()
    
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmers.add(seq[i:i+k])
    
    return kmers

# Symbole M : Ambiguité entre A ou C 
# Symbole S : Ambiguité entre C ou G 
# Symbole K : Ambiguité entre G ou T
# Symbole W : Ambiguité entre A ou T 
# Symbole R : Ambiguité entre A ou G
# Symbole Y : Ambiguité entre C ou T
# Symbole D : Ambiguité entre A, G ou T
# Symbole H : Ambiguité entre A C T
# Après avoir fait toutes les lettres, on peut faire un mappage classique sans se prendre la tête
# N <-> N, M(A ou C) <-> K(G ou T), R(A ou G) <-> Y(C ou T), W(A ou T) <-> W(A ou T), S(C ou G) <-> S(C ou G), V(A,C,G) <-> B(T,C,G), D(A,G,T) <-> H(A,C,T)
def degCanonical(seq):
    return seq.translate(COMP_TABLE)[::-1]
    
def buildMatrix(seq, rcSeq, kmax):
    dictDb = []
    bases_present = set(seq) | set(rcSeq)
    lx = {base: i for i, base in enumerate(sorted(bases_present))}
    dimX = 4

    for k in range(1, kmax + 1):
        matrix = np.zeros((dimX, dimX), dtype= np.int8) # Astuce de rat

        # Seq
        if len(seq) >= k:
            for i in range(len(seq) - k + 1):
                pref = seq[i:i+k-1]
                suff = seq[i+1:i+k]

                if pref in lx and suff in lx:
                    matrix[lx[pref], lx[suff]] = 1

        # Reverse Complement Seq
        if len(rcSeq) >= k:
            for i in range(len(rcSeq) - k + 1):
                pref = rcSeq[i:i+k-1]
                suff = rcSeq[i+1:i+k]
                if pref in lx and suff in lx:
                    matrix[lx[pref], lx[suff]] = 1

        dictDb.append((matrix, lx.copy()))
        newLx = {}
        cpt = 0

        if len(seq) >= k:
            for i in range(len(seq) - k + 1):
                node = seq[i:i+k]
                if node not in newLx:
                    newLx[node] = cpt
                    cpt += 1

        if len(rcSeq) >= k:
            for i in range(len(rcSeq) - k + 1):
                node = rcSeq[i:i+k]
                if node not in newLx:
                    newLx[node] = cpt
                    cpt += 1

        lx = newLx
        dimX = len(lx)
        
        if dimX > SIZEMATRIX:
            print(f"Stop at k={k} (matrix too big: {dimX})")
            break

    return dictDb

def findMAWS(dictDb, kmax, seq, rc_seq):
    maws = {}
    
    # k = 1: lettres absentes dans seq ET rc_seq
    lettersPresent = set(seq) | set(rc_seq)
    maws[1] = set(BASES) - lettersPresent
    
    # k >= 2
    for k in range(2, kmax + 1):
        if k - 1 > len(dictDb):
            maws[k] = set()
            continue
        
        matrix, lx = dictDb[k-2]  
        kmers_present = allKmers(seq, k) | allKmers(rc_seq, k)
        mawK = set()
        
        for pref in lx.keys():
            for b in BASES:
                candidate = pref + b
                suffix = candidate[1:]
                
                if candidate in kmers_present:
                    continue
                
                if pref not in lx or suffix not in lx:
                    continue

                if matrix[lx[pref], lx[suffix]] == 0:
                    mawK.add(candidate)
        
        maws[k] = mawK
    return maws

def processSequences(seqs, kmax):
    results = {}

    for name, seq, _ in seqs:
        print(f"Traitement de {name} ...")

        rcSeq = degCanonical(seq)
        dictDb = buildMatrix(seq, rcSeq, kmax)
        results[name] = findMAWS(dictDb, kmax, seq, rcSeq)
    
    return results

def writeTSV(results, outFile):
    with open(outFile, "w") as f:
        for seqName in sorted(results):
            for k in sorted(results[seqName]):
                mawSet = results[seqName][k]
                if not mawSet:
                    continue
                f.write(f"{seqName}\t{k}\t{','.join(sorted(mawSet))}\n")

def main():
    parser = argparse.ArgumentParser(description="MAW Detection Program")
    parser.add_argument("fastaFile")
    parser.add_argument("-k", "--kmax", type=int, required=True)
    args = parser.parse_args()
    
    start = time.perf_counter()
    
    seqs = get_sequences(args.fastaFile)
    results = processSequences(seqs, kmax=args.kmax)
    writeTSV(results, "resultsProgram1.tsv")
    
    end = time.perf_counter()
    print(f"Temps total : {end - start:.3f} seconds")

if __name__ == "__main__":
    main()