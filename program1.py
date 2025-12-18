from readfa import readfq
from xopen import xopen
import argparse
import numpy as np
import time

COMP_TABLE = str.maketrans({'A':'T','T':'A','C':'G','G':'C', 'N':'N', 'M': 'K', 'K' : 'M', 'R':'Y', 'Y':'R', 'W': 'W', 'S' : 'S', 'V' : 'B', 'B':'V', 'D':'H', 'H':'D'})
BASES = ['A', 'C', 'G', 'T']

def get_sequences(fileName):
    sequences = []
    with xopen(fileName) as fasta:
        for name, seq, _ in readfq(fasta):
            sequences.append((name, seq))
    return sequences

def allKmers(seq, k):
    if len(seq) < k:
        return set()
    
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmers.add(seq[i:i+k])
    
    return kmers

def precomputeKmers(seq, rcSeq, kmax):
    kmersByK = {}
    for k in range(1, kmax + 1):
        kmersByK[k] = allKmers(seq, k) | allKmers(rcSeq, k)
    return kmersByK

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
    
def buildTransitions(seq, rcSeq, kmax):
    transitionsByK = {}
    sequences = (seq, rcSeq)

    for k in range(2, kmax + 1):
        trans = {}

        for s in sequences:
            if len(s) < k:
                continue

            for i in range(len(s) - k + 1):
                pref = s[i:i+k-1]
                suff = s[i+1:i+k]

                if pref in trans:
                    trans[pref].add(suff)
                else:
                    trans[pref] = {suff}

        transitionsByK[k] = trans
    return transitionsByK

def findMAWS(transitions_by_k, kmersByK, kmax):
    maws = {}
    
    # k = 1
    maws[1] = set(BASES) - kmersByK[1]

    # k = 2
    for k in range(2, kmax + 1):
        trans = transitions_by_k[k]
        kmersPresent = kmersByK[k]
        kmersPrev = kmersByK[k-1]

        mawK = set()

        for pref in kmersPrev:
            for b in BASES:
                candidate = pref + b
                suff = candidate[1:]

                if candidate in kmersPresent:
                    continue
                if suff not in kmersPrev:
                    continue

                if pref not in trans or suff not in trans[pref]:
                    mawK.add(candidate)

        maws[k] = mawK
    return maws

def processSequences(seqs, kmax):
    results = {}

    for name, seq in seqs:
        print(f"Traitement de {name} ...")

        rcSeq = degCanonical(seq)
        kmersByK = precomputeKmers(seq, rcSeq, kmax)
        transitions = buildTransitions(seq, rcSeq, kmax)
        results[name] = findMAWS(transitions, kmersByK, kmax)

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