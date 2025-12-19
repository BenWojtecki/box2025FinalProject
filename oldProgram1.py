from readfa import readfq
from xopen import xopen
import argparse
import numpy as np
import time
from multiprocessing import Pool, cpu_count

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

def precomputeKmers(seq, rcSeq, kmax):
    kmers_by_k = {}
    for k in range(1, kmax + 1):
        kmers_by_k[k] = allKmers(seq, k) | allKmers(rcSeq, k)
    return kmers_by_k

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
    
def buildMatrix(seq, rcSeq, kmax, kmerByK):
    dictDb = []
    combined = seq + rcSeq
    lx = {base: i for i, base in enumerate(sorted(kmerByK[1]))}
    dimX = len(lx)

    for k in range(1, kmax + 1):
        matrix = np.zeros((dimX, dimX), dtype= np.int8) # Astuce de rat

        # Combined
        if len(combined) >= k:
            for i in range(len(combined) - k + 1):
                pref = combined[i:i+k-1]
                suff = combined[i+1:i+k]

                if pref in lx and suff in lx:
                    matrix[lx[pref], lx[suff]] = 1

        dictDb.append((matrix, lx))
        
        if k < kmax:
            lx = {kmer: i for i, kmer in enumerate(sorted(kmerByK[k+1]))}
            dimX = len(lx)

    return dictDb

def findMAWS(dictDb, kmax, kmersByK):
    maws = {}
    maws[1] = set(BASES) - kmersByK[1]
    
    # k >= 2
    for k in range(2, min(kmax + 1, len(dictDb) + 1)):

        matrix, lx = dictDb[k-2]  
        kmers_present = kmersByK[k]
        mawK = set()
        
        lxSet = set(lx.keys())

        for pref in lxSet:
            i = lx[pref]
            for b in BASES:
                candidate = pref + b
                suffix = candidate[1:]
                
                if candidate in kmers_present:
                    continue
                
                if pref not in lxSet or suffix not in lxSet:
                    continue

                if matrix[i, lx[suffix]] == 0:
                    mawK.add(candidate)
        
        maws[k] = mawK
    return maws

def processOneSequence(args):
    name, seq, kmax = args
    
    print(f"Traitement de {name} ...")
    
    rcSeq = degCanonical(seq)
    kmersByK = precomputeKmers(seq, rcSeq, kmax)
    dictDb = buildMatrix(seq, rcSeq, kmax, kmersByK)
    maws = findMAWS(dictDb, kmax, kmersByK)
    
    return (name, maws)

def processSequences(seqs, kmax, numProcesses=None):
    if numProcesses is None:
        numProcesses = cpu_count()
    
    print(f"Utilisation de {numProcesses} processus pour {len(seqs)} séquence(s)")
    args_list = [(name, seq, kmax) for name, seq, _ in seqs]
    
    if len(seqs) == 1 or numProcesses == 1:
        # Pas besoin de parallélisation
        results = {}
        for args in args_list:
            name, maws = processOneSequence(args)
            results[name] = maws
    else:
        with Pool(processes=numProcesses) as pool:
            result_list = pool.map(processOneSequence, args_list)
        
        results = {name: maws for name, maws in result_list}
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
    parser.add_argument("-p", "--processes", type=int, default=None, 
                        help=f"Number of parallel processes (default: {cpu_count()})")
    args = parser.parse_args()
    
    start = time.perf_counter()
    
    seqs = get_sequences(args.fastaFile)
    results = processSequences(seqs, kmax=args.kmax, numProcesses=args.processes)
    writeTSV(results, "resultsProgram1.tsv")
    
    end = time.perf_counter()
    print(f"Temps total : {end - start:.3f} seconds")

if __name__ == "__main__":
    main()