from readfa import readfq
from xopen import xopen
import argparse
import numpy as np



def get_sequences(fileName):
    sequences = []
    with xopen(fileName) as fasta:
        for name, seq, _ in readfq(fasta):
            sequences.append((name, seq, len(seq)))
    return sequences

def to_dict(matrix, lx, ly):
    newDict = {}
    cpt = 0

    for pref, i in lx.items():
        for suff, j in ly.items():
            if matrix[i, j] == 1:
                newDict[pref + suff] = cpt
                cpt += 1

    return newDict, cpt

def buildMatrix(seq, kmax):
    db = []
    dictDb = []

    lx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    ly = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    dimX = dimY = 4

    for k in range(kmax):
        matrix = np.zeros((dimX, dimY))

        if len(seq) < k + 2:
            db.append(matrix)
            dictDb.append((matrix, lx, ly))
            continue 

        for i in range(len(seq) - (k + 1)):
            w = seq[i:i + k + 2]  
            h = w[:-1]  
            t = w[-1:]

            if h in lx and t in ly:
                matrix[lx[h], ly[t]] = 1

        db.append(matrix)
        dictDb.append((matrix, lx, ly))

        # matrice vide donc dictionnaires non mis à jour
        if np.count_nonzero(matrix) == 0:
            continue

        # sinon mettre à jour 
        ly, dimY = to_dict(matrix, lx, ly)
        lx = ly
        dimX = dimY

    return db, dictDb

def findMAWS(dictDb, kmax, seq):
    maws = {}
    present = {}
    bases = ['A', 'C', 'G', 'T']

    # k == 1 alors les lettres présentes dans la séquence
    kmersSeq = {}

    for k in range(1, kmax + 1):
        kmersSeq[k] = set()
        if len(seq) >= k:
            for i in range(0, len(seq) - k + 1):
                kmer = seq[i:i + k]
                kmersSeq[k].add(kmer)
        else:
            kmersSeq[k] = set()
    
    # Si k = 1
    present[1] = set(seq)

    # Si k > 1
    for k in range(2, kmax + 1):
        if k - 1 < len(dictDb):
            _, lx, _ = dictDb[k - 1]
            present[k] = set(lx.keys())
        else:
            present[k] = set()

    # MAWs
    for k in range(1, kmax + 1):
        maw_set = set()

        if k == 1:
            alphabet = {"A", "C", "G", "T"}
            for a in alphabet:
                if a not in present[1]:
                    maw_set.add(a)

        else:
            all_prev = present[k - 1]
            kmersReal = kmersSeq[k]

            for pref in all_prev:
                for b in bases:
                    candidate = pref + b
                    #Critère de sélection
                    #Ne doit pas être dans la séquence
                    if candidate in kmersReal:
                        continue

                    prefix = candidate[:-1]
                    suffix = candidate[1:]

                    # Facteurs présents dans la séquence
                    if prefix in kmersSeq[k-1] and suffix in kmersSeq[k-1]:
                        maw_set.add(candidate)

        maws[k] = maw_set
    return maws

def degCanonical(seq):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    rc = ''.join(comp[b] for b in reversed(seq))
    return min(seq, rc)

def process_sequences(seqs, kmax):
    results = {}

    for name, seq, L in seqs:
        print(f"Traitement de {name} ...")
        _, dictDb = buildMatrix(seq, kmax)
        maws = findMAWS(dictDb, kmax, seq)

        for k in maws:
            maws[k] = {degCanonical(x) for x in maws[k]}

        results[name] = maws

    return results

def writeTSV(results, outFile):
    with open(outFile, "w") as f:
        for seqName in sorted(results.keys()):
            mawsK = results[seqName]
        
            for k in sorted(mawsK.keys()):
                mawSet = mawsK[k]

                if not mawSet:
                    continue

                sortedMaws = sorted(mawSet)
                mawStr = ",".join(sortedMaws)

                f.write(f"{seqName}\t{k}\t{mawStr}\n")

def main():
    parser = argparse.ArgumentParser(description="Program1")
    parser.add_argument("fastaFile")
    parser.add_argument("-k", "--kmax", type= int, required= True)


    args = parser.parse_args()
    outFile = "results.tsv"
    seqs = get_sequences(args.fastaFile)
    results = process_sequences(seqs, args.kmax)
    writeTSV(results, outFile)



if __name__ == "__main__":
    main()
