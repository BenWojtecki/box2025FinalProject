from readfa import readfq
from xopen import xopen
import argparse
import numpy as np
import time

SIZEMATRIX = 50000
COMP_TABLE = str.maketrans({'A':'T','T':'A','C':'G','G':'C', 'N':'N', 'M': 'K', 'K' : 'M', 'R':'Y', 'Y':'R', 'W': 'W', 'S' : 'S', 'V' : 'B', 'B':'V', 'D':'H', 'H':'D'})


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
        matrix = np.zeros((dimX, dimY), dtype= np.int8) # Astuce de rat comme on aime 

        if len(seq) < k + 2:
            db.append(matrix)
            dictDb.append((matrix, lx, ly))
            continue 

        seqLength = len(seq)

        for i in range(seqLength - (k + 1)):
            w = seq[i:i + k + 2]  
            h = w[:-1]  
            t = w[-1:]

            if h in lx and t in ly:
                matrix[lx[h], ly[t]] = 1

       
        # matrice vide donc dictionnaires non mis à jour
        if not matrix.any() == 0: # Plus rapide que np.nonzero normalement ^^
            continue

        db.append(matrix)
        dictDb.append((matrix, lx, ly))

        if dimX > SIZEMATRIX:
            print(f"Stop at k = {k+1} (matrix to big): {dimX} inputs")
            break

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
    seqLength = len(seq)
    kmersSeq = {}

    for k in range(1, kmax + 1):
        kmersSeq[k] = set()

    for i in range(seqLength):
        word = ""
        maxK = min(kmax, seqLength - i)
        for k in range(1, maxK + 1):
            word += seq[i + k - 1]
            kmersSeq[k].add(word)


    present[1] = set(seq)

    for k in range(2, kmax + 1):
        if k - 1 < len(dictDb):
            _, lx, _ = dictDb[k - 1]
            present[k] = set(lx.keys())
        else:
            present[k] = set()

    # MAWs
    for k in range(1, kmax + 1):
        if k > 1 and not present[k - 1]:
            maws[k] = set()
            continue

        maw_set = set()

        if k == 1:
            maw_set = set(bases) - present[1]

        else:
            prevKmers = present[k - 1]
            kmersReal = kmersSeq[k]

            for pref in prevKmers:
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

# Gros poblème ici car dans le fichier all_ebi, nous avons de N et des M. D'après la doc du format .fa, le N représente une base non déterminé donc pas de soucis pour ça
# Mais le symbole M est une ambiguité entre A ou C, peut être ajouter une lettre X qui représente les complémentaires de A ou C  donc G ou T mais pas sur de cette idée du tout
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
# Maintenant ça devrait tourner
def degCanonical(seq):
    rc = seq.translate(COMP_TABLE)[::-1]
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
    outFile = "resultsProgram1.tsv"
    start = time.perf_counter()
    seqs = get_sequences(args.fastaFile)
    results = process_sequences(seqs, args.kmax)
    writeTSV(results, outFile)
    end = time.perf_counter()
    print(f"Temps total : {end - start:.3f} seconds")


if __name__ == "__main__":
    main()
