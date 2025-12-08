from readfa import readfq
from xopen import xopen
import numpy as np



def get_sequences(file_name):
    sequences = []
    with xopen(file_name) as fasta:
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

    lx = ly = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    dimX = dimY = 4

    for k in range(kmax):
        matrix = np.zeros((dimX, dimY))

        if len(seq) < k + 2:
            db.append(matrix)
            dictDb.append((matrix, lx, ly))
            continue 

        # remplissage de la matrice
        for i in range(len(seq) - (k + 1)):
            w = seq[i:i + k + 2]  
            h = w[:-1]   # préfixe
            t = w[-1:]   # suffixe

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

    # k == 1 alors les lettres présentes dans la séquence
    present[1] = set(seq)

    for k in range(2, kmax + 1):
        if k - 1 < len(dictDb):
            _, lx, _ = dictDb[k - 1]
            present[k] = set(lx.keys())
        else:
            present[k] = set()

    for k in range(1, kmax + 1):
        maw_set = set()

        if k == 1:
            alphabet = {"A", "C", "G", "T"}
            for a in alphabet:
                if a not in present[1]:
                    maw_set.add(a)

        else:
            all_prev = present[k - 1]
            all_curr = present[k]

            for pref in all_prev:
                for suff in all_prev:
                    candidate = pref + suff[-1]

                    if candidate in all_curr:
                        continue

                    prefix = candidate[:-1]
                    suffix = candidate[1:]

                    if prefix in all_prev and suffix in all_prev:
                        maw_set.add(candidate)

        maws[k] = maw_set
    return maws

def degCanonical(seq):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    rc = ''.join(comp[b] for b in reversed(seq))
    return min(seq, rc)

def process_sequences(seqList, kmax, useCanonical=False):
    results = {}

    for name, seq, L in seqList:
        print(f"Traitement de {name} ({L} bp)...")
        _, dictDb = buildMatrix(seq, kmax)
        maws = findMAWS(dictDb, kmax, seq)

        if useCanonical:
            for k in maws:
                maws[k] = {degCanonical(x) for x in maws[k]}

        results[name] = maws

    return results

if __name__ == "__main__":

    seq = [("example", "AT", 2)]
    kmax = 3

    db, dictDb = buildMatrix("AT", kmax)

    print("\n### MAWs ###")
    maws = findMAWS(dictDb, kmax, "AT")
    for k in maws:
        print(f"k = {k} : {maws[k]}")

