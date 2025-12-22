"""
Program 1 - Minimal Absent Words (MAWs)

Definition of MAWs
1) x is absent iff x not in S
2) for every strict substring w of x:
      w in S OR rc(w) in S
"""

from readfa import readfq
from xopen import xopen
import argparse
import time
from array import array

ENC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
DEC = ['A', 'C', 'G', 'T']
BASES = range(4)

def encodeKmersStream(seq, k):
    if len(seq) < k:
        return
    mask = (1 << (2 * k)) - 1
    val = 0
    valid = 0
    
    for c in seq:
        encVal = ENC.get(c, -1)
        if encVal == -1:
            val = 0
            valid = 0
            continue
        val = ((val << 2) | encVal) & mask
        valid += 1
        if valid >= k:
            yield val

def rcBin(x, k):
    rc = 0
    for _ in range(k):
        rc = (rc << 2) | (3 - (x & 3))
        x >>= 2
    return rc

def decode(x, k):
    out = []
    for _ in range(k):
        out.append(DEC[x & 3])
        x >>= 2
    return ''.join(reversed(out))

def getSequences(fileName):
    with xopen(fileName) as fasta:
        for name, seq, _ in readfq(fasta):
            yield name, seq.upper()

def findMawsStream(seq, kmax):
    prevPresent = None
    prevBitArray = None
    
    for k in range(1, kmax + 1):
        kmers = list(encodeKmersStream(seq, k))
        if not kmers:
            break
        
        # Bit array 
        maxKmer = (1 << (2 * k))
        useBitarray = len(kmers) > maxKmer / 100
        
        if useBitarray:
            present = bytearray((maxKmer + 7) // 8)
            for kmer in kmers:
                present[kmer >> 3] |= (1 << (kmer & 7))
            
            def isPresent(x):
                return present[x >> 3] & (1 << (x & 7))
        else:
            present = set(kmers)
            isPresent = present.__contains__

        if k == 1:
            prevPresent = present
            prevBitArray = useBitarray
            continue
        
        if prevPresent is None:
            prevPresent = present
            prevBitArray = useBitarray
            continue
        
        # Vérifie si les k-mers précédents existent
        if (prevBitArray and all(b == 0 for b in prevPresent)) or (not prevBitArray and not prevPresent):
            prevPresent = present
            prevBitArray = useBitarray
            continue
        
        maskSuffix = (1 << (2 * (k - 1))) - 1
        mawK = set()
        
        # Fonction auxiliaire pour vérifier la présence dans l'ensemble précédent
        if prevBitArray:
            def isPrevPresent(x):
                return prevPresent[x >> 3] & (1 << (x & 7))
        else:
            isPrevPresent = prevPresent.__contains__
        
        if prevBitArray:
            # Pour un tableau de bits, on itére sur les bits définis
            maxPrev = (1 << (2 * (k - 1)))
            prevKmers = []
            for i in range(maxPrev):
                if prevPresent[i >> 3] & (1 << (i & 7)):
                    prevKmers.append(i)
        else:
            prevKmers = prevPresent
        
        for prefix in prevKmers:
            for b in BASES:
                x = (prefix << 2) | b
                
                # Condition 1: absent
                if isPresent(x):
                    continue
                
                # Condition 2: minimality
                suffix = x & maskSuffix
                if isPrevPresent(suffix) or isPrevPresent(rcBin(suffix, k - 1)):
                    xcanon = min(x, rcBin(x, k))
                    mawK.add(xcanon)
        
        if mawK:
            yield k, sorted(decode(x, k) for x in mawK)
        
        prevPresent = present
        prevBitArray = useBitarray

def main():
    parser = argparse.ArgumentParser(description="Minimal Absent Words (Program 1)")
    parser.add_argument("fastaFile")
    parser.add_argument("-k", type=int, required=True)
    parser.add_argument("-o", default="resultsProgram1.tsv")
    args = parser.parse_args()
    
    start = time.perf_counter()
    totalMaws = 0
    totalSeqs = 0

    with open(args.o, "w") as out:
        for name, seq in getSequences(args.fastaFile):
            print(f"Processing {name} (length: {len(seq)})")
            totalSeqs +=1
            for k, maws in findMawsStream(seq, args.k):
                out.write(f"{name}\t{k}\t{','.join(maws)}\n")
                totalMaws += (len(maws))
                out.flush()
                
    end = time.perf_counter()
    print(f"\nTotal number of MAWs: {totalMaws}")
    print(f"\nTotal number of Sequences: {totalSeqs}")
    print(f"\nTotal time: {end - start:.3f} seconds")

if __name__ == "__main__":
    main()