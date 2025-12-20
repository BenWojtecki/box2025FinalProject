"""
Program 1 - Minimal Absent Words (MAWs)
Streaming version (write TSV on the fly, per k)
"""

from readfa import readfq
from xopen import xopen
import argparse
import time

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
        if c not in ENC:
            val = 0
            valid = 0
            continue
        val = ((val << 2) | ENC[c]) & mask
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


def get_sequences(fileName):
    with xopen(fileName) as fasta:
        for name, seq, _ in readfq(fasta):
            yield name, seq


def findMaws_stream(seq, kmax):
    """
    Générateur :
    yield (k, [list of MAWs])
    """
    prevPresent = set()

    for k in range(1, kmax + 1):
        present = set(encodeKmersStream(seq, k))

        # k = 1
        if k == 1:
            maws = []
            for b in BASES:
                if b not in present:
                    maws.append(decode(b, 1))
            if maws:
                yield 1, maws
            prevPresent = present
            continue

        if not prevPresent:
            prevPresent = present
            continue

        maskSuffix = (1 << (2 * (k - 1))) - 1
        mawK = set()

        for prefix in prevPresent:
            for b in BASES:
                x = (prefix << 2) | b

                # Condition 1 : absent
                if x in present:
                    continue

                suffix = x & maskSuffix

                # Condition 2 : minimalité
                if (
                    suffix in prevPresent or
                    rcBin(suffix, k - 1) in prevPresent
                ):
                    xcanon = min(x, rcBin(x, k))
                    mawK.add(xcanon)

        if mawK:
            yield k, sorted(decode(x, k) for x in mawK)

        prevPresent = present


def main():
    parser = argparse.ArgumentParser(description="Minimal Absent Words (Program 1, streaming)")
    parser.add_argument("fastaFile")
    parser.add_argument("-k", type=int, required=True)
    parser.add_argument("-o", default="resultsProgram1.tsv")
    args = parser.parse_args()

    start = time.perf_counter()

    with open(args.o, "w") as out:
        for name, seq in get_sequences(args.fastaFile):
            print(f"Processing {name}")

            for k, maws in findMaws_stream(seq, args.k):
                out.write(f"{name}\t{k}\t{','.join(maws)}\n")

    end = time.perf_counter()
    print(f"Total time: {end - start:.3f} seconds")


if __name__ == "__main__":
    main()
