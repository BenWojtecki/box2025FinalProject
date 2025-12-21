# Software report

## Program1

### Overview

The program computes all canonical minimal absent words up to a maximum length `kmax` for each DNA sequence in a FASTA file. The implementation prioritizes clarity, correctness, and reasonable performance for moderate values of `kmax` (approximately 15 kmers). The program is organized as a single Python file structured into four main components:

1. Encoding utilities (binary encoding of DNA sequences)
2. Sequence input (FASTA parsing)
3. MAWs computation
4. Command line interface and output

Each component is implemented as a set of independent functions to improve readability and modularity.


### Data Representation

#### DNA Encoding

DNA bases are encoded using a 2 bit representation: A -> 00, C -> 01, G -> 10, T -> 11. This encoding allows efficient manipulation of kmers as integers and enables fast computation of reverse complements and for findign MAWs. The mappings are implemented via `ENC`: dictionary (base → binary) and `DEC`: list (binary → base). Depending on kmer density, the program dynamically chooses between a python set or bit array (memory efficient for dense kmer spaces). The choice is based on the ratio between the number of observed kmers and the theoretical maximum (4^k). This strategy reduces memory overhead.

## 4. Main functions

The function `encodeKmersStream(seq, k)` generates all kmers from a DNA sequence and encodes them as integers using a rolling approach. Symbols that are not in the DNA alphabet (A, C, G, T) are ignored. Bit masking is used to keep the kmer length fixed, and the function operates in streaming mode to limit memory usage. Reverse complements are computed on the binary representation of kmers using the function `rcBin(x, k)`, which runs in linear time with respect to k. This approach avoids string-based operations. The function `findMawsStream` implements the main algorithm of program. For each value of k from 1 to kmax, it enumerates all kmers present in the sequence, identifies MAW. Each valid MAW is converted to its canonical form. The function is implemented as a generator, which reduces memory usage and allows progressive output.

#### Complexity
Let n be the length of a sequence and k the kmer size. In the worst case, the time complexity is O(n . kmax + 4^k), while the space complexity can reach O(4^k) for dense kmer storage. In practice, performance mainly depends on the sequence length, its nucleotide composition, and the chosen value of kmax.

#### Conclusion
The program uses binary encoding to handle DNA words efficiently. It switches between sets and bit arrays depending on how many k-mers are present. A streaming approach is used to reduce memory usage. Canonical forms are reported to avoid duplicates. Computing reverse complements takes time proportional to the word length. The number of possible k-mers grows quickly, which limits large `kmax` values. The program runs on a single CPU core. 
To conclude, this implementation provides a clear and robust baseline solution for enumerating minimal absent words in DNA sequences. 


## Program2