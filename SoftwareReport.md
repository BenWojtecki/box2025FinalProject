# Software report
by Benjamin Wojtecki and Noé Vincent
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

### Main functions

The function `encodeKmersStream(seq, k)` generates all kmers from a DNA sequence and encodes them as integers using a rolling approach. Symbols that are not in the DNA alphabet (A, C, G, T) are ignored. Bit masking is used to keep the kmer length fixed, and the function operates in streaming mode to limit memory usage. Reverse complements are computed on the binary representation of kmers using the function `rcBin(x, k)`, which runs in linear time with respect to k. This approach avoids string-based operations. The function `findMawsStream` implements the main algorithm of program. For each value of k from 1 to kmax, it enumerates all kmers present in the sequence, identifies MAW. Each valid MAW is converted to its canonical form. The function is implemented as a generator, which reduces memory usage and allows progressive output.

#### Complexity
Let n be the length of a sequence and k the kmer size. In the worst case, the time complexity is $O(n . kmax + 4^k)$, while the space complexity can reach $O(4^k)$ for dense kmer storage. In practice, performance mainly depends on the sequence length, its nucleotide composition, and the chosen value of kmax.

#### Conclusion
The program uses binary encoding to handle DNA words efficiently. It switches between sets and bit arrays depending on how many k-mers are present. A streaming approach is used to reduce memory usage. Canonical forms are reported to avoid duplicates. Computing reverse complements takes time proportional to the word length. The number of possible k-mers grows quickly, which limits large `kmax` values. The program runs on a single CPU core. 
To conclude, this implementation provides a clear and robust baseline solution for enumerating minimal absent words in DNA sequences. 


## Program2

### Overview

The program computes all canonical minimal p-absent words (of lenght lower than `kmax`). The implementation prioritizes clarity, correctness, and reasonable performance for moderate values of kmax. The program is organized as a single Python file structured into four main components:

- DNA handling with function `get_reverse_complement` and `get_canonical`.
- Data extraction with `get_mawset`.
- pMAWs generation with `find_pmaws`
- I/O Handling with the `main` and `setListToTSV` functions.

### Data representation
The trie T is implemented by StringTrie from the python module pygtrie. The different Aho-Corasick automata are implemented via the ahocorasick_rs module. Both allow clarity and performance as they both rely on combat-proven fast code. (Even though pygtrie is 100% python, where ahocorasick_rs partially relies efficient on Rust code).

MAWs and word are represented using strings with characters in "ATCG".

all_pmaws containing the produced pmaws is a list of string list. Where all_pmaws[k] contains the pmaws of lenght k.

absence_map indicating from which sequences a word is absent are stored as bitarray from the python package bitarray.

### Miscelleanous about the code

Generally the code uses few specific python instruction, excepted `get_mawsets` which uses the `yield` instruction which allows to parse the TSV in a stream, this helps reducing memory consumption.

Note that `a |= b` is short for `a = a | b` where `|` represents bitwise OR.

Line 12 to 24 allow the csv package to be able to read the CSV file when its rows exceed standard size.


### Algorithm and complexity
The main function closely follow the algorithm presented in the paper (see in the repo). It's complexity is explained in the "Results" section but it achieves to solve the problem in $O(size(\mathcal{S}) + 4^{k_{max}}.N )$ where N is the number of sequences and $size(\mathcal{S})$ is the sum of the lengths of the sequences. The space complexity is in $O(size(\mathcal{S}) + 4^{k_{max}}.N )$. In reality these bounds show to be extremely pessimistic and even though the memory usage is noticeable for many value of k or p, low values of k (<=8) and p (<=0.5) seem to guarantee reasonable computation time and relevant results.

