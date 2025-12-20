Author1 : Noé Vincent
Author2 : Benjamin Wojtecki
### Final Rapport 

## Program 1

This program comptures minimal absent word for DNA sequences. Given a set of DNA sequences in FASTA format and maximum word length `kmax`, the program enumarates all canonical minimal absent words of length `k ≤ kmax` for each sequence 

#### Input
The program requires FASTA file containing one or more DNA sequences and a maximum length `kmax` for MAWs. 

#### Command arguments
```bash
python program1.py <fastaFile> -k <kmax> [-o output.tsv]
```

| Argument    | Description                                      |
| ----------- | ------------------------------------------------ |
| `fastaFile` | Input FASTA file with DNA sequences              |
| `-k`        | Maximum MAW length (`kmax`)          |
| `-o`        | Output TSV file (default: `resultsProgram1.tsv`) |

#### Output

The program outputs a TSV file with the following columns:

1. Sequence name
2. Length k
3. Comma-separated list of canonical MAWs of length k with lexicographically sorted

Here, some output example 
```
plasmid_001	3	AAC,ACG,TTG
plasmid_001	4	AAGT,CCGA
```


## Program 2
Let's use bloom filters, in particular a variant presented here : 
Personnal implementation of CountBF [countBF: A General-purpose High Accuracy and Space Efficient Counting Bloom Filter, Nayak, Patigiri, CNSM 2021 ]
https://arxiv.org/pdf/2106.04364
bl