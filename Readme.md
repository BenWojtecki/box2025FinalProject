# BOX Project : computing MAWs and pMAWs
by Benjamin Wojtecki and Noé Vincent


## Program 1

This program computes minimal absent word for DNA sequences. Given a set of DNA sequences in FASTA format and maximum word length `kmax`, the program enumerates all canonical minimal absent words of length `k ≤ kmax` for each sequence

#### Requirements

The following packages should be present in the python environment : 
- xopen
- argparse
- time
- array

The programs work fine for python version 3.12 and latest version of packages on 21st December 2025.

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
3. Comma-separated, lexicographically sorted list of canonical MAWs of length k.

Here, some output example 
```
plasmid_001	3	AAC,ACG,TTG
plasmid_001	4	AAGT,CCGA
```


## Program 2
This program Minimal p-Absent Words (pMAWs). Given a set of MAWS grouped by length and DNA sequences in TSV format, a proportion `p` and maximum word length `kmax`, the program enumerates all canonical minimal p-absent words of length `k ≤ kmax` for each size `k`.

#### Requirements

The following packages should be present in the python environment : 
- bitarray
- rich.progress
- pygrtie
- ahocorasick_rs
- time
- csv
- sys

The programs work fine for python version 3.12 and latest version of packages on 21st December 2025.

#### Input
The program requires a TSV file containing MAWs grouped by sequence and length. The TSV file should be formatted in the following way: 

seq0    2   AA,CC
seq0    3   ACT,TAG
seq0    5   AGTAT
seq1    2   CC,AT
seq1    6   ATCGAG
...

#### Command arguments
```bash
python3 program2.py <fastaFile> -k <kmax> -p <p> [-o output.tsv]
```

| Argument    | Description                                      |
| ----------- | ------------------------------------------------ |
| `fastaFile` | Input FASTA file with DNA sequences              |
| `-k`        | Maximum MAW length (`kmax`)          |
| `-o`        | Output TSV file (default: `resultsProgram1.tsv`) |

#### Output

The program outputs a TSV file with the following columns:

1.  Length k
2. Comma-separated lexicographically sorted list of canonical pMAWs of length.

Here, some output example :
```
3	AAC,ACG,TTG
4	AAGT,CCGA
```