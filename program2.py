import pandas as pan
from bitarray import bitarray
import queue
from rich.progress import Progress
import pygtrie
from ahocorasick_rs import AhoCorasick, MatchKind

ENC = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
alphabet = "ACGT"
BASES = range(4)
complement = str.maketrans("ACGT", "TGCA")

# class TrieNode :
#     __slots__ = ["children","absMap","pAbs", "length"]
#     def __init__(self, n_seq):
#         self.children = [None] * 4
#         self.absMap = bitarray(n_seq)
#         self.absMap.setall(0)
#         self.nAbs = 0
#         self.pAbs  = False
#         self.length = 0

def get_reverse_complement(seq):
    return seq.translate(complement)[::-1]

def get_canonical(seq):
    rc = get_reverse_complement(seq)
    return seq if seq < rc else rc

def is_absent_in_seq(word, maw_set):
    """Checks if 'word' contains any MAW from the sequence's set."""
    # We check if any m in maw_set is a substring of 'word'
    for m in maw_set:
        if m in word:
            return True
    return False

def find_pmaws(maw_sets, p, max_k):
    n = len(maw_sets)
    threshold = (p * n)
    print(threshold)
    all_pmaws = [[]for _ in range(2,max_k)]
    
    print(f"Building {n} automata...")
    automata = [AhoCorasick(list(maw_set),matchkind=MatchKind.LeftmostLongest) for maw_set in maw_sets]
    # Starting with length 2 as per definition |x| > 1
    # 'active_words' stores words of length k that are NOT p-absent
    active_candidates = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"] 
    masks = [bitarray(n)for i in range(n)]
    for i in range(n):
        masks[i][i] = 1
    t = pygtrie.StringTrie()

    for k in range(2, max_k + 1):
        with Progress() as p:
            task = p.add_task(f"Checking length {k}, candidates: {len(active_candidates)}",total = len(active_candidates))
            next_candidates = []
            
            for word in active_candidates:
                p.update(task, advance=1)
                # Count how many sequences find 'word' absent 
                # (i.e., 'word' contains one of their MAWs)
                absence_map = t.longest_prefix(word).value

                if absence_map is None:
                    absence_map = bitarray(n)
                    absence_map.setall(0)
                for i in range(n):
                    # find_matches_as_indexes returns the first match found
                    if not absence_map[i]:
                        ac = automata[i]
                        subword = ac.find_matches_as_strings(word)
                        if subword :
                            sw_absence_map = t.longest_prefix(subword[0]).value
                            if sw_absence_map is not None:
                                absence_map|= sw_absence_map
                            absence_map |= masks[i]
                t[word] = absence_map
                if absence_map.count() >= threshold:
                    # It is p-absent! Since its substrings were not p-absent, 
                    # this is a Minimal p-Absent Word (pMAW).
                    all_pmaws[k-2].append(word)
                else:
                    # Not p-absent yet. Extend it to k+1 for the next round.
                    next_candidates.append(word)
            
            if not next_candidates:
                break
                
            # Generate candidates for k+1 by extending only non-p-absent words
            new_active = set()
            for word in next_candidates:
                for char in alphabet:
                    # We check canonical form here to prune the search space by half
                    # Note: This requires checking both prefix and suffix extensions

                    new_active.add(get_canonical(word + char))
                    new_active.add(get_canonical(char + word))
            
            active_candidates = list(new_active)

    return all_pmaws


# def concatStringList(a,b):


def process_data(fp,kmax):
    print("Processing Data")
    df = pan.read_csv(fp,sep='\t', names=["name","k","maws"], usecols = [0,1,2])
    df.dropna(inplace=True)
    df = df[df['k']<=kmax]
    df["maws"] = df["maws"].apply(lambda x: x.split(","))
    df = df.drop("k",axis=1)
    final = list((df.groupby("name").sum())["maws"])
    return final

print(find_pmaws(process_data("resultsProgram1.tsv",5),0.4,5))