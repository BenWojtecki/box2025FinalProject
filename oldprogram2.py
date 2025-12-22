from collections import Counter
import pandas as pan
from customCountBF import countBF
import prime
from math import log,sqrt,ceil,floor
import time
from rich.progress import Progress
import argparse


def naive(fileName,n, kmax, threshold):
    """
    maws given as a list
    """
    results = [set() for _ in range(kmax)]
    dataframe = pan.read_csv(fileName,delimiter="\t",names = ["k","maws"],usecols=[1,2])
    dataframe.dropna(inplace=True)
    maws = []
    with Progress() as p:
        counting = p.add_task("getting kmers", total=len(dataframe))
        for t in dataframe.itertuples():
            p.update(counting,advance = 1)
            k,l = t.k,t.maws 
            maws[0:0] = t.maws.split(",")
        print("counting")
        mawsDF = pan.DataFrame(Counter(maws).items(),columns=("maw","occ"))
        #print(mawsDF)
        print("counting done")
        data = mawsDF[mawsDF["occ"]>=threshold].sort_values(by=["occ","maw"],ascending=[False,True])
        outputing = p.add_task("outputing kmers", total = len(data["maws"]))
        for m in data["maw"]:
            p.update(outputing,advance = 1)
            results[len(m)-1].add(m)
        return results
def initCountBF(n,threshold,eta,err):

    m = ceil(-(n*log(err))/pow(log(2),2))
    k = floor((m/n)*pow(log(2),3))
    tmp = sqrt(m/256)# From KmerCo
    i = 0
    while prime.prime_numbers[i] <= tmp:
        i+=1

    return countBF(k, prime.prime_numbers[i], prime.prime_numbers[i+3],eta, threshold=threshold)

def CountBFromList(maws,n,kmax=1,threshold = 0,eta = 8,err = 0.01):
    a = initCountBF(n,threshold,eta,err)
    results = [set() for _ in range(kmax)]
    for m in maws:
        if m not in results[0]:
            res = a.checkinsert(m)
            if res[0] :
                #print(m,a.query(m))
                results[0].add(m)

def CountBFFromTSV(fileName,n, kmax, threshold = 0,eta = 8,err = 0.01):
    a =initCountBF(n,threshold,eta,err)
    falsePositives = 0
    results = [set() for _ in range(kmax)]
    dataframe = pan.read_csv(fileName,delimiter="\t",names = ["k","maws"],usecols=[1,2])
    dataframe.dropna(inplace=True)
    with Progress() as p:
        counting = p.add_task("Counting kmers", total = n)
        for t in dataframe.itertuples():
            k,l = t.k,t.maws.split(",")
            #print(l)
            for m in l:
                p.update(counting,advance = 1)
                if m not in results[k-1]:
                    b1,b2 = a.checkinsert(m)
                    if b1 : 
                        results[k-1].add(m)
                        #print("===",m,a.query(m))
                    if b2 : 
                        falsePositives +=1
    print(f"FP:{falsePositives}")
    return results

def setListToTSV(outFile,kmax,setlist):
    with open(outFile, "w") as f:
        f.write("k\tp-maws\n")
        for i  in range(kmax):
            mawSet = sorted(list(setlist[i]))
            f.write(f"{i+1}\t{','.join(mawSet)}\n")



#maws = ["a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","a","b","b","b","c","z","z","z","z"]
# maws = ["a","a","a","a","a","a","a","a","a","b","b","b","c","z","z","z","z"]

# n = len(maws)
# # print(naive(maws,2))
# a = CountBFromList(maws,n,threshold=4)
def test_BF():
    setListToTSV("resultsprogram2.tsv",10,CountBFFromTSV("resultsProgram1.tsv",204447422,kmax=10,threshold=int(2826*0.5),eta=5,err=0.001))

def test_naive():
    setListToTSV("resultsprogram2naive.tsv",10,naive("resultsProgram1.tsv",204447422,kmax=10,threshold=int(2826*0.5)))



def main():
    parser = argparse.ArgumentParser(description="Minimal Absent Words (Program 1)")
    parser.add_argument("-t", type=int, required=True)
    #parser.add_argument("-o", default="resultsProgram1.tsv")
    args = parser.parse_args()
    
    start = time.perf_counter()
    
    if args.t == 0:
        test_naive()
    else:
        test_BF()
                
    end = time.perf_counter()
    print(f"\nTotal time: {end - start:.3f} seconds")

if __name__ == "__main__":
    main()