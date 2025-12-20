from collections import Counter
import pandas as pan
from customCountBF import countBF
import prime
from math import log,sqrt,ceil,floor
from time import perf_counter

def naive(fileName, kmax, threshold):
    """
    maws given as a list
    """
    results = [set() for _ in range(kmax)]
    dataframe = pan.read_csv(fileName,delimiter="\t",names = ["k","maws"],usecols=[1,2])
    dataframe.dropna(inplace=True)
    maws = []
    for t in dataframe.itertuples():
        k,l = t.k,t.maws 
        maws[0:0] = t.maws.split(",")
    mawsDF = pan.DataFrame(Counter(maws).items(),columns=("maw","occ"))
    #print(mawsDF)
    data = mawsDF[mawsDF["occ"]>=threshold].sort_values(by=["occ","maw"],ascending=[False,True])
    for m in data["maw"]:
        results[len(m)-1].add(m)
    return results
def initCountBF(n,threshold,err):

    m = ceil(-(n*log(err))/pow(log(2),2))
    k = floor((m/n)*pow(log(2),3))
    tmp = sqrt(m/256)# From KmerCo
    i = 0
    while prime.prime_numbers[i] <= tmp:
        i+=1

    return countBF(k, prime.prime_numbers[i], prime.prime_numbers[i+3], threshold=threshold)

def CountBFromList(maws,n,kmax=1,threshold = 0,err = 0.01):
    a = initCountBF(n,threshold,err)
    results = [set() for _ in range(kmax)]
    for m in maws:
        if m not in results[0]:
            res = a.checkinsert(m)
            if res[0] :
                #print(m,a.query(m))
                results[0].add(m)

def CountBFFromTSV(fileName,n, kmax, threshold = 0,err = 0.01):
    a =initCountBF(n,threshold,err)
    falsePositives = 0
    results = [set() for _ in range(kmax)]
    dataframe = pan.read_csv(fileName,delimiter="\t",names = ["k","maws"],usecols=[1,2])
    dataframe.dropna(inplace=True)
    for t in dataframe.itertuples():
        k,l = t.k,t.maws.split(",")
        #print(l)
        for m in l:
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

tmp3 = perf_counter()
setListToTSV("resultsprogram2.tsv",2,CountBFFromTSV("resultsProgram1.tsv",449,kmax=2,threshold=25,err=0.001))
tmp4 = perf_counter()

print(tmp4-tmp3)

tmp = perf_counter()
setListToTSV("resultsprogram2naive.tsv",2,naive("resultsProgram1.tsv",kmax=2,threshold=25))
tmp2 = perf_counter()

print(tmp2-tmp)


