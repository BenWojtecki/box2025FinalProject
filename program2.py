from collections import Counter
import pandas as pan
from customCountBF import countBF
import prime
from math import log,sqrt,ceil

def naive(maws, threshold):
    """
    maws given as a list
    """
    mawsDF = pan.DataFrame(Counter(maws).items(),columns=("maw","occ"))

    return mawsDF[mawsDF["occ"]>threshold].sort_values(by=["occ","maw"],ascending=[False,True])

def initCountBF(n,threshold,err):

    m = ceil(-(n*log(err))/pow(log(2),2))
    tmp = sqrt(m)//256 # From KmerCo
    k = ceil((m/n)*log(2))
    i = 0
    while prime.prime_numbers[i] <= tmp:
        i+=1

    return countBF(k, prime.prime_numbers[i], prime.prime_numbers[i+1], threshold=threshold)

def usingBF(maws,n,threshold,err):
    a = initCountBF(n,threshold,err)
    for m in maws:
        b = a.checkinsert(m)
        if b :
            print(m,a.query(m))
maws = ["a","a","a","a","b","b","b","c","z","z","z","z"]
n = len(maws)
print(naive(maws,2))
usingBF(maws,n,2,0.01)