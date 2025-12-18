# Personnal implementation of CountBF [countBF: A General-purpose High Accuracy and Space Efficient Counting Bloom Filter, Nayak, Patigiri, CNSM 2021 ]
#
import mmh3
import numpy as np 
import random as rand
from time import time

class countBF():
    def __init__(self,x,y,eta,k):
        """
        x,y is the shape of the array
        eta is the number of counters in a single cell
        x, eta and y should be pairwise coprime
        k is the number of hash function used.
        """
        rand.seed(time())
        self.bloomFilter = np.zeros(shape=(x,y), dtype = np.uint64)
        self._alpha = 64//eta
        self._x = x
        self._y = y 
        self._eta = eta 
        self._k = k

        self._powers = None
        self._compute_powers()

        self._seeds = [ int(rand.randbytes(4).hex(),16) for _ in range(self._k)]
        self._masks = lambda x,i : (x//self._powers[i])%self._powers[1]


    

    def _compute_powers(self):
        result = np.zeros(shape = (self._eta), dtype=np.uint64)
        increment = 0
        for i in range(self._eta):
            result[i] = 2**increment
            increment = self._alpha+increment
        self._powers = result



    def insert(self,kmer):
        for i in range(self._k):
            h = mmh3.hash(kmer,seed = self._seeds[i], signed = False)
            i,j,l = h%self._x,h%self._y,h%self._eta
            #numpy handles itself the eventual overflows
            self.bloomFilter[i][j] += self._powers[l]

    def query(self,kmer):
        mini_value = float("inf")
        for i in range(self._k):
            h = mmh3.hash(kmer,seed = self._seeds[i], signed = False)
            i,j,l = h%self._x,h%self._y,h%self._eta
            #numpy handles itself the eventual overflows
            read  =  self._masks(self.bloomFilter[i][j],l)
            if  read<mini_value:
                mini_value =read
        return mini_value

    def checkinsert(self,kmer,threshold):
        mini_value = float("inf")
        for i in range(self._k):
            h = mmh3.hash(kmer,seed = self._seeds[i], signed = False)
            i,j,l = h%self._x,h%self._y,h%self._eta
            #numpy handles itself the eventual overflows
            self.bloomFilter[i][j] += self._powers[l]
            read  =  self._masks(self.bloomFilter[i][j],l)
            if  read<mini_value:
                mini_value =read
        return mini_value>=threshold 
            