# Personnal implementation of CountBF [countBF: A General-purpose High Accuracy and Space Efficient Counting Bloom Filter, Nayak, Patigiri, CNSM 2021 ]
#

import numpy as np 

class countBF():
    def __init__(self,x,y,eta):
        """
        x, eta and y should be pairwise coprime
        """
        self.bloomFilter = np.zeros(shape=(x,y), dtype = np.uint64)
        self._alpha = 64//eta
        self._x = x
        self._y = y 
        self._eta = eta 
        self._powers = None
        
        self._compute_powers()
        print(self._powers)

        self._masks = []
        for i in range(self._eta):
            f = (lambda x=x : (x//self._powers[i])%self._alpha)
            print(f"====={i}")
            print(f(5))
            self._masks.append(f)

    def _compute_powers(self):
        result = np.zeros(shape = (self._eta), dtype=np.uint32)
        result[0] = 1
        for i in range(1,self._eta):
            result[i] = (self._alpha*result[i-1])
        self._powers = result

a = countBF(1,2,8)
print(a._masks[0](5))
# f = lambda x : (x//a._powers[0])%a._alpha
# print((5//a._powers[0])%a._alpha)
# print(a._masks[0])
