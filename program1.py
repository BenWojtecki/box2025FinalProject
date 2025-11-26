from readfa import readfq
import numpy as np


def get_sequences(file_name):
    sequences = []
    with xopen(file_name) as fasta : 
        for name,seq,_ in readfq(fasta):
            sequences.append((name,seq,len(seq)))
    return sequences

def to_dict(mat,lx,ly,dimx,dimy):
    print(lx)
    res = {}
    cpt = 0
    label_x = list(lx.keys())
    label_y = list(ly.keys())
    for i in range(dimx):
        for j in range(dimy):
            if mat[i,j] == 1:
                res[label_x[i]+label_y[j]] = cpt
                cpt+=1
    return res,cpt

def main(sequences,kmax):
    db =[]
    alphabet = {"A":0,"T":1,"C":2,"G":3}
    epsilon = {"":0}
    dimx,dimy = (4,1)

    lx = alphabet
    ly = epsilon

    initialized = False

    for k in range (kmax):
        matrix = np.zeros((dimx,dimy))

        for i in range(len(sequences)-k):
            w =sequences[i:i+k+1]
            h = w[0]
            t = w[1:]
            matrix[lx[h],ly[t]] = 1
            
        db.append(matrix)

        print(lx,ly)
        print(matrix)
        print("===================")
        if not initialized:
            lx,dimx = to_dict(db[0],lx,ly,dimx,dimy)

        ly,dimy = to_dict(db[-1],lx,ly,dimx,dimy)
        
        
        
        initialized = True

    return ly


print(main("ATAT",4))