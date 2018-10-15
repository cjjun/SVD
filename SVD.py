import numpy as np 
import math


def Diag(a):
    row,col = np.shape(a)
    U = np.zeros_like(a,dtype = 'float64')
    S = np.zeros((1,row),dtype = 'float64')

    iszero = 0
    for i in range(row):
        x =  np.random.randn(row,1)

        for j in range(100):
            for k in range(i):
                x -= np.reshape( U[:,k] * U[:,k].transpose().dot(x) , (row,1) )

            x = np.dot(a,x)
            if np.linalg.norm(x) < 1e-12:
                a = a + np.eye((row))
                iszero = 1
                x = np.random.randn(row,1)
            x /= np.linalg.norm(x)
        if not iszero:
            S[0,i] = x.transpose().dot(a).dot(x)
        else:
            S[0,i] = 0

        U[:,i] = x[:,0]
    return S,U

def SVD(A):
    row,col = np.shape(A)
    MU = A.dot(A.transpose())
    MV = A.transpose().dot(A)

    s,U = Diag(MU)
    s,V = Diag(MV)

    S = np.zeros_like(A,dtype = 'float64')
    for i in range( min(col,row) ):
            S[i,i] = math.sqrt(s[0,i])

    av = A.dot(V)
   
    for i in range( min(col,row) ):
        if abs(S[i,i]) > 1e-12:
            U[:,i] = av[:,i] / S[i,i]
    return U,S,V

if __name__ == '__main__':
    a = np.array([[0,0.5,0.25],[0.5,0.25,0.125],[0.25,0.125,0.5*0.125]]) 
    # a = a.transpose().dot(a)
    [u,s,v] = SVD(a)
    [u,s,v] = np.linalg.svd(a)
    print(u)
    print(s)
    print(v)
    # print(u.dot(s).dot(v.transpose()))
