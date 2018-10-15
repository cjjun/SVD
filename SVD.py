import numpy as np 
import math

# Diagonize a real symmetric matrix(n*n), return eigenvalues vector(1*n) and orthogonal matrix(n*n)
def Diag(a):                                        
    row,col = np.shape(a)
    U = np.zeros_like(a,dtype = 'float64')
    S = np.zeros((1,row),dtype = 'float64')
    iter = 100
    iszero = 0
    for i in range(row):                        # Use the iteration x_{n+1} = ||Ax_{n}|| , x = x / ||x|| to find x with respect to eigenvalue with largest norm 
        x =  np.random.randn(row,1)

        for j in range(iter):                   # Gram-Smith normalization.
            for k in range(i):
                x -= np.reshape( U[:,k] * U[:,k].transpose().dot(x) , (row,1) )

            x = np.dot(a,x)
            if np.linalg.norm(x) < 1e-12:       # The remaining eigenvalues are zeros, use B = A + I to undertake iteration
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

# The svd implementation with A(m*n)
def SVD(A):                             
    row,col = np.shape(A)
    MU = A.dot(A.transpose())           # Seperately calculate left/right eigenvector matrix
    MV = A.transpose().dot(A)

    s,U = Diag(MU)
    s,V = Diag(MV)

    S = np.zeros_like(A,dtype = 'float64')  # S is the sigma matrix(m*n) with singular values diagonally. 
    for i in range( min(col,row) ):
            S[i,i] = math.sqrt(s[0,i])

    av = A.dot(V)                           # The left and right eigenvectors might not be paired, with a '-' difference. AV = US, we use AV and S to adjust the minus signal in U
   
    for i in range( min(col,row) ):
        if abs(S[i,i]) > 1e-12:             # If the singular values is zero, then no necessity to adjust U. But calculating U is necessary since singular value zeros can't help getting U(:,i)
            U[:,i] = av[:,i] / S[i,i]
    return U,S,V

if __name__ == '__main__':
    a = np.array([[0,0.5,0.25],[0.5,0.25,0.125],[0.25,0.125,0.5*0.125]]) 
    # a = a.transpose().dot(a)
    [u,s,v] = SVD(a)
    print(u)
    print(s)
    print(v)
    # print(u.dot(s).dot(v.transpose()))
